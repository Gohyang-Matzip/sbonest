"""Seeded Metropolis sampling of the fitted model under the supplied absolute sigma.

``init.mcmc`` draws from exp(-chi2/2) with flat priors inside the fit's bounds,
so the intervals are posterior intervals for the selected model and fixed
inputs. Each chain is a random walk whose proposal is the local covariance of
the free parameters, scaled during burn-in toward a 23% acceptance rate and
frozen afterwards. Chains are independent functions of (seed, chain index), so
serial, pooled and resumed runs give identical samples. Agreement between
chains (split R-hat) and the effective sample size are reported; they detect
poor mixing, not model error.
"""
import copy
import numbers

import numpy as np

from sb_analysis import _restore, _state

KEEP = 1000  # at most this many retained draws per chain


def validate_mcmc(settings):
    """Validate init.mcmc; returns normalized steps, seed, chains, burn_in and confidence."""
    if not isinstance(settings, dict) or set(settings) - {'steps', 'seed', 'chains', 'burn_in', 'confidence'}:
        raise ValueError('mcmc allows steps, seed, chains, burn_in and confidence')

    def integer(name, default, minimum):
        value = settings.get(name, default)
        if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
            raise ValueError(f'mcmc.{name} must be an integer of at least {minimum}')
        return value

    steps = integer('steps', None, 100)
    seed = integer('seed', None, 0)
    chains = integer('chains', 4, 2)
    burn = integer('burn_in', steps // 4, 0)
    if burn >= steps:
        raise ValueError('mcmc.burn_in must be smaller than mcmc.steps')
    confidence = settings.get('confidence', .95)
    if (isinstance(confidence, bool) or not isinstance(confidence, numbers.Real)
            or not np.isfinite(confidence) or not 0 < confidence < 1):
        raise ValueError('mcmc.confidence must be between 0 and 1')
    return {'steps': steps, 'seed': seed, 'chains': chains, 'burn_in': burn, 'confidence': float(confidence)}


def _proposal_factor(model, covariance):
    free = np.asarray(model.free, dtype=int)
    local = np.asarray(covariance, dtype=float)[np.ix_(free, free)]
    try:
        if not np.isfinite(local).all():
            raise np.linalg.LinAlgError
        return np.linalg.cholesky(local)
    except np.linalg.LinAlgError:
        raise ValueError('mcmc needs a finite positive-definite local covariance of the free parameters; '
                         'fix weakly determined parameters or those at bounds') from None


def mcmc_chain(model, p, covariance, settings, index):
    """Run chain ``index`` from the fit; returns its JSON-safe row with thinned post-burn-in draws.

    ``samples`` holds free-parameter vectors in the row's ``parameter_names`` order. Shared by the
    serial loop and pool workers.
    """
    p = np.asarray(p, dtype=float)
    free = np.asarray(model.free, dtype=int)
    lower, upper = np.asarray(model.lower)[free], np.asarray(model.upper)[free]
    factor = _proposal_factor(model, covariance)
    rng = np.random.default_rng(np.random.SeedSequence([settings['seed'], index]))
    full = p.copy()

    def half_chi2(x):
        if np.any(x < lower) or np.any(x > upper):
            return np.inf
        full[free] = x
        residual = np.asarray(model.errFunc(full), dtype=float)
        value = .5 * float(residual @ residual)
        return value if np.isfinite(value) else np.inf

    state = _state(model)
    try:
        x = p[free].copy()
        energy = half_chi2(x)
        if not np.isfinite(energy):
            raise ValueError('mcmc cannot start from a non-finite chi2')
        # Disperse the starts by one local-covariance draw; keep the fit if it leaves the bounds.
        start = x + factor @ rng.normal(size=len(free))
        start_energy = half_chi2(start)
        if np.isfinite(start_energy):
            x, energy = start, start_energy
        scale = 2.38 / np.sqrt(len(free))
        steps, burn = settings['steps'], settings['burn_in']
        thin = max(1, (steps - burn) // KEEP)
        samples, accepted, window = [], 0, 0
        for step in range(steps):
            candidate = x + scale * (factor @ rng.normal(size=len(free)))
            candidate_energy = half_chi2(candidate)
            if np.log(rng.uniform()) < energy - candidate_energy:
                x, energy = candidate, candidate_energy
                window += 1
                accepted += step >= burn
            if step < burn and (step + 1) % 50 == 0:
                scale *= np.exp(window / 50 - .234)
                window = 0
            if step >= burn and (step - burn) % thin == 0:
                samples.append(x.tolist())
    finally:
        _restore(model, state)
    return {'index': int(index), 'accept_rate': accepted / (steps - burn), 'proposal_scale': float(scale),
            'thin': thin, 'parameter_names': [model.parameter_names[i] for i in free], 'samples': samples}


def _split_rhat(chains):
    """Split R-hat of an array (chains, draws); None when it is undefined."""
    half = chains.shape[1] // 2
    if half < 2:
        return None
    parts = np.concatenate([chains[:, :half], chains[:, half:2 * half]])
    within = parts.var(axis=1, ddof=1).mean()
    if within <= 0:
        return None
    between = half * parts.mean(axis=1).var(ddof=1)
    return float(np.sqrt(((half - 1) / half * within + between / half) / within))


def _effective_size(chains):
    """Effective sample size from the chain-averaged autocorrelation (initial positive sequence)."""
    count, draws = chains.shape
    centered = chains - chains.mean(axis=1, keepdims=True)
    variance = float(np.mean(centered**2))
    if variance <= 0 or draws < 4:
        return None
    total = 0.
    for lag in range(1, draws - 1, 2):
        pair = sum(float(np.mean(centered[:, k:] * centered[:, :draws - k])) * (draws - k) / draws
                   for k in (lag, lag + 1)) / variance
        if pair <= 0:
            break
        total += pair
    return float(count * draws / (1 + 2 * total))


def mcmc_fit(model, p, covariance, settings, *, completed=None, on_complete=None, pool=None):
    """Sample every chain (resuming an ordered prefix) and summarize the pooled draws.

    The returned summary holds per-parameter posterior quantiles, split R-hat and
    effective sample size, plus kex and pB for two-state models; the draws stay in
    the rows passed to ``on_complete`` (the checkpoint), not in the summary.
    """
    settings = validate_mcmc(settings)
    p = np.asarray(p, dtype=float)
    names = list(model.parameter_names)
    free = [int(i) for i in model.free]
    _proposal_factor(model, covariance)
    rows = copy.deepcopy([] if completed is None else completed)
    if not isinstance(rows, list) or len(rows) > settings['chains']:
        raise ValueError('Invalid completed mcmc records')
    for index, row in enumerate(rows):
        draws = np.asarray(row.get('samples') if isinstance(row, dict) else None, dtype=float)
        if (not isinstance(row, dict) or row.get('index') != index or draws.ndim != 2
                or draws.shape[1] != len(free) or not np.isfinite(draws).all()):
            raise ValueError('Completed mcmc records must be an ordered prefix of valid chains')
    pending = list(range(len(rows), settings['chains']))
    progress = None
    if pending:
        from sb_parallel import Progress

        progress = Progress('mcmc', settings['chains'], enabled=bool(getattr(model, 'verbose', False)), done=len(rows))

    def accept(row):
        rows.append(row)
        if on_complete is not None:
            on_complete(copy.deepcopy(row))
        progress.step()

    if pending and pool is not None:
        from sb_parallel import mcmc_task

        init = model.init_config() if hasattr(model, 'init_config') else model.config['init']
        cfg = {key: value for key, value in init.items()
               if key not in ('initial', 'multistart', 'profile', 'bootstrap', 'profile_interval', 'mcmc')}
        payloads = [{'p': p.tolist(), 'covariance': np.asarray(covariance, dtype=float).tolist(), 'cfg': cfg,
                     'settings': settings, 'index': index} for index in pending]
        for row in pool.map(mcmc_task, payloads):
            accept(row)
    else:
        for index in pending:
            accept(mcmc_chain(model, p, covariance, settings, index))
    draws = np.asarray([row['samples'] for row in rows], dtype=float)  # chains, draws, free
    alpha = (1 - settings['confidence']) / 2
    columns = {names[i]: draws[:, :, k] for k, i in enumerate(free)}
    if 'kbc' not in names and {'kab', 'kba'} & set(columns):
        kab, kba = (columns.get(name, np.full(draws.shape[:2], p[names.index(name)])) for name in ('kab', 'kba'))
        columns.update(kex=kab + kba, pB=kab / (kab + kba))
    parameters = {}
    for name, values in columns.items():
        lo, median, hi = np.quantile(values, [alpha, .5, 1 - alpha]).tolist()
        parameters[name] = {'lower': lo, 'median': median, 'upper': hi, 'mean': float(values.mean()),
                            'sd': float(values.std(ddof=1)), 'rhat': _split_rhat(values),
                            'ess': _effective_size(values)}
    rates = [row['accept_rate'] for row in rows]
    warnings = []
    mixing = [name for name, item in parameters.items() if item['rhat'] is None or item['rhat'] > 1.05]
    if mixing:
        warnings.append('MCMC chains disagree (split R-hat above 1.05 or undefined) for: ' + ', '.join(mixing)
                        + '; run more steps before using these intervals.')
    thin_ess = [name for name, item in parameters.items() if item['ess'] is None or item['ess'] < 100]
    if thin_ess:
        warnings.append('MCMC effective sample size below 100 for: ' + ', '.join(thin_ess)
                        + '; interval limits are unstable.')
    if min(rates) < .05 or max(rates) > .7:
        warnings.append(f'MCMC acceptance rates {min(rates):.2g} to {max(rates):.2g} are outside 0.05 to 0.7; '
                        'the local covariance is a poor proposal here.')
    return {**settings, 'parameter_names': list(columns), 'draws_per_chain': int(draws.shape[1]),
            'accept_rates': rates, 'parameters': parameters, 'warnings': warnings,
            'interpretation': 'Posterior quantiles for exp(-chi2/2) with flat priors inside the fit bounds, '
                              'conditional on the selected model, supplied absolute sigma and fixed inputs; '
                              'excludes model mismatch. Draws are kept in the checkpoint records mcmc-N.'}
