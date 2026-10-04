"""Seeded parametric bootstrap conditional on the fitted model and absolute sigma."""
import copy
import numbers

import numpy as np

from sb_analysis import _restore, _state
from sb_report import derived_errors


def validate_bootstrap(settings):
    """Validate init.bootstrap settings; returns normalized replicates, seed and confidence."""
    if not isinstance(settings, dict) or set(settings) - {'replicates', 'seed', 'confidence'}:
        raise ValueError('bootstrap allows replicates, seed and confidence')
    count, seed = settings.get('replicates'), settings.get('seed')
    confidence = settings.get('confidence', .95)
    if isinstance(count, bool) or not isinstance(count, int) or count < 1:
        raise ValueError('bootstrap.replicates must be a positive integer')
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError('bootstrap.seed must be an explicit nonnegative integer')
    if (isinstance(confidence, bool) or not isinstance(confidence, numbers.Real)
            or not np.isfinite(confidence) or not 0 < confidence < 1):
        raise ValueError('bootstrap.confidence must be between 0 and 1')
    return {'replicates': count, 'seed': seed, 'confidence': float(confidence)}


def bootstrap_draws(model, p, settings, index):
    """Independent Gaussian observations for replicate ``index``, in data order.

    SeedSequence([seed, index]) makes each draw independent of resume position
    and of the process that performs the fit.
    """
    data = model._fit_data if getattr(model, '_fit_data', None) is not None else model._prepare_data()
    parameters = model.seParam(np.asarray(p, dtype=float))
    means = [np.asarray(model.calc(parameters, i, offsets, es), dtype=float)
             for i, es, offsets, _, _ in data]
    if any(not np.isfinite(values).all() for values in means):
        raise ValueError('Nonfinite bootstrap generating predictions')
    rng = np.random.default_rng(np.random.SeedSequence([settings['seed'], index]))
    return [mean + rng.normal(size=mean.shape) * sigma
            for (_, _, _, _, sigma), mean in zip(data, means)]


def sample_row(model, fitted, covariance, index):
    """JSON-safe summary of one refit: parameters, errors, kex/pB, chi2, bounds."""
    names = list(model.parameter_names)
    fitted = np.asarray(fitted, dtype=float)
    if not model.result.success or not np.isfinite(fitted).all():
        raise RuntimeError('Bootstrap fit did not converge to finite parameters')
    residual = np.asarray(model.errFunc(fitted), dtype=float)
    if not np.isfinite(residual).all():
        raise RuntimeError('Nonfinite bootstrap residual')
    total = fitted[names.index('kab')] + fitted[names.index('kba')]
    if total <= 0:
        raise ValueError('Invalid bootstrap exchange rates')
    two_state = 'kbc' not in names
    rate_indices = [names.index('kab'), names.index('kba')]
    errors = (derived_errors(fitted[rate_indices], covariance[np.ix_(rate_indices, rate_indices)])
              if two_state else {'kex': None, 'pB': None})
    free = set(int(i) for i in model.free)
    return {'index': int(index), 'success': True, 'parameters': fitted.tolist(),
            'stderr': [float(np.sqrt(value)) if np.isfinite(value) and value >= 0 else None
                       for value in np.diag(covariance)],
            'derived_se': errors, 'kex': float(total) if two_state else None,
            'pB': float(fitted[names.index('kab')]/total) if two_state else None,
            'chi2': float(residual @ residual), 'message': str(model.result.message),
            'at_bounds': [name for i, name in enumerate(names) if i in free and (
                np.isclose(fitted[i], model.lower[i], rtol=1e-5, atol=1e-8)
                or np.isclose(fitted[i], model.upper[i], rtol=1e-5, atol=1e-8))]}


def bootstrap_replicate(model, p, cfg, index, observed, *, analyses=None):
    """Fit one replicate's observations from the baseline; restore the model's data.

    Returns the JSON-safe sample row. Shared by the serial loop and pool workers.
    ``analyses(model, fitted, covariance)`` may run further analyses on the
    replicate's data while it is in place and return a dict merged into the row.
    """
    p = np.asarray(p, dtype=float)
    data = model._prepare_data()
    if len(observed) != len(data) or any(np.shape(values) != np.shape(es.int) for values, (_, es, *_) in zip(observed, data)):
        raise ValueError('Bootstrap observations do not match the data layout')
    original = [(es, es.int) for _, es, *_ in data]
    old_cache, old_verbose = getattr(model, '_fit_data', None), model.verbose
    row = {'index': int(index), 'success': False, 'parameters': None,
           'stderr': None, 'derived_se': None, 'kex': None, 'pB': None,
           'chi2': None, 'message': '', 'at_bounds': []}
    try:
        model.verbose = False
        for (es, _), values in zip(original, observed):
            es.int = np.asarray(values, dtype=float).tolist()
        model._fit_data = None
        try:
            fitted, covariance = model.fit(p0=p.copy(), fitting_config=cfg)
            row.update(sample_row(model, fitted, covariance, index))
            if analyses is not None:
                row.update(analyses(model, fitted, covariance))
        except Exception as exc:
            row['message'] = f'{type(exc).__name__}: {exc}'
    finally:
        for es, values in original:
            es.int = values
        model._fit_data = old_cache
        model.verbose = old_verbose
    return row


def bootstrap_fit(model, p, settings, *, completed=None, on_complete=None, pool=None):
    """Refit independent Gaussian draws, retaining every outcome and original state.

    SeedSequence([seed, index]) makes each draw independent of resume position.
    Completion callbacks run outside the optimizer exception handler so failed
    persistence and user interrupts cannot be mistaken for failed fits.
    With ``pool``, replicates are fitted in worker processes in index order from
    observations drawn in this process.
    """
    settings = validate_bootstrap(settings)
    p = np.asarray(p, dtype=float)
    names = list(model.parameter_names)
    if p.shape != (len(names),) or not np.isfinite(p).all():
        raise ValueError('Invalid bootstrap baseline parameter vector')
    samples = copy.deepcopy([] if completed is None else completed)
    if not isinstance(samples, list) or len(samples) > settings['replicates']:
        raise ValueError('Invalid completed bootstrap records')
    for index, row in enumerate(samples):
        if (not isinstance(row, dict) or row.get('index') != index
                or not isinstance(row.get('success'), bool)):
            raise ValueError('Completed bootstrap records must be an ordered prefix')
        if row['success']:
            values = np.asarray(row.get('parameters'), dtype=float)
            if values.shape != p.shape or not np.isfinite(values).all():
                raise ValueError('Invalid completed bootstrap parameters')
    state = _state(model)
    old_cache = getattr(model, '_fit_data', None)
    free = set(int(i) for i in model.free)
    varied = [name for i, name in enumerate(names) if i in free]
    init = model.init_config() if hasattr(model, 'init_config') else model.config['init']
    cfg = {key: value for key, value in init.items()
           if key not in ('initial', 'multistart', 'profile', 'bootstrap', 'profile_interval')}
    pending = list(range(len(samples), settings['replicates']))
    progress = None
    if pending:
        from sb_parallel import Progress

        progress = Progress('bootstrap', settings['replicates'],
                            enabled=bool(getattr(model, 'verbose', False)), done=len(samples))

    def accept(row):
        samples.append(row)
        if on_complete is not None:
            on_complete(copy.deepcopy(row))
        if progress is not None:
            progress.step()

    try:
        if pending and pool is not None:
            from sb_parallel import bootstrap_task

            payloads = [{'p': p.tolist(), 'cfg': cfg, 'index': index,
                         'observed': [values.tolist() for values in bootstrap_draws(model, p, settings, index)]}
                        for index in pending]
            for row in pool.map(bootstrap_task, payloads):
                accept(row)
        else:
            for index in pending:
                observed = bootstrap_draws(model, p, settings, index)
                accept(bootstrap_replicate(model, p, cfg, index, observed))
    finally:
        model._fit_data = old_cache
        _restore(model, state)
    good = [row for row in samples if row['success']]
    alpha = (1-settings['confidence'])/2
    intervals = {}
    derived_names = ['kex', 'pB'] if 'kbc' not in names else []
    for name in [*names, *derived_names]:
        index = names.index(name) if name in names else None
        values = [row['parameters'][index] if index is not None else row[name] for row in good]
        fixed = index not in free if index is not None else not any(names[i] in ('kab', 'kba') for i in free)
        if name == 'pB':
            fixed = fixed or any(names.index(rate) not in free and p[names.index(rate)] == 0
                                 for rate in ('kab', 'kba'))
        bounds = np.quantile(values, [alpha, .5, 1-alpha]).tolist() if values else [None]*3
        intervals[name] = dict(zip(('lower', 'median', 'upper'), bounds))
        intervals[name].update(n_success=len(good), fixed=fixed)
    warnings = []
    if len(good) < 100:
        warnings.append('Fewer than 100 successful replicates: percentile intervals are exploratory and tail estimates are unstable.')
    if len(good) != len(samples):
        warnings.append('Intervals use successful replicates only; fitting failures can bias the distribution.')
    if any(row['at_bounds'] for row in good):
        warnings.append('Some replicates reach bounds; inspect boundary frequencies and identifiability.')
    return {**settings, 'successful': len(good), 'parameter_names': names,
            'varying_parameters': varied, 'samples': samples, 'intervals': intervals,
            'warnings': warnings,
            'interpretation': 'Parametric Gaussian percentile intervals conditional on the generating model, supplied absolute sigma and fixed inputs; excludes model mismatch. Fixed intervals are assumptions.'}
