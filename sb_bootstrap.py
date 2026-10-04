"""Seeded parametric bootstrap conditional on the fitted model and absolute sigma."""
import copy
import numbers

import numpy as np

from sb_analysis import _restore, _state
from sb_report import derived_errors


def validate_bootstrap(settings):
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


def bootstrap_fit(model, p, settings, *, completed=None, on_complete=None):
    """Refit independent Gaussian draws, retaining every outcome and original state.

    SeedSequence([seed, index]) makes each draw independent of resume position.
    Completion callbacks run outside the optimizer exception handler so failed
    persistence and user interrupts cannot be mistaken for failed fits.
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
    old_verbose = model.verbose
    data = model._prepare_data()
    original = [(es, es.int) for _, es, *_ in data]
    free = set(int(i) for i in model.free)
    varied = [name for i, name in enumerate(names) if i in free]
    cfg = {key: value for key, value in model.config['init'].items()
           if key not in ('initial', 'multistart', 'profile', 'bootstrap')}
    parameters = model.seParam(p)
    means = [np.asarray(model.calc(parameters, i, offsets, es), dtype=float)
             for i, es, offsets, _, _ in data]
    if any(not np.isfinite(values).all() for values in means):
        raise ValueError('Nonfinite bootstrap generating predictions')
    try:
        model.verbose = False
        for index in range(len(samples), settings['replicates']):
            rng = np.random.default_rng(np.random.SeedSequence([settings['seed'], index]))
            for (_, es, _, _, sigma), mean in zip(data, means):
                es.int = (mean + rng.normal(size=mean.shape) * sigma).tolist()
            model._fit_data = None
            row = {'index': index, 'success': False, 'parameters': None,
                   'stderr': None, 'derived_se': None, 'kex': None, 'pB': None,
                   'chi2': None, 'message': '', 'at_bounds': []}
            try:
                fitted, covariance = model.fit(p0=p.copy(), fitting_config=cfg)
                if not model.result.success or not np.isfinite(fitted).all():
                    raise RuntimeError('Bootstrap fit did not converge to finite parameters')
                residual = np.asarray(model.errFunc(fitted), dtype=float)
                if not np.isfinite(residual).all():
                    raise RuntimeError('Nonfinite bootstrap residual')
                total = fitted[names.index('kab')] + fitted[names.index('kba')]
                if total <= 0:
                    raise ValueError('Invalid bootstrap exchange rates')
                rate_indices = [names.index('kab'), names.index('kba')]
                errors = derived_errors(fitted[rate_indices], covariance[np.ix_(rate_indices, rate_indices)])
                row.update(success=True, parameters=fitted.tolist(),
                           stderr=[float(np.sqrt(value)) if np.isfinite(value) and value >= 0 else None
                                   for value in np.diag(covariance)],
                           derived_se=errors, kex=float(total),
                           pB=float(fitted[names.index('kab')]/total),
                           chi2=float(residual @ residual), message=str(model.result.message),
                           at_bounds=[name for i, name in enumerate(names) if i in free and (
                               np.isclose(fitted[i], model.lower[i], rtol=1e-5, atol=1e-8)
                               or np.isclose(fitted[i], model.upper[i], rtol=1e-5, atol=1e-8))])
            except Exception as exc:
                row['message'] = f'{type(exc).__name__}: {exc}'
            samples.append(row)
            if on_complete is not None:
                on_complete(copy.deepcopy(row))
    finally:
        for es, values in original:
            es.int = values
        model._fit_data = old_cache
        model.verbose = old_verbose
        _restore(model, state)
    good = [row for row in samples if row['success']]
    alpha = (1-settings['confidence'])/2
    intervals = {}
    for name in [*names, 'kex', 'pB']:
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
