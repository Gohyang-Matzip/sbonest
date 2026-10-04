"""Optional process pool for Sideband work; serial by default, numerically identical.

Workers are started with the spawn method so every platform behaves alike. Each
worker builds one SidebandModel from the configuration and keeps its prepared
data cached. Tasks are plain functions that reuse the serial code paths, so a
row produced in a worker is byte-identical to the row the main process would
produce. Mapping preserves submission order, which keeps checkpoint records an
ordered prefix even when later items finish first.
"""
import concurrent.futures
import multiprocessing
import os
import time

import numpy as np

_WORKER = {}


def validate_workers(value):
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError('workers must be a positive integer')
    return value


def _init_worker(config, config_dir):
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
        os.environ.setdefault(name, '1')
    os.environ.setdefault('MPLBACKEND', 'Agg')
    from sbfit import SidebandModel

    model = SidebandModel(config, config_dir)
    model.verbose = False
    model._fit_data = model._prepare_data()
    _WORKER['model'] = model


def _model():
    try:
        return _WORKER['model']
    except KeyError:
        raise RuntimeError('Worker pool task ran outside an initialized worker') from None


def residual_task(p_full):
    """Full residual vector at a full parameter vector; nonfinite values propagate."""
    model = _model()
    if model._fit_data is None:
        model._fit_data = model._prepare_data()
    values = np.asarray(model.errFunc(np.asarray(p_full, dtype=float)), dtype=float)
    if not np.isfinite(values).all():
        raise ValueError('Non-finite Sideband residual')
    return values


def warm_task(_):
    """Force a worker to finish its initializer; returns its process id."""
    _model()
    return os.getpid()


def predict_task(payload):
    """Model predictions for the data blocks assigned to this task."""
    model = _model()
    if model._fit_data is None:
        model._fit_data = model._prepare_data()
    P = model.seParam(np.asarray(payload['p'], dtype=float))
    data = model._fit_data
    return [np.asarray(model.calc(P, data[k][0], data[k][2], data[k][1]), dtype=float)
            for k in payload['blocks']]


def balanced_assignment(sizes, tasks):
    """Assign block indices to at most ``tasks`` groups with balanced point counts."""
    tasks = max(1, min(int(tasks), len(sizes)))
    groups = [[] for _ in range(tasks)]
    loads = [0] * tasks
    for index in sorted(range(len(sizes)), key=lambda k: -sizes[k]):
        target = loads.index(min(loads))
        groups[target].append(index)
        loads[target] += sizes[index]
    return [sorted(group) for group in groups if group]


def multistart_task(payload):
    from sb_analysis import fit_attempt

    model = _model()
    return fit_attempt(model, np.asarray(payload['start'], dtype=float), payload['source'],
                       payload['index'], payload['cfg'])


def profile_task(payload):
    from sb_analysis import profile_point

    model = _model()
    p = np.asarray(payload['p'], dtype=float)
    # Publish the base fit's bounds and free indices exactly as restore_fit does.
    cfg = {key: value for key, value in payload['cfg'].items() if key != 'initial'}
    model.prepare_fit(p0=p, fitting_config=cfg)
    model._fit_data = model._prepare_data()
    try:
        return profile_point(model, p, payload['name'], float(payload['target']),
                             float(payload['baseline']))
    finally:
        model.errFunc(p)


def bootstrap_task(payload):
    from sb_bootstrap import bootstrap_replicate

    model = _model()
    observed = [np.asarray(values, dtype=float) for values in payload['observed']]
    return bootstrap_replicate(model, np.asarray(payload['p'], dtype=float), payload['cfg'],
                               payload['index'], observed)


def study_task(payload):
    """One coverage-study replicate (draws, refit, optional intervals and inner bootstrap)."""
    from validate_uncertainty import study_replicate

    model = _model()
    return study_replicate(model, np.asarray(payload['p'], dtype=float), payload['cfg'], payload['settings'],
                           payload['index'], interval_settings=payload['interval_settings'],
                           inner_settings=payload['inner_settings'])


class WorkerPool:
    """Process pool bound to one configuration; use as a context manager."""

    def __init__(self, workers, config, config_dir='.'):
        self.workers = validate_workers(workers)
        self.config, self.config_dir = config, str(config_dir)
        self.executor = None

    def __enter__(self):
        self.executor = concurrent.futures.ProcessPoolExecutor(
            max_workers=self.workers, mp_context=multiprocessing.get_context('spawn'),
            initializer=_init_worker, initargs=(self.config, self.config_dir))
        # Start every worker now: the executor would otherwise spawn them one by one
        # on demand, which made the first Jacobian several times slower.
        futures = [self.executor.submit(warm_task, index) for index in range(self.workers)]
        try:
            self.started = sorted({future.result() for future in futures})
        except BaseException:
            self.__exit__(None, None, None)
            raise
        self._assignments = {}
        return self

    def __exit__(self, *_):
        if self.executor is not None:
            self.executor.shutdown(wait=True, cancel_futures=True)
            self.executor = None

    def map(self, task, items):
        """Yield results in submission order; the first failure propagates."""
        if self.executor is None:
            raise RuntimeError('Worker pool is not active')
        return self.executor.map(task, list(items))

    def evaluate_many(self, vectors):
        return list(self.map(residual_task, [np.asarray(v, dtype=float) for v in vectors]))

    def predict_blocks(self, p, sizes):
        """Predictions for every data block of the worker models at full vector p.

        Blocks are grouped so each worker receives a similar number of points;
        the result lists one array per block in data order.
        """
        key = tuple(sizes)
        if key not in self._assignments:
            self._assignments[key] = balanced_assignment(list(sizes), self.workers)
        groups = self._assignments[key]
        p = np.asarray(p, dtype=float)
        results = self.map(predict_task, [{'p': p, 'blocks': group} for group in groups])
        predictions = [None] * len(sizes)
        for group, values in zip(groups, results):
            for index, array in zip(group, values):
                predictions[index] = array
        return predictions


class Progress:
    """Print elapsed time and a remaining-time estimate for a counted loop."""

    def __init__(self, label, total, *, enabled=True, done=0):
        self.label, self.total, self.enabled = label, int(total), enabled
        self.done, self.started = int(done), time.monotonic()
        self.new = 0

    def step(self):
        self.done += 1
        self.new += 1
        if not self.enabled:
            return
        elapsed = time.monotonic() - self.started
        per_item = elapsed / self.new
        remaining = per_item * (self.total - self.done)
        print(f'[{self.label}] {self.done}/{self.total} done, {per_item:.1f} s/item, '
              f'elapsed {elapsed:.0f} s, ETA {remaining:.0f} s', flush=True)
