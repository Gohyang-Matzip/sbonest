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
