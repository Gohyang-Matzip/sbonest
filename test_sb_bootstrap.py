"""Parametric bootstrap checks against an independent Gaussian-mean model."""
from types import SimpleNamespace
import json
from pathlib import Path
import tempfile

import numpy as np
from numpy.testing import assert_allclose
from scipy.optimize import OptimizeResult

from sb_bootstrap import bootstrap_fit, validate_bootstrap


class MeanModel:
    """Exact MLE for four Gaussian observations; no nonlinear solver needed."""
    parameter_names = ['kab', 'kba']

    def __init__(self):
        self.es = SimpleNamespace(int=[1., 1., 1., 1.], intstd=[.2]*4,
                                  offset=[0., 1., 2., 3.])
        self.dataset = SimpleNamespace(res=[SimpleNamespace(active=True, estSpecs=[self.es])])
        self.config = {'init': {'vary': ['kab']}}
        self.free = np.array([0])
        self.lower = np.array([0., 1e-8])
        self.upper = np.array([np.inf, np.inf])
        self.initial_parameters = np.array([1., 8.])
        self._fit_data = None
        self.verbose = False
        self.calls = 0

    def _prepare_data(self):
        return [(0, self.es, np.asarray(self.es.offset), np.asarray(self.es.int), np.asarray(self.es.intstd))]

    def seParam(self, p):
        return p

    def calc(self, p, index, offsets, es):
        return np.full(len(offsets), p[0])

    def errFunc(self, p):
        residual = (np.asarray(self.es.int)-p[0])/.2
        self.chi2 = float(residual @ residual)
        return residual

    def fit(self, p0=None, fitting_config=None):
        self.calls += 1
        p = np.array([np.mean(self.es.int), 8.])
        self.nvar, self.npar, self.dof = 4, 1, 3
        self.rank, self.condition = 1, 1.
        self.errFunc(p)
        self.result = OptimizeResult(success=True, nfev=1, message='Exact Gaussian mean', status=1)
        return p, np.diag([.01, 0.])


def check_reproducibility_and_resume():
    model = MeanModel()
    p, _ = model.fit()
    original_result, original_y = model.result, list(model.es.int)
    settings = {'replicates': 80, 'seed': 765, 'confidence': .95}
    first = bootstrap_fit(model, p, settings)
    assert first['successful'] == 80
    assert first['intervals']['kba']['fixed']
    assert first['intervals']['kba']['lower'] == first['intervals']['kba']['upper'] == 8.
    assert first['intervals']['kab']['lower'] < 1. < first['intervals']['kab']['upper']
    assert_allclose(np.std([x['parameters'][0] for x in first['samples']], ddof=1), .1, rtol=.2)
    assert model.result is original_result and model.es.int == original_y
    completed = []

    def stop(row):
        completed.append(row)
        if len(completed) == 3:
            raise KeyboardInterrupt

    try:
        bootstrap_fit(model, p, settings, on_complete=stop)
    except KeyboardInterrupt:
        pass
    else:
        raise AssertionError('Interruption swallowed')
    assert model.result is original_result and model.es.int == original_y
    before = model.calls
    resumed = bootstrap_fit(model, p, settings, completed=completed)
    assert model.calls-before == 77
    assert resumed == first
    different = bootstrap_fit(model, p, {**settings, 'seed': 766})
    assert different['samples'][0]['parameters'] != first['samples'][0]['parameters']


def check_invalid_and_failure():
    for setting in ({}, {'replicates': True, 'seed': 1}, {'replicates': 2},
                    {'replicates': 2, 'seed': -1}, {'replicates': 2, 'seed': 1, 'confidence': 1}):
        try:
            validate_bootstrap(setting)
        except ValueError:
            pass
        else:
            raise AssertionError(f'Invalid settings accepted: {setting}')
    model = MeanModel()
    p, _ = model.fit()
    model.fit = lambda **kwargs: (_ for _ in ()).throw(RuntimeError('Intentional fit failure'))
    result = bootstrap_fit(model, p, {'replicates': 2, 'seed': 1})
    assert result['successful'] == 0 and result['warnings']
    assert all(not row['success'] and row['parameters'] is None for row in result['samples'])
    assert result['intervals']['kab']['lower'] is None
    assert model.es.int == [1.]*4


def check_fixed_zero_population():
    model = MeanModel()
    model.fit()
    model.free = np.array([1])
    model.config['init']['vary'] = ['kba']
    samples = [{'index': i, 'success': True, 'parameters': [0., rate],
                'stderr': [0., .1], 'derived_se': {'kex': .1, 'pB': 0.},
                'kex': rate, 'pB': 0., 'chi2': 1., 'message': 'converged', 'at_bounds': []}
               for i, rate in enumerate((7., 9.))]
    result = bootstrap_fit(model, np.array([0., 8.]), {'replicates': 2, 'seed': 1}, completed=samples)
    assert result['intervals']['pB']['fixed'], 'Fixed zero kab makes pB identically zero'
    assert not result['intervals']['kex']['fixed']


def check_real_cli_bootstrap():
    from sbfit import run_config
    from test_sb_output import small_config

    folder = Path(tempfile.mkdtemp(prefix='sbonest-bootstrap-cli-'))
    config = small_config(folder)
    config['init']['bootstrap'] = {'replicates': 3, 'seed': 789}
    run_config(config, no_pdf=True)
    result = json.loads((folder / 'fit_result.json').read_text())
    assert result['bootstrap']['successful'] == 3
    assert len(list((folder / 'fit_checkpoint').glob('bootstrap-*.json'))) == 3
    assert result['bootstrap']['intervals']['v1n_scale']['lower'] > .8
    assert result['bootstrap']['intervals']['v1n_scale']['upper'] < 1.2
    assert result['bootstrap']['intervals']['kab']['fixed']
    assert result['chi2'] > 0 and result['n_points'] == 294


if __name__ == '__main__':
    check_reproducibility_and_resume()
    check_invalid_and_failure()
    check_fixed_zero_population()
    check_real_cli_bootstrap()
    print('PASS: bootstrap Gaussian reference, seed, resume, failures and state restoration')
