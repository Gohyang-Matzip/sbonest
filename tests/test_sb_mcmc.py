"""MCMC checks: an exact Gaussian posterior, reproducibility, resume and serial/pool identity."""
# ruff: noqa: E402 -- Limit numerical libraries before importing them.
from pathlib import Path
import _env  # noqa: F401

import json
import tempfile

import numpy as np
from numpy.testing import assert_allclose

from sb_mcmc import mcmc_fit, validate_mcmc
from test_sb_bootstrap import MeanModel


def check_validation():
    assert validate_mcmc({'steps': 400, 'seed': 1}) == {'steps': 400, 'seed': 1, 'chains': 4, 'burn_in': 100,
                                                        'confidence': .95}
    for bad in ({}, {'steps': 400}, {'steps': 50, 'seed': 1}, {'steps': 400, 'seed': -1}, {'steps': 400, 'seed': True},
                {'steps': 400, 'seed': 1, 'chains': 1}, {'steps': 400, 'seed': 1, 'burn_in': 400},
                {'steps': 400, 'seed': 1, 'confidence': 1.}, {'steps': 400, 'seed': 1, 'walkers': 3}, []):
        try:
            validate_mcmc(bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f'Accepted invalid mcmc settings {bad!r}')


def check_gaussian_posterior_and_resume():
    """Four observations of 1 with sigma 0.2: the posterior of the mean is N(1, 0.1^2)."""
    model = MeanModel()
    p, covariance = model.fit()
    original = model.result
    settings = {'steps': 6000, 'seed': 31, 'chains': 4}
    rows = []
    first = mcmc_fit(model, p, covariance, settings, on_complete=rows.append)
    item = first['parameters']['kab']
    assert_allclose([item['median'], item['mean']], 1., atol=.01)
    assert_allclose(item['sd'], .1, rtol=.08)
    assert_allclose([item['lower'], item['upper']], [1 - .196, 1 + .196], atol=.02)
    assert item['rhat'] < 1.02 and item['ess'] > 400 and not first['warnings'], first
    assert all(.1 < rate < .7 for rate in first['accept_rates'])
    assert set(first['parameters']) == {'kab', 'kex', 'pB'} and first['draws_per_chain'] == len(rows[0]['samples'])
    assert model.result is original
    json.dumps(first, allow_nan=False)
    # The same seed reproduces every draw; a completed prefix is reused, not resampled.
    again = []
    assert mcmc_fit(model, p, covariance, settings, on_complete=again.append) == first and again == rows
    resumed = []
    assert mcmc_fit(model, p, covariance, settings, completed=rows[:2], on_complete=resumed.append) == first
    assert resumed == rows[2:]
    assert mcmc_fit(model, p, covariance, {**settings, 'seed': 32})['parameters'] != first['parameters']
    for broken in ([{**rows[0], 'index': 1}], [{**rows[0], 'samples': [[float('nan')]]}], rows + rows):
        try:
            mcmc_fit(model, p, covariance, settings, completed=broken)
        except ValueError:
            pass
        else:
            raise AssertionError('Accepted an invalid completed mcmc prefix')
    # A singular local covariance cannot define the proposal.
    try:
        mcmc_fit(model, p, np.zeros((2, 2)), settings)
    except ValueError as exc:
        assert 'positive-definite' in str(exc)
    else:
        raise AssertionError('Singular covariance accepted')
    # Too few steps are flagged rather than silently reported.
    short = mcmc_fit(model, p, covariance, {'steps': 100, 'seed': 31, 'chains': 2})
    assert any('effective sample size' in text for text in short['warnings'])


def check_run_serial_pool_identity_and_report():
    from sb_report import regenerate_report
    from sbfit import run_config
    from test_sb_output import small_config

    folder = Path(tempfile.mkdtemp(prefix='sbonest-mcmc-'))
    results = []
    for name, workers in (('serial', 1), ('pooled', 2)):
        cfg = small_config(folder)
        cfg['Project Name'] = str(folder / name)
        cfg['init']['vary'] = ['kab', 'kba', 'v1n_scale']
        cfg['init']['mcmc'] = {'steps': 120, 'seed': 5, 'chains': 2}
        run_config(cfg, no_pdf=True, workers=workers)
        results.append(json.loads((folder / f'{name}_result.json').read_text()))
        chains = sorted((folder / f'{name}_checkpoint').glob('mcmc-*.json'))
        assert len(chains) == 2
        results[-1]['chains'] = [json.loads(path.read_text())['value'] for path in chains]
    serial, pooled = results
    assert serial['mcmc'] == pooled['mcmc'] and serial['chains'] == pooled['chains']
    assert set(serial['mcmc']['parameters']) == {'kab', 'kba', 'v1n_scale', 'kex', 'pB'}
    assert 'samples' not in serial['mcmc']
    paths = regenerate_report(folder / 'serial_result.json', folder / 'report')
    text = Path(paths['summary_txt']).read_text()
    assert 'MCMC: 2 chains x 120 steps' in text and 'R-hat=' in text
    # Resume publishes the completed outputs without sampling again.
    cfg['Project Name'] = str(folder / 'serial')
    run_config(cfg, no_pdf=True, resume=True)
    assert json.loads((folder / 'serial_result.json').read_text())['mcmc'] == serial['mcmc']
    bad = small_config(folder)
    bad['Project Name'] = str(folder / 'bad')
    bad['init']['mcmc'] = {'steps': 10, 'seed': 5}
    try:
        run_config(bad, no_pdf=True)
    except ValueError as exc:
        assert 'mcmc.steps' in str(exc) and not (folder / 'bad_checkpoint').exists()
    else:
        raise AssertionError('Invalid mcmc settings reached the fit')


if __name__ == '__main__':
    check_validation()
    check_gaussian_posterior_and_resume()
    check_run_serial_pool_identity_and_report()
    print('PASS: MCMC validation, exact Gaussian posterior, resume, serial/pool identity and report')
