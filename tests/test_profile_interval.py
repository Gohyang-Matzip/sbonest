"""Likelihood-ratio profile intervals: settings, root finding, replay and failure handling."""
# ruff: noqa: E402 -- Limit numerical libraries before importing them.
from pathlib import Path
import _env  # noqa: F401


import contextlib
import io
import json
import tempfile
from unittest.mock import patch

import numpy as np

import sb_analysis
from sb_analysis import _interval_settings, profile_intervals
from sbfit import SidebandModel, run_config
from test_sb_output import small_config


def check_settings():
    names = ['kab', 'kba', 'v1n_scale', 'A1.R2b']
    assert _interval_settings({'parameters': ['kex']}, names)['confidence'] == 0.95
    cases = [{}, {'parameters': []}, {'parameters': ['kex', 'kex']}, {'parameters': ['A1.R2b']},
             {'parameters': ['kex'], 'confidence': 1.}, {'parameters': ['kex'], 'max_evaluations': 3},
             {'parameters': ['kex'], 'relative_tolerance': 0.}, {'parameters': ['kex'], 'max_doublings': 0},
             {'parameters': ['kex'], 'other': 1}]
    for bad in cases:
        try:
            _interval_settings(bad, names)
        except ValueError:
            pass
        else:
            raise AssertionError(f'Accepted invalid profile_interval settings {bad}')
    try:
        _interval_settings({'parameters': ['v1n_scale']}, ['kab', 'kba', 'v1n[0]'])
    except ValueError:
        pass
    else:
        raise AssertionError('v1n_scale interval accepted without scale RF mode')


def fitted_model(folder):
    cfg = small_config(folder)
    cfg['init']['vary'] = ['kab', 'kba', 'v1n_scale', 'A1.R2b']
    model = SidebandModel(cfg)
    model.verbose = False
    p, covariance = model.fit(fitting_config=cfg['init'])
    model._fit_data = model._prepare_data()
    return cfg, model, p, covariance


def check_intervals_match_local_errors():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-interval-'))
    cfg, model, p, covariance = fitted_model(folder)
    saved = []
    out = profile_intervals(model, p, covariance, {'parameters': ['kex', 'pB', 'v1n_scale']},
                            on_complete=lambda name, row: saved.append((name, row)))
    assert out['threshold_delta_chi2'] > 3.84 and out['base_chi2'] > 0
    for name in ('kex', 'pB', 'v1n_scale'):
        result = out[name]
        assert result['success'], result
        assert result['lower'] < result['estimate'] < result['upper']
        half = 1.959963984540054 * result['local_se']
        for side in ((result['estimate'] - result['lower']), (result['upper'] - result['estimate'])):
            assert 0.7 < side / half < 1.4, (name, side / half)
        # The located crossings sit at the threshold within the refit tolerance.
        crossings = [row for row in result['points'] if row['target'] in (result['lower'], result['upper'])]
        assert len(crossings) == 2
        assert all(abs(row['delta_chi2'] - out['threshold_delta_chi2']) < 0.05 for row in crossings)
        assert result['evaluations'] == len(result['points']) == len([r for n, r in saved if n == name])
    assert not out['warnings']
    assert np.array_equal(np.asarray(model.errFunc(p)), np.asarray(model.errFunc(p)))
    json.dumps(out, allow_nan=False)
    # Replaying completed rows evaluates nothing and reproduces the same intervals.
    completed = {name: [row for n, row in saved if n == name] for name in ('kex', 'pB', 'v1n_scale')}
    with patch.object(sb_analysis, 'profile_point', side_effect=AssertionError('replay refitted')):
        again = profile_intervals(model, p, covariance, {'parameters': ['kex', 'pB', 'v1n_scale']}, completed=completed)
    assert {k: v for k, v in again.items()} == out
    # A completed row that does not match the deterministic sequence is rejected.
    corrupted = {'kex': [dict(completed['kex'][0], target=completed['kex'][0]['target'] + 1.)]}
    try:
        profile_intervals(model, p, covariance, {'parameters': ['kex']}, completed=corrupted)
    except ValueError as exc:
        assert 'evaluation sequence' in str(exc)
    else:
        raise AssertionError('Mismatched completed interval rows were accepted')


def check_open_intervals_and_failures():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-interval-open-'))
    cfg, model, p, covariance = fitted_model(folder)
    tiny = np.asarray(covariance) * 1e-6
    out = profile_intervals(model, p, tiny, {'parameters': ['kex'], 'max_doublings': 1, 'max_evaluations': 6})
    result = out['kex']
    assert not result['success'] and result['lower'] is None and result['upper'] is None
    assert 'no threshold crossing' in result['lower_message'] and any('one-sided or open' in w for w in out['warnings'])
    assert result['evaluations'] == 4
    # A failed nuisance refit is recorded, not raised, and later parameters still run.
    with patch.object(sb_analysis, 'profile_point', side_effect=lambda *a, **k: {
            'target': float(a[3]), 'start': None, 'parameters': None, 'success': False, 'chi2': None,
            'delta_chi2': None, 'below_base_minimum': False, 'nuisance_parameters': [],
            'optimization_performed': True, 'nfev': None, 'status': None, 'message': 'boom'}):
        failed = profile_intervals(model, p, covariance, {'parameters': ['kex', 'pB']})
    assert all('boom' in failed[name]['message'] and not failed[name]['success'] for name in ('kex', 'pB'))
    assert len(failed['warnings']) == 2 and failed['pB']['evaluations'] == 1
    # Without a finite local error the search starts from a 10% step and warns.
    nan_cov = np.full_like(np.asarray(covariance), np.nan)
    out = profile_intervals(model, p, nan_cov, {'parameters': ['v1n_scale']})
    assert any('No finite local standard error' in w for w in out['warnings']) and out['v1n_scale']['success']


def check_run_config_and_report():
    from sb_report import regenerate_report

    folder = Path(tempfile.mkdtemp(prefix='sbonest-interval-run-'))
    cfg = small_config(folder)
    cfg['init']['vary'] = ['kab', 'kba', 'v1n_scale', 'A1.R2b']
    cfg['init']['profile_interval'] = {'parameters': ['kex'], 'confidence': 0.9}
    with contextlib.redirect_stdout(io.StringIO()):
        run_config(cfg, no_pdf=True)
    result = json.loads((folder / 'fit_result.json').read_text())
    interval = result['profile_intervals']['kex']
    assert interval['success'] and result['profile_intervals']['confidence'] == 0.9
    records = sorted(path.name for path in (folder / 'fit_checkpoint').glob('profile_interval-kex-*.json'))
    assert len(records) == interval['evaluations']
    with patch.object(SidebandModel, 'fit', side_effect=AssertionError('resume refitted')), \
            patch.object(sb_analysis, 'least_squares', side_effect=AssertionError('resume refit profile')):
        run_config(cfg, no_pdf=True, resume=True)
    paths = regenerate_report(folder / 'fit_result.json', folder / 'report')
    text = Path(paths['summary_txt']).read_text()
    assert 'Profile interval kex' in text and f'{interval["lower"]:.8g}' in text
    assert Path(paths['pdf']).stat().st_size > 1000


if __name__ == '__main__':
    check_settings()
    check_intervals_match_local_errors()
    check_open_intervals_and_failures()
    check_run_config_and_report()
    print('PASS: profile interval settings, root finding, replay, failures, run and report')
