"""Executable checks for Sideband uncertainty, provenance and lossless outputs."""
import copy
import csv
import hashlib
import json
import tempfile
from pathlib import Path
from unittest.mock import patch

import numpy as np
from numpy.testing import assert_allclose

from run import load_config
from sbfit import SidebandModel, run_config
from sb_report import provenance

ROOT = Path(__file__).resolve().parent


def small_config(folder):
    cfg = load_config(ROOT / 'example/sideband_auto_H/two_RF.json')
    cfg['Project Name'] = str(folder / 'fit')
    cfg['residues'] = [{'name': n, 'flag': 'on' if n == 'A1' else 'off'}
                       for n in ('A1', 'G2', 'S3')]
    cfg['init'].update(kex={'min': 300, 'max': 300, 'nsteps': 1},
                       pB={'min': .05, 'max': .05, 'nsteps': 1},
                       vary=['v1n_scale'])
    cfg['init']['initial'] = dict(zip(
        ('A1.peak_ppm', 'A1.dw_ppm', 'A1.R1', 'A1.R2a', 'A1.R2b'),
        (120., 3., 1.5, 12., 15.)))
    return cfg


def check_covariance():
    cfg = small_config(Path(tempfile.mkdtemp(prefix='sbonest-cov-')))
    model = SidebandModel(cfg)
    p, cov = model.fit()
    # An independent correlated-rate covariance catches diagonal-only propagation.
    cov[:2, :2] = [[4., -1.5], [-1.5, 9.]]
    info = model.diagnostics(p, cov)
    assert 'covariance' in info, 'Full covariance is missing from the result'
    assert info['parameter_order'] == model.parameter_names
    assert_allclose(info['covariance'], cov)
    rate_gradient = np.array([p[1], -p[0]]) / (p[0] + p[1]) ** 2
    assert_allclose(info['derived_se']['kex'], np.sqrt(10.))
    assert_allclose(info['derived_se']['pB'], np.sqrt(rate_gradient @ cov[:2, :2] @ rate_gradient))
    cov[:2, :2] = np.nan
    unavailable = model.diagnostics(p, cov)
    assert unavailable['derived_se'] == {'kex': None, 'pB': None}
    assert unavailable['covariance'][0][0] is None
    json.dumps(unavailable, allow_nan=False)
    cov[:] = 0
    assert model.diagnostics(p, cov)['derived_se'] == {'kex': 0., 'pB': 0.}


def check_provenance_hashes():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-hashes-'))
    cfg = small_config(folder)
    # Metadata reads the effective waveform after per-dataset overrides.
    waveform = folder / 'pulse.json'
    waveform.write_text('{"waveform": "first bytes"}\n')
    cfg['sideband']['datasets'] = [{'waveform_json': 'pulse.json'}, {}]
    first = provenance(cfg, folder)
    assert first['waveforms'] == [{'path': str(waveform.resolve()),
                                 'sha256': hashlib.sha256(waveform.read_bytes()).hexdigest()}]
    for entry in first['datasets']:
        assert entry['sha256'] == hashlib.sha256(Path(entry['path']).read_bytes()).hexdigest()
    waveform.write_text('{"waveform": "changed bytes"}\n')
    second = provenance(cfg, folder)
    assert second['waveforms'][0]['sha256'] != first['waveforms'][0]['sha256']
    assert second['config_sha256'] == first['config_sha256']
    cfg['init']['initial']['A1.R1'] = 1.6
    assert provenance(cfg, folder)['config_sha256'] != first['config_sha256']


def check_partial_fit_parity():
    cfg = load_config(ROOT / 'example/sideband_auto_H/two_RF.json')
    model = SidebandModel(cfg)
    reference = json.loads((ROOT / 'results/auto_H_refit/fits/two_RF_result.json').read_text())
    initial = np.array([reference['parameters'][name]['value'] for name in model.parameter_names])
    cfg['init']['initial'] = {}
    cfg['init']['vary'] = [name for name in model.parameter_names if name != 'A1.R1H']
    # Independent prior optimizer path: SciPy evaluates each 3-point column.
    with patch('sbfit._block_jacobian', return_value='3-point'):
        expected, expected_cov = model.fit(initial)
    expected_chi2, expected_rank = model.chi2, model.rank
    actual, actual_cov = model.fit(initial)
    assert_allclose(actual, expected, rtol=1e-8, atol=1e-8)
    assert_allclose(model.chi2, expected_chi2, rtol=1e-8, atol=1e-8)
    assert model.rank == expected_rank
    assert_allclose(actual_cov, expected_cov, rtol=1e-7, atol=1e-9, equal_nan=True)


def check_cli_outputs():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-output-'))
    cfg = small_config(folder)
    run_config(cfg, no_pdf=False)
    info = json.loads((folder / 'fit_result.json').read_text())
    assert 'provenance' in info, 'General CLI provenance is missing'
    meta = info['provenance']
    assert meta['elapsed_s'] > 0
    assert set(meta['packages']) >= {'numpy', 'scipy', 'optimalcontrol-nmr'}
    assert len(meta['datasets']) == 2
    assert all(len(row['sha256']) == 64 for row in meta['datasets'])
    assert 'sbfit.py' in meta['source_sha256']
    assert len(meta['config_sha256']) == 64
    with (folder / 'fit_predictions.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == info['n_points']
    assert {row['residue'] for row in rows} == {'A1'}
    assert_allclose(sum(float(row['residual_sigma']) ** 2 for row in rows), info['chi2'], rtol=1e-12)
    for row in rows:
        assert_allclose((float(row['observed']) - float(row['predicted'])) / float(row['sigma']),
                        float(row['residual_sigma']), rtol=1e-12, atol=1e-12)
    for name in ('fit.pdf', 'fit_data.pdf'):
        assert (folder / name).read_bytes().startswith(b'%PDF')
    before = {str(p.relative_to(folder)): p.read_bytes() for p in folder.rglob("*") if p.is_file()}
    try:
        run_config(cfg, no_pdf=True)
    except FileExistsError:
        pass
    else:
        raise AssertionError('Existing result was overwritten')
    assert before == {str(p.relative_to(folder)): p.read_bytes() for p in folder.rglob("*") if p.is_file()}
    blocked = copy.deepcopy(cfg)
    blocked['Project Name'] = str(folder / 'reserved')
    (folder / 'reserved_predictions.csv').write_text('preserve me')
    try:
        run_config(blocked, no_pdf=True)
    except FileExistsError:
        pass
    else:
        raise AssertionError('Existing CSV was not protected')
    print('Output artifacts:', folder)


def check_symlink_protection():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-symlink-'))
    cfg = small_config(folder)
    target = folder / 'must-not-be-created.txt'
    (folder / 'fit_result.txt').symlink_to(target)
    try:
        run_config(cfg, no_pdf=True)
    except FileExistsError:
        pass
    else:
        raise AssertionError('Dangling output symlink was followed')
    assert not target.exists()
    assert not (folder / 'fit_result.json').exists()


def check_cli_analyses():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-analysis-output-'))
    cfg = small_config(folder)
    cfg['init']['multistart'] = {'starts': [{'v1n_scale': 1.1}],
                                'random_starts': 1, 'seed': 4321}
    cfg['init']['profile'] = {'v1n_scale': [1.06, 1.08, 1.10]}
    run_config(cfg, no_pdf=True)
    info = json.loads((folder / 'fit_result.json').read_text())
    assert 'multistart' in info, 'CLI did not execute requested restarts'
    assert len(info['multistart']) == 3
    assert sum(row['selected'] for row in info['multistart']) == 1
    selected = next(row for row in info['multistart'] if row['selected'])
    assert_allclose(selected['chi2'], info['chi2'], rtol=1e-12)
    assert len(info['profiles']['v1n_scale']) == 3
    assert all(row['success'] for row in info['profiles']['v1n_scale'])
    assert info['profiles']['base_chi2'] == info['chi2']
    failed = copy.deepcopy(cfg)
    failed['Project Name'] = str(folder / 'failed')
    failed['init']['max_nfev'] = 1
    try:
        run_config(failed, no_pdf=True)
    except RuntimeError:
        pass
    else:
        raise AssertionError('Unconverged attempts were accepted')
    failure = json.loads((folder / 'failed_result.json').read_text())
    assert failure['success'] is False
    assert len(failure['multistart']) == 3
    assert all(not row['success'] for row in failure['multistart'])
    assert 'provenance' in failure
    assert not (folder / 'failed_predictions.csv').exists()


if __name__ == '__main__':
    for check in (check_covariance, check_provenance_hashes, check_cli_outputs, check_symlink_protection,
                  check_cli_analyses, check_partial_fit_parity):
        check()
        print('PASS', check.__name__)
