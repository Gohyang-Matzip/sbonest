"""Residual and sigma diagnostics on synthetic residual rows and a real small fit."""
# ruff: noqa: E402 -- Limit numerical libraries before importing them.
import os

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(name, '1')
os.environ.setdefault('MPLBACKEND', 'Agg')

import contextlib
import io
import json
import math
from pathlib import Path
import tempfile

import numpy as np

from sb_diagnostics import diagnostics_lines, lag1_autocorrelation, residual_diagnostics, runs_test
from sbfit import run_config
from test_sb_output import small_config


def rows_from(values, *, residue='A1', dataset=0):
    return [{'residue': residue, 'dataset_index': dataset, 'offset_ppm': float(i), 'residual_sigma': float(v)}
            for i, v in enumerate(values)]


def check_tests():
    alternating = runs_test([1, -1] * 20)
    assert alternating['runs'] == 40 and alternating['p_value'] < 1e-6 and alternating['z'] > 0
    clustered = runs_test([1] * 20 + [-1] * 20)
    assert clustered['runs'] == 2 and clustered['p_value'] < 1e-6 and clustered['z'] < 0
    assert runs_test([1, 1, 1])['p_value'] is None and runs_test([])['runs'] == 0
    rng = np.random.default_rng(1)
    random = runs_test(rng.normal(size=400))
    assert 0.01 < random['p_value'] <= 1.0
    assert abs(lag1_autocorrelation(np.sin(np.linspace(0, 20, 200)))) > 0.9
    assert lag1_autocorrelation([1., 1.]) is None and lag1_autocorrelation([2., 2., 2.]) is None


def check_synthetic_rows():
    rng = np.random.default_rng(7)
    good = rows_from(rng.normal(size=300)) + rows_from(rng.normal(size=300), residue='G2', dataset=1)
    cov = np.diag([4., 9.])
    report = residual_diagnostics(good, 2, cov)
    assert abs(report['reduced_chi2'] - 1.) < 3 * report['reduced_chi2_expected_sd'] and not report['warnings']
    assert report['rescaled_stderr'] == [math.sqrt(4. * report['reduced_chi2']), math.sqrt(9. * report['reduced_chi2'])]
    assert set(report['by_residue']) == {'A1', 'G2'} and set(report['by_dataset']) == {'0', '1'}
    assert len(report['blocks']) == 2 and report['dof'] == 598
    json.dumps(report, allow_nan=False)
    text = '\n'.join(diagnostics_lines(report, ['kab', 'kba']))
    assert 'reduced chi2' in text and 'Rescaled standard errors' in text and 'kab:' in text
    small_sigma = residual_diagnostics(rows_from(rng.normal(scale=2., size=400)), 3)
    assert any('too small' in w for w in small_sigma['warnings']) and small_sigma['sigma_scale_estimate'] > 1.7
    large_sigma = residual_diagnostics(rows_from(rng.normal(scale=.5, size=400)), 3)
    assert any('too large' in w for w in large_sigma['warnings'])
    structured = residual_diagnostics(rows_from(np.sin(np.linspace(0, 6 * np.pi, 300)) * 1.3 + rng.normal(scale=.2, size=300)), 3)
    assert any('Runs test' in w for w in structured['warnings']) and any('autocorrelation' in w for w in structured['warnings'])
    heavy = residual_diagnostics(rows_from(np.r_[rng.normal(size=400), [5., -6., 7., 8., -9.]]), 2)
    assert any('beyond 3 sigma' in w for w in heavy['warnings'])
    for bad in ((rows_from([1., 2.]), 2), ([], 1)):
        try:
            residual_diagnostics(*bad)
        except ValueError:
            pass
        else:
            raise AssertionError('Insufficient degrees of freedom accepted')


def check_fit_integration():
    from sb_report import regenerate_report

    folder = Path(tempfile.mkdtemp(prefix='sbonest-diagnostics-'))
    cfg = small_config(folder)
    cfg['init']['vary'] = ['kab', 'kba', 'v1n_scale', 'A1.R2b']
    with contextlib.redirect_stdout(io.StringIO()):
        run_config(cfg, no_pdf=True)
    result = json.loads((folder / 'fit_result.json').read_text())
    report = result['residual_diagnostics']
    assert report['n_points'] == result['n_points'] and report['dof'] == result['dof']
    assert abs(report['chi2'] - result['chi2']) < 1e-9 * result['chi2']
    kab = result['parameters']['kab']
    assert abs(kab['stderr_rescaled'] - kab['stderr'] * math.sqrt(report['reduced_chi2'])) < 1e-12
    assert result['parameters']['A1.R1']['stderr_rescaled'] is None
    text = (folder / 'fit_result.txt').read_text()
    assert 'Residual diagnostics' in text and 'Rescaled standard errors' in text
    paths = regenerate_report(folder / 'fit_result.json', folder / 'report')
    summary = json.loads(Path(paths['summary_json']).read_text())
    assert abs(summary['residual_diagnostics']['reduced_chi2'] - report['reduced_chi2']) < 1e-12
    assert 'overall runs test' in Path(paths['summary_txt']).read_text()


if __name__ == '__main__':
    check_tests()
    check_synthetic_rows()
    check_fit_integration()
    print('PASS: runs test, autocorrelation, sigma scale warnings, fit integration and report')
