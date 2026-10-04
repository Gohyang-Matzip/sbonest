"""Parity and safety checks for the optional worker pool.

Every parallel result must equal the serial result exactly: the pool only
changes which process evaluates a residual, never the arithmetic.
"""
# ruff: noqa: E402 -- Limit numerical libraries before importing them.
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]  # repository root: flat modules live there
sys.path.insert(0, str(ROOT))

import os

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(name, '1')
os.environ.setdefault('MPLBACKEND', 'Agg')

import contextlib
import io
import json
import subprocess
import tempfile
from unittest.mock import patch

import numpy as np
from numpy.testing import assert_array_equal

from sb_analysis import identifiability, local_jacobian
from sb_parallel import Progress, WorkerPool, validate_workers
from sbfit import SidebandModel, run_config
from test_sb_output import small_config



def analysis_config(folder, name):
    cfg = small_config(folder)
    cfg['Project Name'] = str(folder / name)
    cfg['init']['vary'] = ['kab', 'kba', 'v1n_scale', 'A1.R2b', 'A1.R1H', 'A1.R2H']
    cfg['init']['multistart'] = {'starts': [{'kab': 20., 'kba': 380.}], 'random_starts': 1, 'seed': 7}
    cfg['init']['profile'] = {'v1n_scale': [1.06], 'kex': [290.]}
    cfg['init']['bootstrap'] = {'replicates': 2, 'seed': 11}
    return cfg


def strip(result):
    return {key: value for key, value in result.items()
            if key not in ('provenance', 'checkpoint', 'config')}


def check_validation_and_progress():
    for bad in (0, -1, True, 1.5, '2'):
        try:
            validate_workers(bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f'Accepted invalid worker count {bad!r}')
    assert validate_workers(3) == 3
    stream = io.StringIO()
    with contextlib.redirect_stdout(stream):
        progress = Progress('demo', 3, enabled=True, done=1)
        progress.step()
        progress.step()
        Progress('silent', 2, enabled=False).step()
    lines = stream.getvalue().splitlines()
    assert len(lines) == 2 and lines[0].startswith('[demo] 2/3 done') and 'ETA' in lines[1]
    assert all('silent' not in line for line in lines)


def check_jacobian_and_identifiability_parity():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-parallel-jac-'))
    cfg = small_config(folder)
    cfg['init']['vary'] = ['kab', 'kba', 'v1n_scale', 'A1.R2b', 'A1.R1H']
    model = SidebandModel(cfg)
    p0 = model.prepare_fit()
    model._fit_data = model._prepare_data()
    serial = local_jacobian(model, p0)
    report = identifiability(model, p0)
    with WorkerPool(2, cfg) as pool:
        parallel = local_jacobian(model, p0, pool=pool)
        pooled = identifiability(model, p0, pool=pool)
        try:
            pool.evaluate_many([np.full(len(p0), np.nan)])
        except ValueError:
            pass  # The same ValueError the serial residual raises for this vector.
        else:
            raise AssertionError('Worker failure did not propagate')
        # The pool is still usable after a task error.
        again = pool.evaluate_many([p0, p0])
    assert_array_equal(serial, parallel)
    assert serial.shape == (294, 5)
    assert_array_equal(again[0], again[1])
    assert_array_equal(again[0], model.errFunc(p0))
    assert report == pooled
    assert report['jacobian_rank'] <= 5 and report['n_free'] == 5
    assert set(report['expected_se']) == set(cfg['init']['vary'])
    assert report['relative_se']['A1.R1H'] is None or report['relative_se']['A1.R1H'] >= 0
    json.dumps(report, allow_nan=False)
    # A zero-sensitivity column is reported, not hidden.
    with patch.object(SidebandModel, 'errFunc', lambda self, p: np.ones(294)):
        flat = identifiability(model, p0)
    assert flat['jacobian_rank'] == 0 and set(flat['zero_sensitivity']) == set(cfg['init']['vary'])
    assert set(flat['weak_parameters']) == set(cfg['init']['vary'])


def check_run_config_parity_and_resume():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-parallel-run-'))
    serial_cfg, parallel_cfg = analysis_config(folder, 'serial'), analysis_config(folder, 'parallel')
    stream = io.StringIO()
    with contextlib.redirect_stdout(stream):
        run_config(serial_cfg, no_pdf=True, workers=1)
        run_config(parallel_cfg, no_pdf=True, workers=2)
    serial = json.loads((folder / 'serial_result.json').read_text())
    parallel = json.loads((folder / 'parallel_result.json').read_text())
    assert strip(serial) == strip(parallel), 'Parallel analysis differs from the serial analysis'
    assert serial['provenance']['workers'] == 1 and parallel['provenance']['workers'] == 2
    assert len(serial['multistart']) == 3 and all(len(serial['profiles'][k]) == 1 for k in ('v1n_scale', 'kex'))
    assert serial['bootstrap']['successful'] <= 2 and len(serial['bootstrap']['samples']) == 2
    output = stream.getvalue()
    assert '[bootstrap] 2/2 done' in output and '[profile] 2/2 done' in output and '[multistart] 3/3 done' in output
    names = sorted(p.name for p in (folder / 'parallel_checkpoint').glob('*.json'))
    assert names == sorted(p.name for p in (folder / 'serial_checkpoint').glob('*.json'))
    assert {'attempt-1.json', 'attempt-2.json', 'profile-kex-0.json', 'bootstrap-1.json'} <= set(names)
    for name in names:
        a = json.loads((folder / 'serial_checkpoint' / name).read_text())
        b = json.loads((folder / 'parallel_checkpoint' / name).read_text())
        if name not in ('manifest.json', 'provenance.json', 'exports.json'):
            assert a == b, f'Checkpoint record {name} differs between serial and parallel runs'
    # Worker count is not part of the checkpoint identity: a parallel run resumes serially.
    with patch.object(SidebandModel, 'fit', side_effect=AssertionError('Resume refitted')):
        run_config(parallel_cfg, no_pdf=True, resume=True, workers=1)
        run_config(serial_cfg, no_pdf=True, resume=True, workers=3)


def check_cli_flags():
    env = dict(os.environ, PYTHONPATH=str(ROOT))
    config = ROOT / 'example/sideband_auto_H/two_RF.json'
    for args in (['--workers', '0'], ['--identifiability'], ['--check', '--resume']):
        run = subprocess.run([sys.executable, 'run.py', str(config), *args], cwd=ROOT, env=env,
                             capture_output=True, text=True)
        assert run.returncode == 2, (args, run.stderr)
    run = subprocess.run([sys.executable, 'run.py', str(config), '--check', '--identifiability',
                          '--no-pdf', '--workers', '2'], cwd=ROOT, env=env, capture_output=True, text=True)
    assert run.returncode in (0, 1), run.stderr
    summary = json.loads(run.stdout)
    report = summary['identifiability']
    assert report['n_free'] == 24 and len(report['scaled_singular_values']) == 24
    assert any('initial point' in warning for warning in summary['warnings'])
    assert not any(Path(path).exists() for path in summary['outputs'])


if __name__ == '__main__':
    check_validation_and_progress()
    check_jacobian_and_identifiability_parity()
    check_run_config_parity_and_resume()
    check_cli_flags()
    print('PASS: worker pool parity, identifiability diagnostics and CLI flags')
