"""Unified command line and packaging metadata."""
# ruff: noqa: E402 -- Limit numerical libraries before importing them.
import os

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(name, '1')
os.environ.setdefault('MPLBACKEND', 'Agg')

import importlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import tomllib

from sb_compare import _residue_config
from test_sb_compare import write_synthetic

ROOT = Path(__file__).resolve().parent


def run_cli(*args, cwd=None):
    env = dict(os.environ, PYTHONPATH=str(ROOT))
    return subprocess.run([sys.executable, str(ROOT / 'sb_cli.py'), *args], capture_output=True, text=True,
                          env=env, cwd=cwd or ROOT)


def check_packaging_metadata():
    meta = tomllib.loads((ROOT / 'pyproject.toml').read_text())
    assert meta['project']['name'] == 'sbonest' and meta['project']['scripts'] == {'sbonest': 'sb_cli:main'}
    modules = meta['tool']['setuptools']['py-modules']
    missing = [m for m in modules if not (ROOT / f'{m}.py').exists()]
    assert not missing, missing
    for required in ('run', 'sbfit', 'sideband', 'sb_cli', 'sb_server', 'sb_import', 'sb_design', 'sb_compare'):
        assert required in modules, required
    module_name, function = meta['project']['scripts']['sbonest'].split(':')
    assert callable(getattr(importlib.import_module(module_name), function))
    requirements = (ROOT / 'requirements.txt').read_text().split() + (ROOT / 'requirements-sideband.txt').read_text().split()
    for dependency in meta['project']['dependencies']:
        base = dependency.split('>=')[0]
        assert any(req.startswith(base) for req in requirements), dependency


def check_commands():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-cli-'))
    assert run_cli('version').stdout.startswith('sbonest')
    assert run_cli('--help').returncode == 0
    for passthrough in (('import-bruker', '--help'), ('benchmark', '--help')):
        run = run_cli(*passthrough)
        assert run.returncode == 0 and 'usage' in run.stdout, passthrough
    demo = run_cli('init-demo', '--out', str(folder / 'demo'))
    assert demo.returncode == 0 and (folder / 'demo' / 'fit.json').exists()
    path = write_synthetic(folder / 'small', {'A1': (15., 285.), 'G2': (15., 285.), 'S3': (15., 285.)})
    config = _residue_config(json.loads(path.read_text()), ['A1', 'G2', 'S3'], {'A1'})
    config['init']['vary'] = ['kab', 'kba', 'v1n_scale', 'A1.R2b']
    config['Project Name'] = str(folder / 'small' / 'cli_fit')
    path.write_text(json.dumps(config))
    check = run_cli('check', str(path), '--no-pdf', '--identifiability')
    assert check.returncode == 0 and json.loads(check.stdout)['identifiability']['n_free'] == 4
    fit = run_cli('fit', str(path), '--no-pdf', '--workers', '1')
    assert fit.returncode == 0, fit.stderr[-1500:]
    result = json.loads((folder / 'small' / 'cli_fit_result.json').read_text())
    assert result['success']
    again = run_cli('fit', str(path), '--no-pdf')
    assert again.returncode == 1 and 'already exists' in again.stderr
    resume = run_cli('resume', str(path), '--no-pdf')
    assert resume.returncode == 0 and 'Resumed completed outputs' in resume.stdout
    report = run_cli('report', str(folder / 'small' / 'cli_fit_result.json'), '--out', str(folder / 'small' / 'cli_report'))
    assert report.returncode == 0 and (folder / 'small' / 'cli_report.pdf').exists()
    design = {'config': 'small/fit.json',
              'truth': {name: item['value'] for name, item in result['parameters'].items()},
              'scenarios': [{'name': 'one', 'datasets': [{'v1n_hz': 30., 'T': .03, 'sigma': .004,
                                                          'offsets_rel_ppm': {'min': -.3, 'max': .3, 'n': 7}}]}]}
    (folder / 'design.json').write_text(json.dumps(design))
    designed = run_cli('design', str(folder / 'design.json'), '--out', str(folder / 'design_out'))
    assert designed.returncode == 0, designed.stderr[-1500:]
    assert (folder / 'design_out' / 'design.txt').exists()
    bad = run_cli('fit', str(path), '--workers', '0')
    assert bad.returncode == 2


if __name__ == '__main__':
    check_packaging_metadata()
    check_commands()
    print('PASS: packaging metadata and unified command line')
