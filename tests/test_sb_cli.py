"""Unified command line and packaging metadata."""
# ruff: noqa: E402 -- Limit numerical libraries before importing them.
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]  # repository root: flat modules live there
sys.path.insert(0, str(ROOT))

import os

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(name, '1')
os.environ.setdefault('MPLBACKEND', 'Agg')

import importlib
import contextlib
import io
import json
import subprocess
import tempfile
import tomllib
from unittest.mock import patch

from sb_compare import _residue_config
from test_sb_compare import write_synthetic



def run_cli(*args, cwd=None):
    env = dict(os.environ, PYTHONPATH=str(ROOT))
    return subprocess.run([sys.executable, str(ROOT / 'sb_cli.py'), *args], capture_output=True, text=True,
                          env=env, cwd=cwd or ROOT, timeout=120)


def check_version_without_resources():
    """Version hashes match provenance without reading configuration or datasets."""
    import sb_cli
    from sb_report import provenance
    from importlib.metadata import PackageNotFoundError

    folder = Path(tempfile.mkdtemp(prefix='sbonest-cli-version-'))
    path = write_synthetic(folder, {'A1': (15., 285.), 'G2': (15., 285.), 'S3': (15., 285.)})
    expected = provenance(json.loads(path.read_text()), folder)['source_sha256']
    original = Path.read_bytes

    def source_only(path):
        assert path.parent == ROOT and path.name in expected, path
        return original(path)

    for installed in (True, False):
        output = io.StringIO()
        with patch.object(Path, 'read_text', side_effect=AssertionError('resource access')), \
                patch.object(Path, 'read_bytes', source_only), \
                patch('sb_report.provenance', side_effect=AssertionError('provenance access')), \
                patch('importlib.metadata.version', return_value='1.2.0',
                      side_effect=None if installed else PackageNotFoundError), \
                contextlib.redirect_stdout(output):
            assert sb_cli.main(['version']) == 0
        lines = output.getvalue().splitlines()
        assert lines[0] == ('sbonest 1.2.0' if installed else
                            'sbonest (not installed as a package; running from source)')
        hashes = dict(line.split(': ', 1) for line in lines[1:])
        assert list(hashes) == list(expected)
        assert hashes == expected


def check_compare_and_serve_contracts():
    """Preserve script argument conversion, delegation, defaults, and rejection."""
    import sb_cli

    config = {'init': {'Method': 'Sideband'}}
    with patch('sb_cli._config', return_value=(config, ROOT)), \
            patch('sb_compare.compare_models', return_value={}) as models, \
            patch('sb_compare.run_comparison', return_value={}) as residues:
        assert sb_cli.main(['compare', 'config.json', '--out', 'out', '--models',
                            'Sideband', 'Sideband_3st_Linear', '--h-ppm-c', 'A1=7.1',
                            'G2=8.2', '--pdf', '--workers', '2']) == 0
        models.assert_called_once_with(config, ROOT, 'out',
                                       models=('Sideband', 'Sideband_3st_Linear'),
                                       h_ppm_c={'A1': 7.1, 'G2': 8.2}, no_pdf=False, workers=2)
        assert sb_cli.main(['compare', 'config.json', '--out', 'out', '--models', 'Sideband']) == 0
        assert models.call_args.kwargs['h_ppm_c'] is None
        assert sb_cli.main(['compare', 'config.json', '--out', 'out']) == 0
        residues.assert_called_once_with(config, ROOT, 'out', no_pdf=True, workers=1)
        models.reset_mock()
        for shift in ('A1=bad', 'A1', '=7.1', 'A1='):
            with contextlib.redirect_stderr(io.StringIO()):
                try:
                    sb_cli.main(['compare', 'config.json', '--out', 'out', '--models',
                                 'Sideband', '--h-ppm-c', shift])
                except SystemExit as exc:
                    assert exc.code == 1
                else:
                    raise AssertionError(shift)
        models.assert_not_called()
    with patch('sb_server.main', return_value=0) as serve:
        assert sb_cli.main(['serve']) == 0
        serve.assert_called_with(['--host', '127.0.0.1', '--port', '5050'])
        assert sb_cli.main(['serve', '--host', 'localhost', '--port', '5057',
                            '--token', 'qa-token', '--max-age-days', '30']) == 0
        serve.assert_called_with(['--host', 'localhost', '--port', '5057',
                                  '--token', 'qa-token', '--max-age-days', '30.0'])
    with patch('sb_server.app.run') as run, patch('sb_server.archive_expired') as archive, \
            contextlib.redirect_stderr(io.StringIO()):
        try:
            sb_cli.main(['serve', '--max-age-days', '0'])
        except SystemExit as exc:
            assert exc.code == 2
        else:
            raise AssertionError('invalid age accepted')
        run.assert_not_called()
        archive.assert_not_called()


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
    # The committed API reference matches the generated text, and the run workflow
    # names stay importable from sbfit after the split into sb_run.
    check = subprocess.run([sys.executable, str(ROOT / 'scripts/generate_api_reference.py'), '--check'],
                           capture_output=True, text=True, cwd=ROOT)
    assert check.returncode == 0, check.stderr
    import sbfit
    import sb_run
    assert sbfit.run_config is sb_run.run_config and sbfit.check_config is sb_run.check_config
    assert 'sb_run' in modules


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
    check_version_without_resources()
    check_compare_and_serve_contracts()
    check_commands()
    print('PASS: packaging metadata and unified command line')
