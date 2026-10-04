"""Retained wheel/sdist acceptance; requires test-only build and twine tooling."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import tempfile
import tomllib
import zipfile

ROOT = Path(__file__).resolve().parent
NUMERICAL_ENV = dict(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                     VECLIB_MAXIMUM_THREADS='1', MPLBACKEND='Agg',
                     PYTHONDONTWRITEBYTECODE='1', PYTHONNOUSERSITE='1',
                     PYTHONSAFEPATH='1')


def _hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _snapshot(path):
    return {str(p.relative_to(path)): _hash(p) for p in sorted(path.rglob('*')) if p.is_file()}


def check_distribution(out):
    """Build and exercise both artifact types, retaining commands and failed outputs."""
    out = Path(out).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, **NUMERICAL_ENV)
    env.pop('PYTHONPATH', None)
    env.pop('PYTHONHOME', None)
    assert 'PYTHONPATH' not in env
    manifest = {'root': str(ROOT), 'out': str(out), 'python': sys.executable,
                'environment': {**NUMERICAL_ENV, 'PYTHONPATH': None, 'PYTHONHOME': None},
                'commands': [], 'installations': {}, 'status': 'running',
                'cleanup': {'policy': 'Await every child; archive workspaces even on failure.',
                            'active_pid': None, 'workspaces': []}}

    def save():
        (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')

    def run(command, cwd, expected=0):
        assert 'PYTHONPATH' not in env
        index = len(manifest['commands'])
        record = {'argv': list(map(str, command)), 'cwd': str(cwd),
                  'expected': expected, 'exit': None, 'log': f'command-{index:02d}.log'}
        manifest['commands'].append(record)
        save()
        print(f'RUN {index}: {record["argv"]}', flush=True)
        # Register ownership before waiting; close the child on timeout/interruption.
        with (out / record['log']).open('w') as log:
            process = subprocess.Popen(record['argv'], cwd=cwd, env=env,
                                       stdout=log, stderr=subprocess.STDOUT,
                                       start_new_session=os.name == 'posix')
            manifest['cleanup']['active_pid'] = process.pid
            save()
            try:
                record['exit'] = process.wait(timeout=1200)
            finally:
                if process.poll() is None:
                    if os.name == 'posix':
                        os.killpg(process.pid, signal.SIGTERM)
                    else:
                        process.terminate()
                    try:
                        process.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        if os.name == 'posix':
                            os.killpg(process.pid, signal.SIGKILL)
                        else:
                            process.kill()
                        process.wait()
                record['exit'] = process.returncode
                manifest['cleanup']['active_pid'] = None
                save()
        assert (record['exit'] == 0 if expected == 0 else record['exit'] != 0), record
        print(f'EXIT {index}: {record["exit"]}', flush=True)
        return (out / record['log']).read_text()

    save()
    try:
        metadata = tomllib.loads((ROOT / 'pyproject.toml').read_text())
        modules = metadata['tool']['setuptools']['py-modules']
        packages = metadata['tool']['setuptools']['packages']
        version = metadata['project']['version']
        assert re.search(r'^version:\s*' + re.escape(version) + r'\s*$',
                         (ROOT / 'CITATION.cff').read_text(), re.M)
        released = re.findall(r'^## (\d+\.\d+\.\d+)\b', (ROOT / 'CHANGELOG.md').read_text(), re.M)
        assert released[0] == version
        manifest['version'] = version
        # Build outside the live source directory so setuptools cannot dirty it.
        stage = out / 'build-source'
        stage.mkdir()
        for name in [*(name + '.py' for name in modules), 'pyproject.toml', 'README.md',
                     'LICENSE', 'CITATION.cff', 'CHANGELOG.md']:
            shutil.copy2(ROOT / name, stage / name)
        for package in packages:
            shutil.copytree(ROOT / package, stage / package,
                            ignore=shutil.ignore_patterns('__pycache__'))
        manifest['source_sha256'] = _snapshot(stage)
        canonical = ROOT / 'example' / 'sideband_auto_H'
        canonical_config = json.loads((canonical / 'two_RF.json').read_text())
        data_hashes = {Path(name).name: _hash(canonical / name)
                       for name in canonical_config['datasets']}
        manifest['demo_data_sha256'] = data_hashes
        dist = out / 'dist'
        run([sys.executable, '-m', 'build', '--outdir', dist], stage)
        wheels = list(dist.glob('*.whl'))
        sdists = list(dist.glob('*.tar.gz'))
        assert len(wheels) == len(sdists) == 1
        run([sys.executable, '-m', 'twine', 'check', '--strict', *wheels, *sdists], out)
        manifest['artifacts'] = {p.name: _hash(p) for p in wheels + sdists}

        for kind, artifact in [('wheel', wheels[0]), ('sdist', sdists[0])]:
            base = out / kind
            base.mkdir()
            venv = base / 'venv'
            run([sys.executable, '-m', 'venv', venv], out)
            python = venv / ('Scripts/python.exe' if os.name == 'nt' else 'bin/python')
            cli = python.parent / ('sbonest.exe' if os.name == 'nt' else 'sbonest')
            install_artifact = artifact
            if kind == 'sdist':
                rebuilt = base / 'rebuilt'
                run([python, '-m', 'pip', 'wheel', '--no-deps', '--no-cache-dir',
                     '--wheel-dir', rebuilt, artifact], out)
                rebuilt_wheels = list(rebuilt.glob('*.whl'))
                assert len(rebuilt_wheels) == 1
                install_artifact = rebuilt_wheels[0]
            run([python, '-m', 'pip', 'install', '--no-cache-dir', install_artifact], out)
            run([python, '-m', 'pip', 'check'], out)
            run([python, '-m', 'pip', 'freeze'], out)
            workspace = Path(tempfile.mkdtemp(prefix=f'sbonest-{kind}-')).resolve()
            assert not workspace.is_relative_to(ROOT)
            cleanup = {'original': str(workspace), 'archived': str(base / 'workspace'), 'closed': False}
            manifest['cleanup']['workspaces'].append(cleanup)
            installation = {'venv': str(venv), 'python': str(python), 'cli': str(cli),
                            'cwd': str(workspace), 'artifact': str(install_artifact),
                            'artifact_sha256': _hash(install_artifact)}
            manifest['installations'][kind] = installation
            save()
            try:
                probe = '''import importlib, importlib.metadata, json, os, pathlib, sys
root, venv, names, expected_version = sys.argv[1:]
root, venv = pathlib.Path(root), pathlib.Path(venv)
assert 'PYTHONPATH' not in os.environ
assert not pathlib.Path.cwd().is_relative_to(root)
assert all(not pathlib.Path(p).resolve().is_relative_to(root) or
           pathlib.Path(p).resolve().is_relative_to(venv) for p in sys.path)
origins = {}
for name in json.loads(names):
    module = importlib.import_module(name)
    origin = pathlib.Path(module.__file__).resolve()
    assert origin.is_relative_to(venv), (name, str(origin))
    origins[name] = str(origin)
assert importlib.metadata.version('sbonest') == expected_version
print(json.dumps({'origins': origins, 'sys_path': sys.path,
                  'prefix': sys.prefix, 'cwd': str(pathlib.Path.cwd()),
                  'PYTHONPATH_present': 'PYTHONPATH' in os.environ}))
'''
                installation['isolation'] = json.loads(run(
                    [python, '-c', probe, ROOT, venv, json.dumps(modules + packages), version], workspace))
                with zipfile.ZipFile(install_artifact) as archive:
                    for name in modules:
                        assert hashlib.sha256(archive.read(name + '.py')).hexdigest() == manifest['source_sha256'][name + '.py']
                output = run([cli, 'version'], workspace)
                lines = output.strip().splitlines()
                assert lines[0] == f'sbonest {version}'
                source_hashes = dict(line.split(': ', 1) for line in lines[1:])
                assert 'sb_diagnostics.py' in source_hashes
                assert all(manifest['source_sha256'][name] == value for name, value in source_hashes.items())
                installation['version_source_sha256'] = source_hashes
                run([cli, 'init-demo', '--out', 'demo'], workspace)
                demo = workspace / 'demo'
                assert {name: _hash(demo / 'data' / name) for name in data_hashes} == data_hashes
                before = _snapshot(demo)
                run([cli, 'init-demo', '--out', 'demo'], workspace, expected=1)
                assert _snapshot(demo) == before
                installation['duplicate_demo_sha256'] = before
                checked = json.loads(run([cli, 'check', 'demo/fit.json', '--no-pdf'], workspace))
                assert checked['valid'] is True
                run([cli, 'fit', 'demo/fit.json', '--no-pdf', '--workers', '1'], workspace)
                result = json.loads((demo / 'fit_result.json').read_text())
                assert result['success'] is True
                assert result['n_points'] == 882 and result['n_parameters'] == 24
                assert result['provenance']['source_sha256'] == source_hashes
                before = _snapshot(demo)
                run([cli, 'fit', 'demo/fit.json', '--no-pdf', '--workers', '1'], workspace, expected=1)
                assert _snapshot(demo) == before
                installation['duplicate_fit_sha256'] = before
                run([cli, 'resume', 'demo/fit.json', '--no-pdf', '--workers', '1'], workspace)
                resumed = json.loads((demo / 'fit_result.json').read_text())
                assert resumed['success'] is True and resumed['chi2'] == result['chi2']
                run([cli, 'report', 'demo/fit_result.json', '--out', 'demo/report'], workspace)
                assert (demo / 'report_summary.txt').stat().st_size > 0
                assert (demo / 'report_summary.json').stat().st_size > 0
                assert (demo / 'report.pdf').stat().st_size > 0
                (workspace / 'malformed.json').write_text('{ malformed json\n')
                before = _snapshot(workspace)
                run([cli, 'check', 'malformed.json', '--no-pdf'], workspace, expected=1)
                assert _snapshot(workspace) == before
                installation['outputs_sha256'] = before
                installation['status'] = 'passed'
            finally:
                # Preserve calculation-time paths/hashes; do not rewrite checkpoints.
                shutil.move(str(workspace), base / 'workspace')
                cleanup['closed'] = True
                save()
        manifest['status'] = 'passed'
        print('PASS wheel and sdist installed workflows', flush=True)
    except BaseException as exc:
        manifest['status'] = 'failed'
        manifest['error'] = f'{type(exc).__name__}: {exc}'
        raise
    finally:
        save()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', required=True, help='New retained evidence directory')
    check_distribution(parser.parse_args().out)
