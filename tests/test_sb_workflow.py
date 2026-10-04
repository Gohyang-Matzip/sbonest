"""Saved-result reports reconcile CSV evidence and never refit or overwrite."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]  # repository root: flat modules live there
sys.path.insert(0, str(ROOT))

import copy
import csv
from importlib import resources
import json
import subprocess
import tempfile
from unittest.mock import patch

import numpy as np
from numpy.testing import assert_allclose



def cli(*args):
    return subprocess.run([sys.executable, str(ROOT / 'sb_workflow.py'), *map(str, args)],
                          capture_output=True, text=True, timeout=120)


def check_demo():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-workflow-'))
    out = folder / 'new-demo'
    before = {p: p.read_bytes() for p in (ROOT / 'example/sideband_auto_H').glob('*')}
    result = cli('init-demo', '--out', out)
    assert result.returncode == 0, result.stderr
    cfg = json.loads((out / 'fit.json').read_text())
    assert cfg['Project Name'] == str(out / 'fit')
    original = json.loads((ROOT / 'example/sideband_auto_H/two_RF.json').read_text())
    original['Project Name'] = str(out / 'fit')
    original['datasets'] = ['data/noisy_25.txt', 'data/noisy_100.txt']
    assert cfg == original
    assert len(cfg['datasets']) == 2
    for name, source in zip(cfg['datasets'], ('noisy_25.txt', 'noisy_100.txt')):
        assert not Path(name).is_absolute()
        assert (out / name).read_bytes() == (ROOT / 'results/peakwise_H_fit' / source).read_bytes()
    assert before == {p: p.read_bytes() for p in before}
    saved = {p.relative_to(out): p.read_bytes() for p in out.rglob('*') if p.is_file()}
    assert cli('init-demo', '--out', out).returncode != 0
    assert saved == {p.relative_to(out): p.read_bytes() for p in out.rglob('*') if p.is_file()}
    missing = folder / 'missing'
    dangling = folder / 'dangling'
    dangling.symlink_to(missing, target_is_directory=True)
    assert cli('init-demo', '--out', dangling).returncode != 0
    assert dangling.is_symlink() and not missing.exists()
    print('Demo artifacts:', out)


def check_demo_resources():
    from sb_workflow import init_demo
    bundled = resources.files('sbonest_data').joinpath('sideband_auto_H')
    names = {'two_RF.json', 'noisy_25.txt', 'noisy_100.txt'}
    assert {entry.name for entry in bundled.iterdir()} == names
    canonical = json.loads((ROOT / 'example/sideband_auto_H/two_RF.json').read_text())
    canonical['datasets'] = ['noisy_25.txt', 'noisy_100.txt']
    assert json.loads(bundled.joinpath('two_RF.json').read_bytes()) == canonical
    for name in canonical['datasets']:
        assert bundled.joinpath(name).read_bytes() == (ROOT / 'results/peakwise_H_fit' / name).read_bytes()
    folder = Path(tempfile.mkdtemp(prefix='sbonest-demo-resources-'))
    with patch('sb_workflow.__file__', str(folder / 'no-checkout/sb_workflow.py')):
        assert init_demo(folder / 'portable').is_file()
    read_bytes = type(bundled).read_bytes
    for name in ('two_RF.json', 'noisy_25.txt', 'noisy_100.txt'):
        def fail_read(path):
            if path.name == name:
                raise OSError('Bundled input unavailable')
            return read_bytes(path)
        out = folder / name / 'must-not-exist'
        with patch.object(type(bundled), 'read_bytes', fail_read):
            try:
                init_demo(out)
            except OSError:
                pass
            else:
                raise AssertionError(f'Accepted unreadable resource {name}')
        assert not out.parent.exists()
    print('Resource boundary artifacts:', folder)


def fixture(folder):
    info = {'success': True, 'method': 'Sideband', 'n_points': 4, 'n_parameters': 1,
            'dof': 3, 'chi2': 6., 'kex': 300., 'pB': .05, 'v1n_hz': [25., 100.],
            'config': {'datasets': ['a.txt', 'b.txt'],
                       'residues': [{'name': x, 'flag': 'on'} for x in ('A1', 'G2')]},
            'parameters': {'kab': {'value': 15., 'stderr': .5, 'vary': True}},
            'multistart': [{'index': 0, 'success': False, 'chi2': 9., 'selected': False,
                            'status': 0, 'message': 'evaluation limit'},
                           {'index': 1, 'success': True, 'chi2': 6., 'selected': True,
                            'status': 1, 'message': 'converged'}],
            'profiles': {'base_chi2': 6., 'warnings': ['synthetic profile'],
                         'v1n_scale': [{'target': .9, 'success': True, 'chi2': 5.,
                                        'delta_chi2': -1., 'message': 'lower'},
                                       {'target': 1., 'success': False, 'chi2': None,
                                        'delta_chi2': None, 'message': 'failed'},
                                       {'target': 1.1, 'success': True, 'chi2': 7.,
                                        'delta_chi2': 1., 'message': 'converged'}]},
            'bootstrap': {'replicates': 3, 'seed': 42, 'confidence': .95, 'successful': 2,
                          'parameter_names': ['kab'], 'varying_parameters': ['kab'],
                          'samples': [{'index': 0, 'success': True, 'parameters': [14.],
                                       'stderr': [.4], 'chi2': 4., 'message': 'converged'},
                                      {'index': 1, 'success': False, 'parameters': None,
                                       'stderr': None, 'chi2': None, 'message': 'failed'},
                                      {'index': 2, 'success': True, 'parameters': [16.],
                                       'stderr': [.6], 'chi2': 7., 'message': 'converged'}],
                          'intervals': {'kab': {'lower': 14.05, 'median': 15., 'upper': 15.95,
                                               'n_success': 2, 'fixed': False}},
                          'warnings': ['Synthetic bootstrap: too few replicates.']},
            'warnings': ['Hand-calculated synthetic report fixture.']}
    rows = []
    for residue, ds, residual in [('A1', 0, -1.), ('A1', 1, 1.), ('G2', 0, 2.), ('G2', 1, 0.)]:
        rows.append({'residue': residue, 'dataset_index': ds, 'field_mhz': 60.,
                     'saturation_s': .2, 'v1n_nominal_hz': [25., 100.][ds],
                     'v1n_hz': [25., 100.][ds], 'offset_ppm': 120.,
                     'observed': .5 + residual * .1, 'predicted': .5, 'sigma': .1,
                     'residual_sigma': residual})
    result = folder / 'fit_result.json'
    predictions = folder / 'fit_predictions.csv'
    result.write_text(json.dumps(info))
    write_rows(predictions, rows)
    return result, predictions, info, rows


def write_rows(path, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def check_report():
    from sb_report import regenerate_report
    folder = Path(tempfile.mkdtemp(prefix='sbonest-report-'))
    result, predictions, info, _ = fixture(folder)
    before = result.read_bytes(), predictions.read_bytes()
    with patch('sbfit.SidebandModel.fit', side_effect=AssertionError('Report must not refit')):
        paths = regenerate_report(result, folder / 'report')
    summary = json.loads(Path(paths['summary_json']).read_text())
    assert_allclose(summary['overall']['chi2'], 6.)
    assert_allclose(summary['overall']['residual_mean'], .5)
    assert_allclose(summary['overall']['residual_rms'], np.sqrt(1.5))
    assert summary['overall']['n_points'] == 4
    assert summary['by_residue']['A1'] == {'n_points': 2, 'chi2': 2., 'residual_mean': 0., 'residual_rms': 1.}
    assert_allclose(summary['by_dataset']['0']['chi2'], 5.)
    assert len(summary['by_residue_dataset']) == 4
    assert summary['parameters'] == info['parameters']
    assert summary['multistart'] == info['multistart']
    assert summary['profiles'] == info['profiles']
    assert summary['bootstrap'] == info['bootstrap']
    assert Path(paths['pdf']).read_bytes().startswith(b'%PDF')
    assert b'/Type /Page' in Path(paths['pdf']).read_bytes()
    assert '-1' in Path(paths['summary_txt']).read_text()
    assert before == (result.read_bytes(), predictions.read_bytes())
    assert cli('report', result, '--out', folder / 'cli_report').returncode == 0
    print('Report artifacts:', folder)


def check_invalid_report_no_writes():
    from sb_report import regenerate_report
    folder = Path(tempfile.mkdtemp(prefix='sbonest-report-invalid-'))
    result, predictions, original, rows = fixture(folder)
    mutations = [lambda info, data: info.update(chi2=7),
                 lambda info, data: info.update(n_points=5),
                 lambda info, data: data[0].update(residue='wrong'),
                 lambda info, data: data[0].update(dataset_index=2),
                 lambda info, data: data[0].update(sigma=0),
                 lambda info, data: data[0].update(observed=float('nan')),
                 lambda info, data: data[0].update(residual_sigma=1),
                 lambda info, data: data[0].update(v1n_hz=24),
                 lambda info, data: data[0].pop('field_mhz'),
                 lambda info, data: info.update(success=False),
                 lambda info, data: info['parameters']['kab'].update(stderr=-1),
                 lambda info, data: info['profiles']['v1n_scale'][0].update(delta_chi2=1),
                 lambda info, data: info['multistart'][1].update(index=0),
                 lambda info, data: info['bootstrap'].update(successful=3),
                 lambda info, data: info['bootstrap']['samples'][0].update(parameters=[None]),
                 lambda info, data: info['bootstrap']['intervals']['kab'].update(lower=20)]
    for index, mutation in enumerate(mutations):
        info, data = copy.deepcopy(original), copy.deepcopy(rows)
        mutation(info, data)
        result.write_text(json.dumps(info))
        # Varying CSV keys deliberately exercises missing-column validation.
        with predictions.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(data[0]), extrasaction='ignore')
            writer.writeheader()
            writer.writerows(data)
        out = folder / f'invalid_{index}' / 'report'
        try:
            regenerate_report(result, out)
        except ValueError:
            pass
        else:
            raise AssertionError(f'Accepted invalid saved result {index}')
        assert not out.parent.exists()
    result, predictions, _, _ = fixture(folder)
    for suffix in ('_summary.json', '_summary.txt', '.pdf'):
        prefix = folder / ('blocked' + str(len(list(folder.iterdir()))))
        target = folder / (prefix.name + '-must-not-exist')
        protected = Path(str(prefix) + suffix)
        protected.symlink_to(target)
        before = set(folder.iterdir())
        assert cli('report', result, '--out', prefix).returncode != 0
        assert before == set(folder.iterdir())
        assert protected.is_symlink() and not target.exists()


def check_saved_fit_report():
    from sb_report import regenerate_report
    from sbfit import run_config
    from test_sb_output import small_config
    folder = Path(tempfile.mkdtemp(prefix='sbonest-saved-fit-report-'))
    run_config(small_config(folder), no_pdf=True)
    result = folder / 'fit_result.json'
    with patch('sbfit.SidebandModel.fit', side_effect=AssertionError('Report must not refit')):
        paths = regenerate_report(result, folder / 'saved')
    original = json.loads(result.read_text())
    summary = json.loads(Path(paths['summary_json']).read_text())
    assert summary['overall']['n_points'] == 294
    assert_allclose(summary['overall']['chi2'], original['chi2'], rtol=1e-12)
    assert_allclose(sum(row['chi2'] for row in summary['by_dataset'].values()), original['chi2'], rtol=1e-12)


if __name__ == '__main__':
    for check in (check_demo, check_demo_resources, check_report, check_invalid_report_no_writes, check_saved_fit_report):
        check()
        print('PASS', check.__name__)
