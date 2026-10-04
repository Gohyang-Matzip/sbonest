"""Independent coverage arithmetic and a real Sideband study CLI check."""
import copy
import importlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile

from numpy.testing import assert_allclose

from test_sb_output import small_config


ROOT = Path(__file__).resolve().parent


def runner():
    assert importlib.util.find_spec('validate_uncertainty') is not None, 'Missing synthetic coverage runner'
    return importlib.import_module('validate_uncertainty')


def check_coverage_arithmetic():
    # Two eligible intervals, exactly one covered; one missing SE, one zero SE,
    # and one failed fit must not enter the coverage denominator.
    sample = {'success': True, 'parameters': [2., 8., 1.], 'stderr': [.1, 0., 0.],
              'derived_se': {'kex': .1, 'pB': .01}, 'kex': 10., 'pB': .2,
              'at_bounds': ['kab']}
    outside = {**sample, 'parameters': [3., 8., 1.], 'kex': 11., 'pB': 3./11,
               'at_bounds': []}
    missing = {**sample, 'stderr': [None, 0., 0.],
               'derived_se': {'kex': None, 'pB': None}, 'at_bounds': []}
    zero = {**sample, 'stderr': [0., 0., 0.],
            'derived_se': {'kex': 0., 'pB': 0.}, 'at_bounds': []}
    failed = {'success': False, 'parameters': None, 'stderr': None, 'at_bounds': []}
    rows = [sample, outside, missing, zero, failed]
    truth = {'kab': 2., 'kba': 8., 'R1': 1.}
    summary = runner().coverage_summary(rows, list(truth), ['kab'], truth)
    assert summary['successful'] == 4 and summary['failed'] == 1
    for name in ('kab', 'kex', 'pB'):
        result = summary['parameters'][name]
        assert result['eligible_se'] == 2 and result['covered'] == 1
        assert result['finite_se'] == 3 and result['zero_se'] == 1
        assert result['unavailable_se'] == 1 and result['fraction'] == .5
        assert_allclose(result['binomial_se'], .3535533905932738)
        assert_allclose(result['wilson_95'], [.0945312057342307, .9054687942657693])
        assert result['at_bounds'] == 1 and result['boundary_fraction'] == .25
    assert summary['parameters']['kba']['fixed']
    assert summary['parameters']['R1']['fixed']
    assert summary['parameters']['kba']['fraction'] is None
    excluded = runner().coverage_summary(rows, list(truth), ['R1'], truth)
    assert excluded['parameters']['kex']['fixed'] and excluded['parameters']['pB']['fixed']
    unavailable = runner().coverage_summary([missing, failed], list(truth), ['kab'], truth)
    assert unavailable['parameters']['kab']['fraction'] is None
    assert unavailable['parameters']['kab']['wilson_95'] is None
    json.dumps(summary, allow_nan=False)


def study_inputs(folder):
    cfg = small_config(folder / 'unused-fit-output')
    truth = {'kab': 15., 'kba': 285., 'v1n_scale': 1.08, 'A1.peak_ppm': 120.,
             'A1.dw_ppm': 3., 'A1.R1': 1.5, 'A1.R2a': 12., 'A1.R2b': 15.,
             'A1.R1H': 2., 'A1.R2H': 25.}
    config_path, truth_path = folder / 'fit.json', folder / 'truth.json'
    config_path.write_text(json.dumps(cfg))
    truth_path.write_text(json.dumps({'truth': truth}))
    return config_path, truth_path, truth


def check_truth_validation_before_writes():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-coverage-truth-'))
    config_path, truth_path, truth = study_inputs(folder)
    for change in ({'A1.R1': 1.6}, {'v1n_scale': 2.}, {'kab': float('nan')}):
        invalid = {**truth, **change}
        truth_path.write_text(json.dumps(invalid))
        try:
            runner().run_study(config_path, truth_path, 2, 7, folder / 'study')
        except ValueError:
            pass
        else:
            raise AssertionError(f'Invalid generating truth accepted: {change}')
        assert not (folder / 'study').exists()
    truth_path.write_text(json.dumps({key: value for key, value in truth.items() if key != 'kab'}))
    try:
        runner().run_study(config_path, truth_path, 2, 7, folder / 'study')
    except ValueError:
        pass
    else:
        raise AssertionError('Incomplete truth accepted')
    assert not (folder / 'study').exists()


def check_real_cli_and_seed():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-coverage-cli-'))
    config_path, truth_path, truth = study_inputs(folder)
    out = folder / 'study'
    command = [sys.executable, str(ROOT / 'validate_uncertainty.py'), '--config', str(config_path),
               '--truth', str(truth_path), '--replicates', '3', '--seed', '812', '--out', str(out)]
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    assert result.returncode == 0, result.stderr
    summary = json.loads((out / 'coverage.json').read_text())
    samples = json.loads((out / 'samples.json').read_text())
    assert summary['successful'] == 3 and summary['failed'] == 0
    assert summary['parameters']['v1n_scale']['eligible_se'] == 3
    assert summary['parameters']['kex']['fixed'] and summary['parameters']['pB']['fixed']
    assert len(samples) == 3 and len(list((out / 'samples').glob('*.json'))) == 3
    settings = json.loads((out / 'study.json').read_text())
    assert settings['truth'] == truth and settings['settings']['seed'] == 812
    assert len(settings['provenance']['datasets']) == 2
    assert 'validate_uncertainty.py' in settings['provenance']['source_sha256']
    assert not (folder / 'unused-fit-output').exists()
    before = {str(p.relative_to(out)): p.read_bytes() for p in out.rglob('*') if p.is_file()}
    blocked = subprocess.run(command, text=True, capture_output=True, check=False)
    assert blocked.returncode != 0
    assert before == {str(p.relative_to(out)): p.read_bytes() for p in out.rglob('*') if p.is_file()}
    repeated = runner().run_study(config_path, truth_path, 3, 812, folder / 'repeat')
    assert repeated['parameters'] == summary['parameters']
    assert json.loads((folder / 'repeat' / 'samples.json').read_text()) == samples
    changed = copy.deepcopy(truth)
    changed['v1n_scale'] = 1.
    truth_path.write_text(json.dumps(changed))
    runner().run_study(config_path, truth_path, 1, 812, folder / 'different-truth')
    other = json.loads((folder / 'different-truth' / 'samples.json').read_text())[0]
    assert abs(other['parameters'][2] - samples[0]['parameters'][2]) > .05
    print('Synthetic coverage artifacts:', out)


if __name__ == '__main__':
    for check in (check_coverage_arithmetic, check_truth_validation_before_writes, check_real_cli_and_seed):
        check()
        print('PASS', check.__name__)
