"""Executable preflight checks: validate before fitting and leave no outputs."""

from pathlib import Path
import sys
from _env import ROOT

import copy
import json
import subprocess
import tempfile
from unittest.mock import patch

from numpy.testing import assert_allclose

import sbfit
from test_sb_output import small_config


def check_optional_validation_before_fit():
    # A malformed optional scan must fail before either the grid or optimizer.
    folder = Path(tempfile.mkdtemp(prefix='sbonest-check-invalid-'))
    for name, settings in (
        ('profile', {'v1n_scale': []}),
        ('profile', {'unknown': [1.]}),
        ('multistart', {'random_starts': 1}),
        ('multistart', {'starts': [{'not_a_parameter': 1.}]}),
        ('bootstrap', {'replicates': 3}),
        ('bootstrap', {'replicates': True, 'seed': 1}),
        ('bootstrap', {'replicates': 3, 'seed': 1, 'confidence': 1.}),
    ):
        cfg = small_config(folder / name)
        cfg['init'][name] = settings
        with patch('sbfit.least_squares', side_effect=AssertionError('Optimizer ran')), \
                patch('sbfit.generate_initial_parameters', side_effect=AssertionError('Grid ran')):
            try:
                sbfit.run_config(cfg, no_pdf=True)
            except ValueError:
                pass
            else:
                raise AssertionError(f'Invalid {name} accepted')
        assert not list(folder.iterdir()), 'Invalid settings created output directories'


def check_preflight_summary_and_no_outputs():
    assert hasattr(sbfit, 'check_config'), 'Missing zero-fit preflight API'
    folder = Path(tempfile.mkdtemp(prefix='sbonest-check-'))
    cfg = small_config(folder / 'new-output')
    cfg['init']['bootstrap'] = {'replicates': 3, 'seed': 17}
    cfg['init']['profile'] = {'v1n_scale': [1., 1.1]}
    before = copy.deepcopy(cfg)
    with patch('sbfit.least_squares', side_effect=AssertionError('Optimizer ran')), \
            patch.object(sbfit.SidebandModel, 'fit', side_effect=AssertionError('Fit ran')):
        summary = sbfit.check_config(cfg)
    assert summary['valid'], summary['errors']
    assert summary['n_points'] == 294 and summary['n_parameters'] == 1
    assert summary['free_parameters'] == ['v1n_scale']
    assert 'A1.R1H' in summary['fixed_parameters']
    assert len(summary['datasets']) == 2
    assert [row['n_points'] for row in summary['datasets']] == [147, 147]
    assert [row['nominal_v1n_hz'] for row in summary['datasets']] == [25., 100.]
    assert all(Path(row['path']).is_absolute() for row in summary['datasets'])
    assert cfg == before and not list(folder.iterdir())
    json.dumps(summary, allow_nan=False)
    Path(cfg['Project Name'] + '_result.json').parent.mkdir()
    conflict = Path(cfg['Project Name'] + '_result.json')
    conflict.symlink_to(folder / 'absent.json')
    blocked = sbfit.check_config(cfg)
    assert not blocked['valid'] and str(conflict) in blocked['output_conflicts']
    assert not (folder / 'absent.json').exists()


def check_prepare_freshness_and_validation():
    assert hasattr(sbfit.SidebandModel, 'prepare_fit'), 'Missing preparation API'
    cfg = small_config(Path(tempfile.mkdtemp(prefix='sbonest-prepare-')))
    model = sbfit.SidebandModel(cfg)
    with patch('sbfit.least_squares', side_effect=AssertionError('Optimizer ran')):
        initial = model.prepare_fit()
    assert_allclose(initial[:3], [15., 285., 1.])
    assert_allclose(initial[3:], [120., 3., 1.5, 12., 15., 2., 25.])
    assert model.free.tolist() == [2] and model._fit_data is None
    model.result = object()
    model.rank = 99
    model._fit_data = model._prepare_data()
    cfg['init']['initial']['A1.R1'] = 3.
    cfg['init']['bounds'] = {'v1n_scale': [1.05, 1.15]}
    try:
        model.prepare_fit()
    except ValueError:
        pass
    else:
        raise AssertionError('Out-of-bounds baseline accepted')
    assert model.result is None and model._fit_data is None
    assert not hasattr(model, 'rank'), 'Old rank survived a new failed preparation'
    assert model.initial_parameters[5] == 3.
    assert_allclose(model.lower[2], 1.05)
    assert_allclose(model.upper[2], 1.15)
    # An explicit restart must still recover from an invalid default baseline.
    cfg['init']['multistart'] = {'starts': [{'v1n_scale': 1.10}]}
    p, _ = sbfit.run_config(cfg, no_pdf=True)
    assert p[5] == 3. and 1.05 <= p[2] <= 1.15


def check_invalid_names_bounds_and_cli():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-check-cli-'))
    cfg = small_config(folder / 'outputs')
    for changes in (
        {'initial': {'missing': 1.}},
        {'bounds': {'v1n_scale': [1.2, .8]}},
        {'vary': ['v1n_scale', 'v1n_scale']},
        {'vary': ['missing']},
        {'max_nfev': 0},
        {'kex': {'min': 300., 'max': 300., 'nsteps': 0}},
    ):
        invalid = copy.deepcopy(cfg)
        invalid['init'].update(changes)
        with patch('sbfit.generate_initial_parameters', side_effect=AssertionError('Grid ran')):
            summary = sbfit.check_config(invalid)
        assert not summary['valid'] and summary['errors'], changes
    config_path = folder / 'config.json'
    config_path.write_text(json.dumps(cfg))
    for script in ('run.py', 'sbfit.py'):
        result = subprocess.run([sys.executable, str(ROOT / script), str(config_path), '--check'],
                                text=True, capture_output=True, cwd=folder, check=False)
        assert result.returncode == 0, result.stderr
        assert json.loads(result.stdout)['valid']
        assert not (folder / 'outputs').exists()
    cfg['init']['profile'] = {'v1n_scale': []}
    config_path.write_text(json.dumps(cfg))
    result = subprocess.run([sys.executable, str(ROOT / 'run.py'), str(config_path), '--check'],
                            text=True, capture_output=True, cwd=folder, check=False)
    assert result.returncode == 1
    assert not json.loads(result.stdout)['valid']


if __name__ == '__main__':
    for check in (check_optional_validation_before_fit, check_preflight_summary_and_no_outputs,
                  check_prepare_freshness_and_validation, check_invalid_names_bounds_and_cli):
        check()
        print('PASS', check.__name__)
