"""Experimental design: expected local errors for planned acquisitions."""
# ruff: noqa: E402 -- Limit numerical libraries before importing them.
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]  # repository root: flat modules live there
sys.path.insert(0, str(ROOT))

import os

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(name, '1')
os.environ.setdefault('MPLBACKEND', 'Agg')

import json
import subprocess
import tempfile

import numpy as np

from sb_design import _offsets, load_design, run_design, select_measurements
from sbfit import SidebandModel
from test_sb_compare import write_synthetic



def truth_for(config_path):
    config = json.loads(Path(config_path).read_text())
    truth = {'kab': 15., 'kba': 285., 'v1n_scale': 1.0, **config['init']['initial']}
    return truth


def write_design(folder, config_path, truth, scenarios):
    Path(folder).mkdir(parents=True, exist_ok=True)
    design = {'config': os.path.relpath(config_path, folder), 'truth': truth, 'scenarios': scenarios}
    path = Path(folder) / 'design.json'
    path.write_text(json.dumps(design, indent=1) + '\n')
    return path


def check_validation():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-design-val-'))
    config_path = write_synthetic(folder / 'base', {'A1': (15., 285.), 'G2': (15., 285.), 'S3': (15., 285.)})
    truth = truth_for(config_path)
    good = [{'name': 'a', 'datasets': [{'v1n_hz': 25., 'T': .03, 'sigma': .004, 'offsets_rel_ppm': {'min': -.4, 'max': .4, 'n': 9}}]}]
    assert load_design(write_design(folder / 'ok', config_path, truth, good))['scenarios'][0]['name'] == 'a'
    cases = [
        ({**truth, 'extra': 1.}, good, 'unknown truth name'),
        ({k: v for k, v in truth.items() if k != 'kba'}, good, 'missing truth name'),
        (truth, [{'name': 'bad name', 'datasets': good[0]['datasets']}], 'scenario name'),
        (truth, [{'name': 'a', 'datasets': []}], 'empty datasets'),
        (truth, [{'name': 'a', 'datasets': [{'v1n_hz': 25., 'T': .03, 'sigma': .004}]}], 'offsets missing'),
        (truth, [{'name': 'a', 'datasets': [{'v1n_hz': 25., 'T': .03, 'sigma': 0., 'offsets_ppm': [1., 2.]}]}], 'zero sigma'),
        (truth, [{'name': 'a', 'datasets': [{'v1n_hz': 25., 'T': .03, 'sigma': .004, 'offsets_ppm': [1., 1.]}]}], 'repeated offsets'),
        (truth, [good[0], good[0]], 'duplicate names'),
    ]
    for index, (bad_truth, scenarios, label) in enumerate(cases):
        path = write_design(folder / f'bad{index}', config_path, bad_truth, scenarios)
        try:
            run_design(path, folder / f'bad{index}' / 'out')
        except ValueError:
            pass
        else:
            raise AssertionError(f'Accepted invalid design: {label}')
        assert not (folder / f'bad{index}' / 'out' / 'design.json').exists()
    assert _offsets({'min': 0., 'max': 1., 'n': 3}, 'x').tolist() == [0., .5, 1.]


def check_expected_errors_and_inputs():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-design-'))
    config_path = write_synthetic(folder / 'base', {'A1': (15., 285.), 'G2': (15., 285.), 'S3': (15., 285.)})
    truth = truth_for(config_path)
    grid = {'min': -.4, 'max': .4, 'n': 13}
    scenarios = [
        {'name': 'two_rf', 'datasets': [{'v1n_hz': 25., 'T': .03, 'sigma': .004, 'offsets_rel_ppm': grid},
                                        {'v1n_hz': 80., 'T': .03, 'sigma': .004, 'offsets_rel_ppm': grid}]},
        {'name': 'one_rf', 'datasets': [{'v1n_hz': 80., 'T': .03, 'sigma': .004, 'offsets_rel_ppm': grid}]},
        {'name': 'two_rf_noisier', 'datasets': [{'v1n_hz': 25., 'T': .03, 'sigma': .008, 'offsets_rel_ppm': grid},
                                                {'v1n_hz': 80., 'T': .03, 'sigma': .008, 'offsets_rel_ppm': grid}]},
    ]
    path = write_design(folder, config_path, truth, scenarios)
    paths = run_design(path, folder / 'out', workers=1)
    summary = json.loads(Path(paths['design_json']).read_text())
    rows = {row['name']: row for row in summary['scenarios']}
    assert rows['two_rf']['n_points'] == 78 and rows['one_rf']['n_points'] == 39
    se = {name: row['identifiability']['derived_se']['kex'] for name, row in rows.items()}
    assert se['two_rf'] is not None and se['two_rf'] < se['one_rf']
    # Doubling sigma doubles every expected standard error exactly (linear statistics).
    ratio = se['two_rf_noisier'] / se['two_rf']
    assert abs(ratio - 2.) < 1e-6, ratio
    assert rows['two_rf']['identifiability']['chi2_at_point'] < 1e-16, 'Noise-free inputs must reproduce the truth'
    # Written inputs are complete SBONEST datasets that reload to the same model.
    config = json.loads((folder / 'out' / 'scenarios' / 'two_rf' / 'design_config.json').read_text())
    model = SidebandModel(config, Path(config_path).parent)
    vector = np.array([truth[name] for name in model.parameter_names])
    assert float(np.sum(np.asarray(model.errFunc(vector)) ** 2)) < 1e-16
    text = Path(paths['design_txt']).read_text()
    assert 'two_rf | 78 | 2 |' in text and Path(paths['design_pdf']).stat().st_size > 1000
    try:
        run_design(path, folder / 'out')
    except FileExistsError:
        pass
    else:
        raise AssertionError('Existing design output was overwritten')
    # Pooled evaluation gives the same report.
    pooled = json.loads(Path(run_design(path, folder / 'out_pool', workers=2)['design_json']).read_text())
    assert [r['identifiability'] for r in pooled['scenarios']] == [r['identifiability'] for r in summary['scenarios']]


def check_truth_result_and_cli():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-design-cli-'))
    result = {'parameters': {name: {'value': value} for name, value in truth_for(
        write_synthetic(folder / 'base', {'A1': (15., 285.), 'G2': (15., 285.), 'S3': (15., 285.)})).items()}}
    (folder / 'truth_result.json').write_text(json.dumps(result))
    design = {'config': 'base/fit.json', 'truth_result': 'truth_result.json',
              'scenarios': [{'name': 's', 'datasets': [{'v1n_hz': 40., 'T': .03, 'sigma': .004,
                                                        'offsets_ppm': [119.6, 119.8, 120., 120.2, 120.4, 110., 115.]}]}]}
    (folder / 'design.json').write_text(json.dumps(design))
    env = dict(os.environ, PYTHONPATH=str(ROOT))
    run = subprocess.run([sys.executable, str(ROOT / 'sb_workflow.py'), 'design', str(folder / 'design.json'),
                          '--out', str(folder / 'cli_out')], capture_output=True, text=True, env=env, cwd=ROOT)
    assert run.returncode == 0, run.stderr[-2000:]
    summary = json.loads((folder / 'cli_out' / 'design.json').read_text())
    assert summary['scenarios'][0]['n_points'] == 21 and summary['truth']['kab'] == 15.


def check_optimization():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-design-opt-'))
    config_path = write_synthetic(folder / 'base', {'A1': (15., 285.), 'G2': (15., 285.), 'S3': (15., 285.)})
    truth = truth_for(config_path)
    grid = {'min': 105., 'max': 135., 'n': 31}
    candidate = {'name': 'cand', 'datasets': [{'v1n_hz': 25., 'T': .03, 'sigma': .004, 'offsets_ppm': grid},
                                              {'v1n_hz': 80., 'T': .03, 'sigma': .004, 'offsets_ppm': grid}]}
    design = {'config': os.path.relpath(config_path, folder), 'truth': truth, 'scenarios': [candidate],
              'optimize': {'scenario': 'cand', 'budget': 16, 'criterion': 'kex', 'min_per_dataset': 3}}
    (folder / 'design.json').write_text(json.dumps(design))
    paths = run_design(folder / 'design.json', folder / 'out')
    summary = json.loads(Path(paths['design_json']).read_text())
    rows = {row['name']: row for row in summary['scenarios']}
    opt = summary['optimization']
    assert set(rows) == {'cand', 'cand_optimized_16', 'cand_uniform_16'}
    selected = opt['selected_offsets_ppm']
    assert sum(len(v) for v in selected.values()) == 16 and all(len(v) >= 3 for v in selected.values())
    assert rows['cand_optimized_16']['n_points'] == 16 * 3
    # The incremental criterion equals the re-evaluated variance of kex for the optimized design.
    se = rows['cand_optimized_16']['identifiability']['derived_se']['kex']
    assert abs(opt['optimized_criterion'] - se**2) < 1e-6 * se**2, (opt['optimized_criterion'], se**2)
    assert abs(opt['candidate_criterion'] - rows['cand']['identifiability']['derived_se']['kex']**2) < 1e-6 * opt['candidate_criterion']
    assert se <= rows['cand_uniform_16']['identifiability']['derived_se']['kex'], 'Optimized design worse than uniform'
    assert opt['criterion_path'][0]['n_measurements'] == 62 and opt['criterion_path'][-1]['n_measurements'] == 16
    assert all(b['criterion'] >= a['criterion'] - 1e-9 for a, b in zip(opt['criterion_path'], opt['criterion_path'][1:]))
    assert 'Optimized design' in Path(paths['design_txt']).read_text()
    # Direct selection with the D criterion and a parameter criterion on a tiny Jacobian.
    rng = np.random.default_rng(0)
    J = rng.normal(size=(12, 3))
    groups = {(0, k): [k, k + 6] for k in range(6)}
    kept, path = select_measurements(J, groups, 4, None, min_per_dataset=0)
    assert len(kept) == 4 and path[-1]['n_measurements'] == 4 and path[-1]['criterion'] >= path[0]['criterion']
    kept, _ = select_measurements(J, groups, 5, np.array([0., 0., 1.]), min_per_dataset=0)
    assert len(kept) == 5
    for bad in ({'scenario': 'cand', 'budget': 1}, {'scenario': 'missing', 'budget': 10},
                {'scenario': 'cand', 'budget': 4, 'min_per_dataset': 3}, {'scenario': 'cand', 'budget': 10, 'criterion': ''}):
        design['optimize'] = bad
        (folder / 'bad.json').write_text(json.dumps(design))
        try:
            run_design(folder / 'bad.json', folder / 'bad_out')
        except ValueError:
            pass
        else:
            raise AssertionError(f'Invalid optimize settings accepted: {bad}')
    relative = dict(design, scenarios=[{'name': 'cand', 'datasets': [{'v1n_hz': 25., 'T': .03, 'sigma': .004,
                                                                     'offsets_rel_ppm': {'min': -.4, 'max': .4, 'n': 9}}]}],
                    optimize={'scenario': 'cand', 'budget': 4})
    (folder / 'rel.json').write_text(json.dumps(relative))
    try:
        run_design(folder / 'rel.json', folder / 'rel_out')
    except ValueError as exc:
        assert 'offsets_ppm' in str(exc)
    else:
        raise AssertionError('Relative candidate grid accepted for optimization')


if __name__ == '__main__':
    check_validation()
    check_expected_errors_and_inputs()
    check_truth_result_and_cli()
    check_optimization()
    print('PASS: design validation, expected errors, synthetic inputs, CLI and optimization')
