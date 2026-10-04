"""Retained canonical full-fit golden and exact same-machine worker checks."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]  # repository root: flat modules live there
sys.path.insert(0, str(ROOT))

import os

for _name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(_name, '1')
os.environ.setdefault('MPLBACKEND', 'Agg')

import argparse
import copy
import csv
import json
import tempfile

import numpy as np
from numpy.testing import assert_allclose

from run import load_config
from sbfit import run_config

REFERENCE = ROOT / 'results/auto_H_refit/fits/two_RF_result.json'


def check_historical(actual, reference):
    """Compare a new successful fit with the immutable historical fields."""
    assert actual['success'] is True, 'success'
    for key, value in (('n_points', 882), ('n_parameters', 24), ('dof', 858),
                       ('jacobian_rank', 24)):
        assert actual[key] == reference[key] == value, key
    names = list(reference['parameters'])
    assert list(actual['parameters']) == actual['parameter_order'] == names, 'parameter_order'
    assert_allclose(actual['chi2'], reference['chi2'], rtol=1e-7, atol=1e-7,
                    err_msg='chi2', equal_nan=False)
    for name in names:
        got, expected = actual['parameters'][name], reference['parameters'][name]
        assert got['vary'] is expected['vary'], f'{name}.vary'
        for field, rtol, atol in (
                ('value', 1e-6, 1e-5 if name in ('A1.R1H', 'S3.R1H') else 1e-8),
                ('stderr', 1e-5, 1e-7)):
            assert (got[field] is None) == (expected[field] is None), f'{name}.{field} null mask'
            if expected[field] is not None:
                assert_allclose(got[field], expected[field], rtol=rtol, atol=atol,
                                err_msg=f'{name}.{field}', equal_nan=False)


def check_exact(actual, expected, label='result'):
    """Compare all nested scientific values, including numeric dtype and bytes."""
    assert type(actual) is type(expected), f'{label}: type'
    if isinstance(actual, dict):
        assert list(actual) == list(expected), f'{label}: keys/order'
        for key in actual:
            check_exact(actual[key], expected[key], f'{label}.{key}')
    elif isinstance(actual, list):
        assert len(actual) == len(expected), f'{label}: shape'
        for index, (left, right) in enumerate(zip(actual, expected)):
            check_exact(left, right, f'{label}[{index}]')
    elif isinstance(actual, (int, float)) and not isinstance(actual, bool):
        left, right = np.asarray(actual), np.asarray(expected)
        assert left.dtype == right.dtype and left.shape == right.shape, f'{label}: dtype/shape'
        assert left.tobytes() == right.tobytes(), f'{label}: bytes'
    else:
        assert actual == expected, label


def check_workers(serial, parallel):
    """Exclude only provenance and output-prefix-derived paths from exact parity."""
    normalized = []
    for result in (serial, parallel):
        result = copy.deepcopy(result)
        result.pop('provenance')
        prefix = result['config']['Project Name']
        assert result['checkpoint'] == prefix + '_checkpoint', 'checkpoint path'
        result['checkpoint'] = '<output>_checkpoint'
        result['config']['Project Name'] = '<output>'
        normalized.append(result)
    check_exact(*normalized)
    for key in ('covariance',):
        left, right = (np.asarray(result[key], dtype=np.float64) for result in (serial, parallel))
        assert left.shape == right.shape == (24, 24), key
        assert left.dtype == right.dtype and left.tobytes() == right.tobytes(), key


def check_predictions(serial_path, parallel_path):
    """Compare prediction labels and every numeric CSV column exactly."""
    tables = []
    for path in (serial_path, parallel_path):
        with path.open(newline='', encoding='utf-8') as stream:
            reader = csv.DictReader(stream)
            tables.append((reader.fieldnames, list(reader)))
    columns, rows = tables[0]
    assert columns == tables[1][0], 'prediction columns'
    assert len(rows) == len(tables[1][1]) == 882, 'prediction count'
    for column in columns:
        if column == 'residue':
            assert [row[column] for row in rows] == [row[column] for row in tables[1][1]]
        else:
            left, right = (np.asarray([row[column] for row in table[1]], dtype=np.float64)
                           for table in tables)
            assert left.dtype == right.dtype and left.shape == right.shape, column
            assert left.tobytes() == right.tobytes(), column


def check_mutations(reference):
    """Prove comparator sensitivity on disposable fixtures before fitting."""
    fixture = copy.deepcopy(reference)
    fixture.update(success=True, parameter_order=list(reference['parameters']),
                   covariance=np.eye(24).tolist(), provenance={},
                   checkpoint=reference['config']['Project Name'] + '_checkpoint')
    # Identity covariance is only a comparator fixture, never a historical oracle.
    check_historical(fixture, reference)
    check_workers(fixture, fixture)
    results = {}
    for label in ('chi2', 'kab.value', 'kab.stderr', 'covariance'):
        changed = copy.deepcopy(fixture)
        if label == 'chi2':
            changed['chi2'] += 1
        elif label == 'covariance':
            changed['covariance'][0][0] += 1
        else:
            changed['parameters']['kab'][label.split('.')[1]] += 1
        try:
            if label == 'covariance':
                check_workers(changed, fixture)
            else:
                check_historical(changed, reference)
        except AssertionError as exc:
            assert label in str(exc), (label, str(exc))
            results[label] = {'rejected': True, 'assertion': str(exc)}
        else:
            raise AssertionError(f'Accepted mutation: {label}')
    check_historical(fixture, reference)
    check_workers(fixture, fixture)
    return results


def check_full_fit(folder):
    """Run the unchanged canonical configuration under both worker counts."""
    reference = json.loads(REFERENCE.read_text())
    mutations = check_mutations(reference)
    (folder / 'task-4-mutation.json').write_text(json.dumps(mutations, indent=2) + '\n')
    results = []
    for workers in (1, 2):
        config = load_config(ROOT / 'example/sideband_auto_H/two_RF.json')
        assert all(Path(path).is_absolute() for path in config['datasets'])
        config['Project Name'] = str(folder / f'workers-{workers}')
        run_config(config, no_pdf=True, workers=workers)
        results.append(json.loads((folder / f'workers-{workers}_result.json').read_text()))
    for result in results:
        check_historical(result, reference)
    check_workers(*results)
    check_predictions(folder / 'workers-1_predictions.csv', folder / 'workers-2_predictions.csv')
    comparison = {'success': True, 'workers': [1, 2], 'historical_chi2': reference['chi2'],
                  'chi2': [result['chi2'] for result in results], 'n_points': 882,
                  'n_parameters': 24, 'dof': 858, 'jacobian_rank': 24,
                  'parameter_order': results[0]['parameter_order'],
                  'historical_tolerances': {'chi2': [1e-7, 1e-7], 'value': [1e-6, 1e-8],
                                            'boundary_value_atol': 1e-5, 'stderr': [1e-5, 1e-7]},
                  'exact_worker_outputs': True, 'exact_prediction_numeric_bytes': True}
    (folder / 'comparison.json').write_text(json.dumps(comparison, indent=2) + '\n')


def main():
    """Run with a fresh optional evidence directory; retain artifacts on failure."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    if args.out is None:
        folder = Path(tempfile.mkdtemp(prefix='sbonest-golden-')).resolve()
    else:
        folder = args.out.expanduser().absolute()
        folder.mkdir(parents=True, exist_ok=False)
    print('Golden artifacts:', folder, flush=True)
    check_full_fit(folder)
    print('PASS: canonical historical golden and exact worker identity')


if __name__ == '__main__':
    main()
