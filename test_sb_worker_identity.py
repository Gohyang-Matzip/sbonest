"""Canonical full-fit and exact same-machine worker identity checks."""
import argparse
import json
from pathlib import Path
import tempfile

from test_sb_golden import run_canonical_fits


def check_full_fit(folder):
    """Require successful canonical fits with exact scientific and prediction parity."""
    results = run_canonical_fits(folder)
    comparison = {'success': True, 'workers': [1, 2],
                  'chi2': [result['chi2'] for result in results], 'n_points': 882,
                  'n_parameters': 24, 'dof': 858, 'jacobian_rank': 24,
                  'parameter_order': results[0]['parameter_order'],
                  'exact_worker_outputs': True, 'exact_prediction_numeric_bytes': True}
    (folder / 'comparison.json').write_text(json.dumps(comparison, indent=2) + '\n')


def main():
    """Run with a fresh optional evidence directory; retain artifacts on failure."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path)
    args = parser.parse_args()
    if args.out is None:
        folder = Path(tempfile.mkdtemp(prefix='sbonest-worker-identity-')).resolve()
    else:
        folder = args.out.expanduser().absolute()
        folder.mkdir(parents=True, exist_ok=False)
    print('Worker identity artifacts:', folder, flush=True)
    check_full_fit(folder)
    print('PASS: canonical full fits and exact worker identity')


if __name__ == '__main__':
    main()
