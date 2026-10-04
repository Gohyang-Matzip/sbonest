"""Global (shared-exchange) versus individual model comparison on small synthetic data."""
# ruff: noqa: E402 -- Limit numerical libraries before importing them.
import os

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(name, '1')
os.environ.setdefault('MPLBACKEND', 'Agg')

import contextlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np

from sb_compare import _residue_config, information_criteria, nested_f_test, run_comparison
from sideband import composite_segments, profile

ROOT = Path(__file__).resolve().parent
SF_N = 600.


def write_synthetic(folder, rates, *, sigma=0.004, seed=3, sizes=(13, 15, 14), nus=(25., 80.)):
    """Three residues, two RF levels; ``rates`` maps label -> (kab, kba)."""
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    residues = [('A1', 120., 3., 6.2, 6.5), ('G2', 110., -2., 7.2, 7.5), ('S3', 115., 4., 9.8, 9.6)]
    rng = np.random.default_rng(seed)
    paths, initial = [], {}
    for nu in nus:
        path = folder / f'rf{int(nu)}.txt'
        lines = [f'{SF_N}', '0.03', f'{nu} 0', '# offset intensity error']
        for (label, peak, dw, ha, hb), size in zip(residues, sizes):
            kab, kba = rates[label]
            offsets = np.linspace(-240., 240., size)
            y = profile(composite_segments(70e-6), offsets, T=.03, nu=nu, kab=kab, kba=kba, dw=dw * SF_N,
                        r1=1.5, r2a=12., r2b=18., ha=(ha - 8.5) * 600., hb=(hb - 8.5) * 600., J=92., r1h=2., r2h=25.)
            y = y + rng.normal(0., sigma, size)
            lines.append(f'# {label} R2a: 12 R2b: 18 dw: {dw}')
            lines.extend(f'{peak + x / SF_N:.15g} {value:.15g} {sigma}' for x, value in zip(offsets, y))
            initial.update({f'{label}.{key}': value for key, value in zip(
                ('peak_ppm', 'dw_ppm', 'R1', 'R2a', 'R2b'), (peak, dw, 1.5, 12., 18.))})
        path.write_text('\n'.join(lines) + '\n')
        paths.append(path.name)
    config = {'Project Name': str(folder / 'fit'), 'datasets': paths,
              'residues': [{'name': r[0], 'flag': 'on'} for r in residues],
              'init': {'Method': 'Sideband', 'kex': {'min': 300., 'max': 300., 'nsteps': 1},
                       'pB': {'min': .05, 'max': .05, 'nsteps': 1}, 'initial': initial, 'max_nfev': 200,
                       'vary': ['kab', 'kba', 'v1n_scale'] + [f'{r[0]}.{k}' for r in residues
                                                              for k in ('peak_ppm', 'dw_ppm', 'R2b')]},
              'sideband': {'decoupling': {'h_larmor_mhz': 600., 'h_carrier_ppm': 8.5, 'p90_s': 70e-6,
                                          'R1H': 2., 'R2H': 25.},
                           'residues': {r[0]: {'h_ppm_a': r[3], 'h_ppm_b': r[4]} for r in residues},
                           'v1n': {'mode': 'scale', 'initial': 1., 'bounds': [.8, 1.2]}}}
    path = folder / 'fit.json'
    path.write_text(json.dumps(config, indent=1) + '\n')
    return path


def check_statistics_and_config_pruning():
    ic = information_criteria(100., 5, 50)
    assert ic['aic'] == 110. and abs(ic['aicc'] - (110. + 60. / 44.)) < 1e-12 and ic['bic'] > ic['aic']
    assert information_criteria(1., 10, 11)['aicc'] == float('inf')
    test = nested_f_test(120., 5, 100., 9, 60)
    assert test['extra_parameters'] == 4 and test['dof'] == 51 and 0 < test['p_value'] < 1
    assert nested_f_test(120., 9, 100., 5, 60)['p_value'] is None
    config = json.loads(write_synthetic(Path(tempfile.mkdtemp(prefix='sbonest-cmp-cfg-')),
                                        {'A1': (15., 285.), 'G2': (15., 285.), 'S3': (15., 285.)}).read_text())
    config['init']['bounds'] = {'A1.R2b': [1., 50.], 'G2.R2b': [1., 50.]}
    config['init']['multistart'] = {'starts': [{'A1.R2b': 20., 'G2.R2b': 21., 'kab': 20.}]}
    config['init']['profile'] = {'kex': [300.]}
    pruned = _residue_config(config, ['A1', 'G2', 'S3'], {'G2'})
    assert [r['flag'] for r in pruned['residues']] == ['off', 'on', 'off']
    assert set(pruned['init']['initial']) == {f'G2.{k}' for k in ('peak_ppm', 'dw_ppm', 'R1', 'R2a', 'R2b')}
    assert pruned['init']['bounds'] == {'G2.R2b': [1., 50.]}
    assert pruned['init']['vary'] == ['kab', 'kba', 'v1n_scale', 'G2.peak_ppm', 'G2.dw_ppm', 'G2.R2b']
    assert pruned['init']['multistart']['starts'] == [{'G2.R2b': 21., 'kab': 20.}]
    assert 'profile' not in pruned['init']
    assert config['init']['vary'][3] == 'A1.peak_ppm', 'Original configuration was modified'


def check_shared_and_distinct_exchange():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-compare-'))
    shared = write_synthetic(folder / 'shared', {'A1': (15., 285.), 'G2': (15., 285.), 'S3': (15., 285.)})
    distinct = write_synthetic(folder / 'distinct', {'A1': (15., 285.), 'G2': (15., 285.), 'S3': (60., 540.)})
    results = {}
    with contextlib.redirect_stdout(io.StringIO()):
        for label, path in (('shared', shared), ('distinct', distinct)):
            config = json.loads(path.read_text())
            config['datasets'] = [str(path.parent / name) for name in config['datasets']]
            paths = run_comparison(config, path.parent, folder / f'{label}_out', no_pdf=True)
            results[label] = json.loads(Path(paths['comparison_json']).read_text())
            assert Path(paths['comparison_pdf']).stat().st_size > 1000
            assert (folder / f'{label}_out' / 'global' / 'fit_checkpoint').is_dir()
            assert all((folder / f'{label}_out' / 'individual' / r / 'fit_result.json').exists()
                       for r in ('A1', 'G2', 'S3'))
    for label, summary in results.items():
        c = summary['comparison']
        assert c is not None, summary['warnings']
        # Each individual fit carries its own kab, kba and v1n_scale: 3 x (3 - 1) extra.
        assert c['individual_sum']['k'] == c['global']['k'] + 6
        assert c['f_test']['extra_parameters'] == 6 and c['n_points'] == 84
        assert c['delta_chi2'] >= -1e-6, 'Global chi2 should not beat the nesting individual model'
    assert results['shared']['comparison']['preferred_by_aicc'] == 'global'
    assert results['shared']['comparison']['f_test']['p_value'] > 0.01
    assert results['distinct']['comparison']['preferred_by_aicc'] == 'individual'
    assert results['distinct']['comparison']['f_test']['p_value'] < 1e-6
    s3 = results['distinct']['individual']['S3']
    # Small noisy data: only the direction of the per-residue estimates is asserted.
    assert s3['kex'] > 450. and results['distinct']['individual']['A1']['kex'] < 400.
    try:
        run_comparison(json.loads(shared.read_text()), shared.parent, folder / 'shared_out')
    except FileExistsError:
        pass
    else:
        raise AssertionError('Existing comparison output was overwritten')


def check_cli():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-compare-cli-'))
    path = write_synthetic(folder, {'A1': (15., 285.), 'G2': (15., 285.), 'S3': (15., 285.)})
    env = dict(os.environ, PYTHONPATH=str(ROOT))
    run = subprocess.run([sys.executable, str(ROOT / 'sb_workflow.py'), 'compare', str(path),
                          '--out', str(folder / 'cli_out'), '--workers', '2'],
                         capture_output=True, text=True, env=env, cwd=ROOT)
    assert run.returncode == 0, run.stderr[-2000:]
    assert (folder / 'cli_out' / 'comparison.txt').read_text().startswith('Sideband model comparison')
    bad = subprocess.run([sys.executable, str(ROOT / 'sb_workflow.py'), 'compare', str(path),
                          '--out', str(folder / 'cli_out2'), '--workers', '0'],
                         capture_output=True, text=True, env=env, cwd=ROOT)
    assert bad.returncode == 2


if __name__ == '__main__':
    check_statistics_and_config_pruning()
    check_shared_and_distinct_exchange()
    check_cli()
    print('PASS: global versus individual comparison statistics, fits and CLI')
