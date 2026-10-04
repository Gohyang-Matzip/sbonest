"""Multi-field relaxation groups and three-state Sideband models on synthetic data."""
# ruff: noqa: E402 -- Limit numerical libraries before importing them.
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]  # repository root: flat modules live there
sys.path.insert(0, str(ROOT))

import os

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(name, '1')
os.environ.setdefault('MPLBACKEND', 'Agg')

import contextlib
import io
import json
import subprocess
import tempfile

import numpy as np
from numpy.testing import assert_allclose

from sbfit import SidebandModel, run_config
from sideband import composite_segments, profile, profile_states, stationary_populations

SEGMENTS = composite_segments(70e-6)


def write_file(path, field, nu, rows, sigma, header):
    lines = [f'{field:.12g}', '0.05', f'{nu:g} 0', '# offset intensity error', header]
    lines.extend(f'{o:.15g} {v:.15g} {sigma:g}' for o, v in rows)
    path.write_text('\n'.join(lines) + '\n')


def two_field_config(folder, *, sigma=0.002, seed=5):
    """One residue at 600 and 800 MHz with field-specific true H and N relaxation."""
    folder.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    relax = {600.: dict(r1=1.5, r2a=12., r2b=15., r1h=2., r2h=25.),
             800.: dict(r1=1.2, r2a=15., r2b=19., r1h=1.5, r2h=32.)}
    datasets, overrides = [], []
    for larmor in (600., 800.):
        sf_n = larmor * 0.101329118
        for nu in (25., 100.):
            offsets = np.linspace(-800., 800., 41)
            r = relax[larmor]
            y = profile(SEGMENTS, offsets, T=.05, nu=nu * 1.05, kab=15., kba=285., dw=3. * sf_n,
                        r1=r['r1'], r2a=r['r2a'], r2b=r['r2b'], ha=(6.2 - 8.5) * larmor, hb=(6.5 - 8.5) * larmor,
                        J=92., r1h=r['r1h'], r2h=r['r2h'])
            y = y + rng.normal(0., sigma, len(y))
            path = folder / f'f{int(larmor)}_rf{int(nu)}.txt'
            write_file(path, sf_n, nu, zip(120. + offsets / sf_n, y), sigma, '# A1 R2a: 13 R2b: 18 dw: 2.8')
            datasets.append(path.name)
            overrides.append({'h_larmor_mhz': larmor})
    config = {'Project Name': str(folder / 'fit'), 'datasets': datasets,
              'residues': [{'name': 'A1', 'flag': 'on'}],
              'init': {'Method': 'Sideband', 'kex': {'min': 200., 'max': 400., 'nsteps': 3},
                       'pB': {'min': .03, 'max': .07, 'nsteps': 3},
                       'initial': {'A1.peak_ppm': 120.02, 'A1.R1': 1.3}, 'max_nfev': 300},
              'sideband': {'decoupling': {'h_carrier_ppm': 8.5, 'p90_s': 70e-6},
                           'datasets': overrides,
                           'residues': {'A1': {'h_ppm_a': 6.2, 'h_ppm_b': 6.5}},
                           'proton_relaxation': {'mode': 'fit'},
                           'nitrogen_relaxation': {'mode': 'per_field'},
                           'v1n': {'mode': 'scale', 'initial': 1., 'bounds': [.8, 1.2]}}}
    (folder / 'fit.json').write_text(json.dumps(config, indent=1) + '\n')
    return config, relax


def check_two_fields_per_field_relaxation():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-two-fields-'))
    config, relax = two_field_config(folder)
    model = SidebandModel(config, folder)
    assert model.field_groups == [600., 800.] and model.dataset_group == [0, 0, 1, 1]
    expected = ['kab', 'kba', 'v1n_scale', 'A1.peak_ppm', 'A1.dw_ppm',
                'A1.R1[0]', 'A1.R2a[0]', 'A1.R2b[0]', 'A1.R1[1]', 'A1.R2a[1]', 'A1.R2b[1]',
                'A1.R1H[0]', 'A1.R2H[0]', 'A1.R1H[1]', 'A1.R2H[1]']
    assert model.parameter_names == expected, model.parameter_names
    # The expanded initial vector copies one relaxation set per group and default H rates.
    p0 = model.prepare_fit()
    assert p0[5:8].tolist() == p0[8:11].tolist() and p0[11:13].tolist() == [2., 25.] == p0[13:15].tolist()
    assert model.physical_bounds('A1.R2a[1]') == (0., np.inf) and model.physical_bounds('A1.dw_ppm') == (-np.inf, np.inf)
    config['init']['vary'] = [name for name in expected if name not in ('A1.R1H[0]', 'A1.R1H[1]')]
    config['init']['initial'].update({'A1.R1H[0]': 2., 'A1.R1H[1]': 1.5})
    # Ungrouped aliases address every group; explicit grouped entries take precedence.
    config['init']['initial'].update({'A1.R2b': 17., 'A1.R2b[1]': 18.})
    config['init']['multistart'] = {'starts': [{'A1.R2a': 14.}]}
    normalized = model.init_config()
    assert normalized['initial']['A1.R2b[0]'] == 17. and normalized['initial']['A1.R2b[1]'] == 18.
    assert normalized['multistart']['starts'] == [{'A1.R2a[0]': 14., 'A1.R2a[1]': 14.}]
    assert model.expand_name('A1.R2a') == ['A1.R2a[0]', 'A1.R2a[1]'] and model.expand_name('A1.R2a[0]') == ['A1.R2a[0]']
    assert model.expand_name('nonsense') == []
    with contextlib.redirect_stdout(io.StringIO()):
        run_config(config, folder, no_pdf=True, workers=1)
    result = json.loads((folder / 'fit_result.json').read_text())
    assert result['success'] and result['states'] == 2 and result['model'] == 'Sideband'
    assert len(result['multistart']) == 2 and all(row['success'] for row in result['multistart'])
    assert [g['h_larmor_mhz'] for g in result['field_groups']] == [600., 800.]
    values = {name: item['value'] for name, item in result['parameters'].items()}
    assert abs(values['kab'] + values['kba'] - 300.) < 15., values
    assert abs(values['v1n_scale'] - 1.05) < 0.01
    for g, larmor in enumerate((600., 800.)):
        assert abs(values[f'A1.R2a[{g}]'] - relax[larmor]['r2a']) < 1.5, (g, values[f'A1.R2a[{g}]'])
        assert abs(values[f'A1.R2b[{g}]'] - relax[larmor]['r2b']) < 3., (g, values[f'A1.R2b[{g}]'])
        assert abs(values[f'A1.R1[{g}]'] - relax[larmor]['r1']) < 0.3
    text = (folder / 'fit_result.txt').read_text()
    assert 'Nitrogen relaxation mode: per_field' in text and 'field group 1: 800 MHz' in text
    assert 'A1.R2H[1] [s-1]' in text
    # Shared nitrogen relaxation across fields keeps the single-field names.
    shared = json.loads((folder / 'fit.json').read_text())
    shared['sideband']['nitrogen_relaxation'] = {'mode': 'shared'}
    names = SidebandModel(shared, folder).parameter_names
    assert 'A1.R1' in names and 'A1.R1H[1]' in names and 'A1.R1[0]' not in names


def three_state_config(folder, *, sigma=0.002, seed=9, fixed_h=True):
    """One residue, linear A<->B<->C exchange with distinct minor-state shifts."""
    folder.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    sf_n = 60.7974708
    larmor = 600.
    K = np.array([[0., 15., 0.], [285., 0., 100.], [0., 50., 0.]])
    shifts = [0., 3. * sf_n, -4. * sf_n]
    h_shifts = [(6.2 - 8.5) * larmor, (6.5 - 8.5) * larmor, (7.1 - 8.5) * larmor]
    datasets = []
    for nu in (25., 60., 120.):
        offsets = np.linspace(-800., 800., 65)
        y = profile_states(SEGMENTS, offsets, T=.05, nu=nu, exchange=K, shifts=shifts, h_shifts=h_shifts,
                           r1=1.5, r2=[12., 15., 20.], J=92., r1h=2., r2h=25.)
        y = y + rng.normal(0., sigma, len(y))
        path = folder / f'rf{int(nu)}.txt'
        write_file(path, sf_n, nu, zip(120. + offsets / sf_n, y), sigma, '# A1 R2a: 12 R2b: 15 dw: 3')
        datasets.append(path.name)
    decoupling = {'h_larmor_mhz': larmor, 'h_carrier_ppm': 8.5, 'p90_s': 70e-6}
    if fixed_h:
        decoupling.update(R1H=2., R2H=25.)
    config = {'Project Name': str(folder / 'fit'), 'datasets': datasets,
              'residues': [{'name': 'A1', 'flag': 'on'}],
              'init': {'Method': 'Sideband_3st_Linear',
                       'initial': {'kab': 12., 'kba': 300., 'kbc': 80., 'kcb': 60., 'A1.peak_ppm': 120.,
                                   'A1.dwC_ppm': -3.5, 'A1.R2c': 18., 'A1.R1': 1.4},
                       'max_nfev': 400},
              'sideband': {'decoupling': decoupling,
                           'residues': {'A1': {'h_ppm_a': 6.2, 'h_ppm_b': 6.5, 'h_ppm_c': 7.1}},
                           'v1n': {'mode': 'fixed'}}}
    (folder / 'fit.json').write_text(json.dumps(config, indent=1) + '\n')
    return config, K


def check_three_state_model():
    assert_allclose(stationary_populations(np.array([[0., 15.], [285., 0.]])), [.95, .05], rtol=0, atol=1e-15)
    folder = Path(tempfile.mkdtemp(prefix='sbonest-three-state-'))
    config, K = three_state_config(folder)
    model = SidebandModel(config, folder)
    assert model.states == 3 and model.rate_names == ['kab', 'kba', 'kbc', 'kcb']
    assert model.parameter_names == ['kab', 'kba', 'kbc', 'kcb', 'A1.peak_ppm', 'A1.dw_ppm', 'A1.dwC_ppm',
                                     'A1.R1', 'A1.R2a', 'A1.R2b', 'A1.R2c']
    assert model.physical_bounds('kcb') == (1e-8, np.inf) and model.physical_bounds('kbc') == (0., np.inf)
    assert model.physical_bounds('A1.dwC_ppm') == (-np.inf, np.inf)
    with contextlib.redirect_stdout(io.StringIO()):
        run_config(config, folder, no_pdf=True)
    result = json.loads((folder / 'fit_result.json').read_text())
    assert result['success'] and result['states'] == 3 and result['model'] == 'Sideband_3st_Linear'
    assert result['kex'] is None and result['pB'] is None and result['derived_se'] == {'kex': None, 'pB': None}
    rates = result['exchange']['rates']
    truth = dict(zip(['kab', 'kba', 'kbc', 'kcb'], [15., 285., 100., 50.]))
    for name, value in truth.items():
        assert abs(rates[name] - value) < 0.25 * value + 3., (name, rates[name])
    populations = result['exchange']['populations']
    expected = stationary_populations(K)
    assert abs(populations['A'] - expected[0]) < 0.02 and abs(populations['C'] - expected[2]) < 0.03
    values = {name: item['value'] for name, item in result['parameters'].items()}
    assert abs(values['A1.dwC_ppm'] + 4.) < 0.2 and abs(values['A1.dw_ppm'] - 3.) < 0.2
    assert abs(values['A1.R2c'] - 20.) < 5.
    text = (folder / 'fit_result.txt').read_text()
    assert 'kbc [s-1]' in text and 'Sideband' in text
    # Reports, bootstrap and the two-state-only analyses behave as documented.
    from sb_report import regenerate_report
    paths = regenerate_report(folder / 'fit_result.json', folder / 'report')
    assert 'kex: unavailable' in Path(paths['summary_txt']).read_text()
    boot = json.loads((folder / 'fit.json').read_text())
    boot['Project Name'] = str(folder / 'boot')
    boot['init']['bootstrap'] = {'replicates': 1, 'seed': 3}
    with contextlib.redirect_stdout(io.StringIO()):
        run_config(boot, folder, no_pdf=True)
    sample = json.loads((folder / 'boot_result.json').read_text())['bootstrap']
    assert sample['samples'][0]['kex'] is None and 'kex' not in sample['intervals'] and 'kcb' in sample['intervals']
    for key, settings in (('profile', {'kex': [300.]}), ('profile_interval', {'parameters': ['pB']})):
        bad = json.loads((folder / 'fit.json').read_text())
        bad['Project Name'] = str(folder / f'bad_{key}')
        bad['init'][key] = settings
        try:
            run_config(bad, folder, no_pdf=True)
        except ValueError as exc:
            assert 'two-state' in str(exc)
        else:
            raise AssertionError(f'{key} of kex/pB accepted for a three-state model')
    # Missing h_ppm_c and a wrong Method are rejected before fitting.
    for mutate in (lambda c: c['sideband']['residues']['A1'].pop('h_ppm_c'),
                   lambda c: c['init'].update(Method='Sideband_3st_Square')):
        bad = json.loads((folder / 'fit.json').read_text())
        mutate(bad)
        try:
            SidebandModel(bad, folder)
        except ValueError:
            pass
        else:
            raise AssertionError('Invalid three-state configuration accepted')


def check_triangle_and_check_cli():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-triangle-'))
    config, _ = three_state_config(folder)
    config['init']['Method'] = 'Sideband_3st_Triangle'
    config['init']['initial'].update(kca=5., kac=1.)
    (folder / 'fit.json').write_text(json.dumps(config, indent=1) + '\n')
    model = SidebandModel(config, folder)
    assert model.rate_names == ['kab', 'kba', 'kbc', 'kcb', 'kca', 'kac']
    P = model.seParam(model.prepare_fit())
    K = model.exchange_matrix(P)
    assert K[2, 0] == 5. and K[0, 2] == 1. and K[1, 2] == 80.
    env = dict(os.environ, PYTHONPATH=str(ROOT))
    run = subprocess.run([sys.executable, str(ROOT / 'run.py'), str(folder / 'fit.json'), '--check',
                          '--identifiability', '--no-pdf'], capture_output=True, text=True, env=env, cwd=ROOT)
    assert run.returncode in (0, 1), run.stderr[-1500:]
    summary = json.loads(run.stdout)
    assert summary['valid'] and summary['model'] == 'Sideband_3st_Triangle' and summary['states'] == 3
    assert summary['identifiability']['derived_se'] == {'kex': None, 'pB': None}
    assert summary['n_parameters'] == 13


if __name__ == '__main__':
    check_two_fields_per_field_relaxation()
    check_three_state_model()
    check_triangle_and_check_cli()
    print('PASS: field groups, per-field relaxation, three-state models and CLI')
