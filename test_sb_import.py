"""Bruker pseudo-2D import on a synthetic processed directory with submatrix layout."""
# ruff: noqa: E402 -- Limit numerical libraries before importing them.
import os

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(name, '1')

import json
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np
from numpy.testing import assert_allclose

from est_data import EstDataSet
from sb_import import convert, find_peaks, peaks_from_reference, read_offsets, read_parameters, read_pseudo2d

ROOT = Path(__file__).resolve().parent
SF = 600.13
SW = 7200.0
SI2, SI1, XDIM2, XDIM1 = 256, 24, 64, 8
OFFSET_PPM = 12.0


def write_bruker(folder, *, byte_order=0, dtype_code=0, seed=1):
    """Synthetic 1H rows: two Lorentzian peaks whose heights follow known profiles."""
    folder.mkdir(parents=True, exist_ok=True)
    ppm = OFFSET_PPM - np.arange(SI2) * (SW / SF) / SI2
    rng = np.random.default_rng(seed)
    offsets_n = np.linspace(110., 130., SI1 - 1)
    profile_a = 1. - 0.6 * np.exp(-0.5 * ((offsets_n - 118.) / 1.2) ** 2)
    profile_g = 1. - 0.4 * np.exp(-0.5 * ((offsets_n - 123.) / 1.5) ** 2)
    rows = np.zeros((SI1, SI2))
    scale = 50000.
    for row in range(SI1):
        fa = 1. if row == 0 else profile_a[row - 1]   # row 0 is the reference
        fg = 1. if row == 0 else profile_g[row - 1]
        rows[row] = (scale * fa / (1 + ((ppm - 8.30) / 0.015) ** 2)
                     + 0.8 * scale * fg / (1 + ((ppm - 7.60) / 0.015) ** 2)
                     + rng.normal(0., 0.004 * scale, SI2))
    nc = -3 if dtype_code == 0 else 0
    stored = rows / 2.0 ** nc
    blocks = stored.reshape(SI1 // XDIM1, XDIM1, SI2 // XDIM2, XDIM2).transpose(0, 2, 1, 3)
    dtype = np.dtype(np.int32 if dtype_code == 0 else np.float64).newbyteorder('>' if byte_order else '<')
    (folder / '2rr').write_bytes(np.ascontiguousarray(blocks).astype(dtype).tobytes())
    (folder / 'procs').write_text('\n'.join([
        '##TITLE= Parameter file', f'##$BYTORDP= {byte_order}', f'##$DTYPP= {dtype_code}',
        f'##$NC_proc= {nc}', f'##$OFFSET= {OFFSET_PPM}', f'##$SF= {SF}', f'##$SI= {SI2}',
        f'##$SW_p= {SW}', f'##$XDIM= {XDIM2}', '##$TITLE= <cest>', '##$LIST= (0..2)', '1 2 3', '##END=']) + '\n')
    (folder / 'proc2s').write_text('\n'.join([
        '##TITLE= Parameter file', f'##$SI= {SI1}', f'##$XDIM= {XDIM1}', '##END=']) + '\n')
    (folder / 'offsets.txt').write_text('# reference row first\n1000\n' + '\n'.join(f'{o:.10g}' for o in offsets_n) + '\n')
    return folder, ppm, offsets_n, profile_a, profile_g, rows


def check_reader_and_parameters():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-bruker-'))
    for byte_order, dtype_code in ((0, 0), (1, 0), (0, 2)):
        pdata, ppm, _, _, _, rows = write_bruker(folder / f'b{byte_order}_d{dtype_code}', byte_order=byte_order, dtype_code=dtype_code)
        data, axis = read_pseudo2d(pdata)
        assert data.shape == (SI1, SI2) and axis['rows'] == SI1 and axis['points'] == SI2
        assert_allclose(axis['ppm'], ppm)
        tolerance = 2.0 ** -3 if dtype_code == 0 else 1e-9
        assert np.max(np.abs(data - rows)) <= tolerance * 1.01, 'Submatrix layout or scaling mismatch'
    params = read_parameters(pdata / 'procs')
    assert params['TITLE'] == 'cest' and params['LIST'] == [1, 2, 3] and params['SF'] == SF
    assert read_offsets(pdata / 'offsets.txt', unit='ppm').shape == (SI1,)
    hz = Path(folder / 'hz.txt')
    hz.write_text('0\n-600\n600\n')
    assert_allclose(read_offsets(hz, unit='hz', carrier_ppm=120., field_mhz=60.), [120., 110., 130.])
    try:
        read_offsets(hz, unit='hz')
    except ValueError:
        pass
    else:
        raise AssertionError('Hz offsets without a carrier were accepted')


def check_convert_and_reload():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-bruker-convert-'))
    pdata, ppm, offsets_n, profile_a, profile_g, _ = write_bruker(folder / 'pdata')
    out = folder / 'cest_25.txt'
    summary = convert(pdata, pdata / 'offsets.txt', out, peaks={'A1': 8.30, 'G2': 7.60}, half_width=0.03,
                      reference_row=0, noise_region=(2.0, 4.0), saturation_s=0.4, v1_hz=25.,
                      dw_ppm={'A1': 2.5, 'G2': -1.}, r2a=11., r2b=19.)
    assert abs(summary['field_mhz'] - SF * 0.101329118) < 1e-9
    assert summary['peaks']['A1']['n_points'] == SI1 - 1 and summary['excluded_rows'] == [0]
    dataset = EstDataSet()
    dataset.addData(str(out))
    assert abs(dataset.fields[0] - summary['field_mhz']) < 1e-8 and dataset.Ts == [0.4] and dataset.v1s == [25.]
    labels = [r.label for r in dataset.res]
    assert labels == ['A1', 'G2']
    for residue, truth in zip(dataset.res, (profile_a, profile_g)):
        es = residue.estSpecs[0]
        assert_allclose(es.offset, offsets_n, rtol=0, atol=1e-5)
        assert np.max(np.abs(np.asarray(es.int) - truth)) < 0.03, residue.label
        assert all(abs(e - es.intstd[0]) < 1e-12 for e in es.intstd) and 0.002 < es.intstd[0] < 0.02
        assert es.initr2a == 11. and es.initr2b == 19.
    assert dataset.res[0].estSpecs[0].initdw == 2.5 and dataset.res[1].estSpecs[0].initdw == -1.
    try:
        convert(pdata, pdata / 'offsets.txt', out, peaks={'A1': 8.30}, half_width=0.03, reference_row=0,
                noise_region=(2.0, 4.0), saturation_s=0.4, v1_hz=25.)
    except FileExistsError:
        pass
    else:
        raise AssertionError('Existing output was overwritten')
    for bad in (dict(reference_row=99), dict(saturation_s=0.), dict(noise_region=(8.299, 8.301)),
                dict(peaks={'A1': 20.})):
        kwargs = dict(peaks={'A1': 8.30}, half_width=0.03, reference_row=0, noise_region=(2.0, 4.0),
                      saturation_s=0.4, v1_hz=25.)
        kwargs.update(bad)
        try:
            convert(pdata, pdata / 'offsets.txt', folder / 'bad.txt', **kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError(f'Invalid conversion accepted: {bad}')
        assert not (folder / 'bad.txt').exists()


def check_cli():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-bruker-cli-'))
    pdata, *_ = write_bruker(folder / 'pdata')
    env = dict(os.environ, PYTHONPATH=str(ROOT))
    run = subprocess.run([sys.executable, str(ROOT / 'sb_import.py'), str(pdata), str(pdata / 'offsets.txt'),
                          '--out', str(folder / 'cli.txt'), '--peak', 'A1=8.30:2.5', '--peak', 'G2=7.6',
                          '--half-width', '0.03', '--reference-row', '0', '--noise-region', '2', '4',
                          '--saturation-s', '0.4', '--v1-hz', '25', '--mode', 'sum'],
                         capture_output=True, text=True, env=env, cwd=ROOT)
    assert run.returncode == 0, run.stderr[-1500:]
    summary = json.loads(run.stdout)
    assert set(summary['peaks']) == {'A1', 'G2'} and (folder / 'cli.txt').exists()
    bad = subprocess.run([sys.executable, str(ROOT / 'sb_import.py'), str(pdata), str(pdata / 'offsets.txt'),
                          '--out', str(folder / 'cli2.txt'), '--peak', 'A=8.3', '--reference-row', '0',
                          '--noise-region', '2', '4', '--saturation-s', '0.4', '--v1-hz', '25'],
                         capture_output=True, text=True, env=env, cwd=ROOT)
    assert bad.returncode == 2 and 'residue number' in bad.stderr


def check_fqlist_peaks_and_qa():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-bruker-qa-'))
    pdata, ppm, offsets_n, profile_a, profile_g, rows = write_bruker(folder / 'pdata')
    # Bruker fq list headers select the unit; O1 lines are ignored; Hz lists need the carrier.
    fq = folder / 'fq1list'
    fq.write_text('bf ppm\n1000\n' + '\n'.join(f'{o:.10g}' for o in offsets_n) + '\n')
    assert_allclose(read_offsets(fq, unit='hz')[1:], offsets_n)
    fq_hz = folder / 'fq2list'
    fq_hz.write_text('sfo hz\nO1 7300.0\n0\n-600\n600\n')
    assert_allclose(read_offsets(fq_hz, unit='ppm', carrier_ppm=120., field_mhz=60.), [120., 110., 130.])
    (folder / 'p.txt').write_text('P\n1.5 2.5\n')
    assert_allclose(read_offsets(folder / 'p.txt', unit='hz'), [1.5, 2.5])
    (folder / 'bad.txt').write_text('1.0\nabc\n')
    try:
        read_offsets(folder / 'bad.txt', unit='ppm')
    except ValueError:
        pass
    else:
        raise AssertionError('Non-numeric offset accepted')
    # Peak extraction from the reference row finds the two synthetic peaks.
    found = peaks_from_reference(pdata, 0, (2.0, 4.0), snr=10.)
    positions = sorted(position for _, position, _ in found)
    assert len(found) == 2 and abs(positions[0] - 7.60) < 0.02 and abs(positions[1] - 8.30) < 0.02, found
    assert [label for label, _, _ in found] == ['P1', 'P2'] and found[0][2] >= found[1][2]
    assert find_peaks(np.zeros(50), np.linspace(0, 1, 50), threshold=1.) == []
    # The QA figure is written next to the conversion and refuses to overwrite.
    out = folder / 'cest.txt'
    summary = convert(pdata, fq, out, peaks={'A1': 8.30, 'G2': 7.60}, half_width=0.03, reference_row=0,
                      noise_region=(2.0, 4.0), saturation_s=0.4, v1_hz=25., qa_path=folder / 'qa.pdf')
    assert Path(summary['qa_pdf']).read_bytes().startswith(b'%PDF-')
    try:
        convert(pdata, fq, folder / 'again.txt', peaks={'A1': 8.30}, half_width=0.03, reference_row=0,
                noise_region=(2.0, 4.0), saturation_s=0.4, v1_hz=25., qa_path=folder / 'qa.pdf')
    except FileExistsError:
        pass
    else:
        raise AssertionError('QA PDF was overwritten')
    env = dict(os.environ, PYTHONPATH=str(ROOT))
    run = subprocess.run([sys.executable, str(ROOT / 'sb_import.py'), str(pdata), str(fq), '--out', str(folder / 'auto.txt'),
                          '--peaks-from-reference', '--reference-row', '0', '--noise-region', '2', '4',
                          '--saturation-s', '0.4', '--v1-hz', '25', '--qa-pdf', str(folder / 'auto_qa.pdf')],
                         capture_output=True, text=True, env=env, cwd=ROOT)
    assert run.returncode == 0, run.stderr[-1500:]
    assert set(json.loads(run.stdout)['peaks']) == {'P1', 'P2'} and (folder / 'auto_qa.pdf').exists()
    neither = subprocess.run([sys.executable, str(ROOT / 'sb_import.py'), str(pdata), str(fq), '--out', str(folder / 'x.txt'),
                              '--reference-row', '0', '--noise-region', '2', '4', '--saturation-s', '0.4', '--v1-hz', '25'],
                             capture_output=True, text=True, env=env, cwd=ROOT)
    assert neither.returncode == 2


if __name__ == '__main__':
    check_reader_and_parameters()
    check_convert_and_reload()
    check_cli()
    check_fqlist_peaks_and_qa()
    print('PASS: Bruker pseudo-2D reader, conversion, reload, CLI, fq lists, peak extraction and QA figure')
