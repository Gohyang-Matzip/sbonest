"""Lossless Sideband reports and calculation-time provenance."""
import csv
from datetime import datetime, timezone
import hashlib
from importlib.metadata import version
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

import numpy as np


def derived_errors(p, covariance):
    """Propagate correlated exchange rates using supplied absolute-sigma covariance."""
    total = p[0] + p[1]
    gradients = (np.ones(2), np.array([p[1], -p[0]]) / total**2)
    rates = covariance[:2, :2]
    return {
        name: float(np.sqrt(max(0., gradient @ rates @ gradient)))
        if np.isfinite(rates).all() else None
        for name, gradient in zip(('kex', 'pB'), gradients)
    }


def provenance(config, config_dir):
    """Snapshot files before fitting, including the actual waveform bytes."""
    root = Path(__file__).resolve().parent

    def digest(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    def file_record(path):
        path = (Path(config_dir) / path).resolve()
        return {'path': str(path), 'sha256': digest(path)}

    sources = ('run.py', 'sbfit.py', 'sideband.py', 'fit.py', 'estmodel.py',
               'est_data.py', 'sb_report.py', 'sb_analysis.py')
    sb = config['sideband']
    waveforms = set()
    for override in sb.get('datasets', [{}] * len(config['datasets'])):
        entry = {**sb['decoupling'], **override}
        if 'waveform_json' in entry:
            waveforms.add(entry['waveform_json'])
    git = {}
    try:
        git['commit'] = subprocess.check_output(
            ['git', '-C', str(root), 'rev-parse', 'HEAD'], text=True,
            stderr=subprocess.DEVNULL).strip()
        git['dirty'] = bool(subprocess.check_output(
            ['git', '-C', str(root), 'status', '--porcelain', '--untracked-files=normal'],
            text=True, stderr=subprocess.DEVNULL).strip())
    except (OSError, subprocess.CalledProcessError):
        git = {'commit': None, 'dirty': None}
    return {
        'started_utc': datetime.now(timezone.utc).isoformat(),
        'python': sys.version,
        'platform': platform.platform(),
        'packages': {name: version(name) for name in
                     ('numpy', 'scipy', 'matplotlib', 'optimalcontrol-nmr')},
        'threads': {name: os.environ.get(name) for name in
                    ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')},
        'git': git,
        'config_sha256': hashlib.sha256(json.dumps(
            config, sort_keys=True, allow_nan=False).encode()).hexdigest(),
        'datasets': [file_record(path) for path in config['datasets']],
        'waveforms': [file_record(path) for path in sorted(waveforms)],
        'source_sha256': {name: digest(root / name) for name in sources},
    }


def prediction_rows(model, p):
    """Full precision observations and signed standardized residuals in input order."""
    parameters = model.seParam(p)
    rows = []
    for index, es, offsets, observed, sigma in model._prepare_data():
        predicted = model.calc(parameters, index, offsets, es)
        rf, _ = model.rf_values(parameters, es)
        for x, y, fitted, error in zip(offsets, observed, predicted, sigma):
            rows.append({
                'residue': model.dataset.res[index].label,
                'dataset_index': es.dataset_index,
                'field_mhz': es.field,
                'saturation_s': es.T,
                'v1n_nominal_hz': es.v1,
                'v1n_hz': rf,
                'offset_ppm': float(x),
                'observed': float(y),
                'predicted': float(fitted),
                'sigma': float(error),
                'residual_sigma': float((y - fitted) / error),
            })
    return rows


def write_predictions(path, rows):
    with Path(path).open('x', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def fit_pdf(path, rows):
    """One page per residue with profile and standardized-residual panels."""
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    with PdfPages(path) as pdf:
        for residue in dict.fromkeys(row['residue'] for row in rows):
            selected = [row for row in rows if row['residue'] == residue]
            fig, (ax, residual_ax) = plt.subplots(
                2, 1, figsize=(8.5, 7), sharex=True,
                gridspec_kw={'height_ratios': [3, 1]}, layout='constrained')
            for dataset in dict.fromkeys(row['dataset_index'] for row in selected):
                points = sorted((row for row in selected if row['dataset_index'] == dataset),
                                key=lambda row: row['offset_ppm'])
                x = [row['offset_ppm'] for row in points]
                first = points[0]
                label = (f"dataset {dataset}: {first['v1n_hz']:.3g} Hz, "
                         f"{first['saturation_s'] * 1000:g} ms")
                artist = ax.errorbar(x, [row['observed'] for row in points],
                                    yerr=[row['sigma'] for row in points],
                                    fmt='o', markersize=3, label=label)
                color = artist[0].get_color()
                ax.plot(x, [row['predicted'] for row in points], color=color)
                residual_ax.plot(x, [row['residual_sigma'] for row in points],
                                 'o-', color=color, markersize=3, linewidth=.6)
            ax.set(title=residue, ylabel='Intensity (I/I0)')
            ax.legend(fontsize='small')
            ax.grid(alpha=.25)
            residual_ax.axhline(0, color='black', linewidth=.7)
            residual_ax.set(xlabel='Chemical shift offset (ppm)', ylabel='Residual / σ')
            residual_ax.grid(alpha=.25)
            pdf.savefig(fig)
            plt.close(fig)
