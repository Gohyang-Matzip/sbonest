"""Lossless Sideband reports and calculation-time provenance."""
import csv
from datetime import datetime, timezone
import hashlib
from io import BytesIO
from importlib.metadata import version
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import textwrap

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
               'est_data.py', 'sb_report.py', 'sb_analysis.py', 'sb_workflow.py',
               'sb_checkpoint.py', 'sb_bootstrap.py')
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


def fit_pdf(path, rows, summary=None):
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
        if summary is not None:
            _summary_pdf_pages(pdf, summary)


PREDICTION_COLUMNS = ('residue', 'dataset_index', 'field_mhz', 'saturation_s',
                      'v1n_nominal_hz', 'v1n_hz', 'offset_ppm', 'observed',
                      'predicted', 'sigma', 'residual_sigma')


def _finite(value, label, *, nonnegative=False):
    if (isinstance(value, bool) or not isinstance(value, (int, float))
            or not np.isfinite(value) or (nonnegative and value < 0)):
        raise ValueError(f'{label} must be a finite' + (' nonnegative' if nonnegative else '') + ' number')
    return float(value)


def _load_report_inputs(result_path, predictions_path):
    """Reconcile saved evidence before creating directories or report files."""
    def reject_constant(value):
        raise ValueError(f'Non-finite JSON value: {value}')

    info = json.loads(result_path.read_text(encoding='utf-8'), parse_constant=reject_constant)
    if not isinstance(info, dict) or info.get('method') != 'Sideband' or info.get('success', True) is not True:
        raise ValueError('Reporting requires a successful saved Sideband result')
    expected_chi2 = _finite(info.get('chi2'), 'Result chi2', nonnegative=True)
    count = info.get('n_points')
    if isinstance(count, bool) or not isinstance(count, int) or count <= 0:
        raise ValueError('Result n_points must be a positive integer')
    parameters = info.get('parameters')
    if not isinstance(parameters, dict) or not parameters:
        raise ValueError('Result parameters are missing')
    for name, item in parameters.items():
        if not isinstance(item, dict):
            raise ValueError(f'Invalid parameter {name}')
        _finite(item.get('value'), f'{name}.value')
        if item.get('stderr') is not None:
            _finite(item['stderr'], f'{name}.stderr', nonnegative=True)
        if not isinstance(item.get('vary'), bool):
            raise ValueError(f'{name}.vary must be boolean')
    config = info.get('config', {})
    if not isinstance(config, dict) or not isinstance(config.get('datasets'), list) or not config['datasets']:
        raise ValueError('Saved result must identify its datasets')
    # Fitted parameter names also identify active residues omitted from config.residues.
    residues = {name.rsplit('.', 1)[0] for name in parameters if '.' in name}
    if not residues:
        entries = config.get('residues', [])
        if not isinstance(entries, list) or not all(isinstance(row, dict) for row in entries):
            raise ValueError('Invalid saved residue metadata')
        residues = {row.get('name') for row in entries if row.get('flag') == 'on'}
    if not residues or not all(isinstance(name, str) and name for name in residues):
        raise ValueError('Saved result must identify its fitted residues')
    rf_values = info.get('v1n_hz')
    if not isinstance(rf_values, list) or len(rf_values) != len(config['datasets']):
        raise ValueError('Saved RF values must align with datasets')
    for value in rf_values:
        _finite(value, 'v1n_hz', nonnegative=True)
    rows, conditions = [], {}
    with predictions_path.open(newline='', encoding='utf-8') as stream:
        reader = csv.DictReader(stream)
        columns = reader.fieldnames or []
        if len(set(columns)) != len(columns) or not set(PREDICTION_COLUMNS).issubset(columns):
            raise ValueError('Predictions CSV is missing required columns or repeats a column')
        for number, raw in enumerate(reader, 2):
            if None in raw or any(value is None for value in raw.values()):
                raise ValueError(f'Predictions CSV row {number} has the wrong number of columns')
            try:
                row = {key: float(raw[key]) for key in PREDICTION_COLUMNS[2:]}
                row.update(residue=raw['residue'], dataset_index=int(raw['dataset_index']))
            except (ValueError, TypeError) as exc:
                raise ValueError(f'Invalid numeric value in predictions row {number}') from exc
            if not all(np.isfinite(row[key]) for key in PREDICTION_COLUMNS[2:]):
                raise ValueError(f'Non-finite predictions row {number}')
            ds = row['dataset_index']
            if row['residue'] not in residues or not 0 <= ds < len(config['datasets']):
                raise ValueError(f'Unknown residue or dataset in predictions row {number}')
            if (row['sigma'] <= 0 or row['field_mhz'] <= 0 or row['saturation_s'] < 0
                    or row['v1n_nominal_hz'] < 0 or row['v1n_hz'] < 0):
                raise ValueError(f'Invalid uncertainty or acquisition values in predictions row {number}')
            expected = (row['observed'] - row['predicted']) / row['sigma']
            if not np.isclose(expected, row['residual_sigma'], rtol=1e-10, atol=1e-10):
                raise ValueError(f'Inconsistent standardized residual in predictions row {number}')
            if not np.isclose(row['v1n_hz'], rf_values[ds], rtol=1e-10, atol=1e-10):
                raise ValueError(f'CSV RF differs from saved result in row {number}')
            acquisition = tuple(row[key] for key in PREDICTION_COLUMNS[2:6])
            if ds in conditions and not np.allclose(acquisition, conditions[ds], rtol=1e-12, atol=1e-12):
                raise ValueError(f'Inconsistent dataset acquisition values in row {number}')
            conditions[ds] = acquisition
            rows.append(row)
    if len(rows) != count:
        raise ValueError(f'Predictions row count {len(rows)} differs from result n_points {count}')
    if {row['residue'] for row in rows} != residues:
        raise ValueError('Predictions are missing fitted residues')
    chi2 = float(np.sum(np.square([row['residual_sigma'] for row in rows])))
    if not np.isfinite(chi2) or not np.isclose(chi2, expected_chi2, rtol=1e-9, atol=1e-9):
        raise ValueError(f'Predictions chi2 {chi2} differs from saved result {expected_chi2}')
    _validate_analyses(info)
    return info, rows


def _validate_analyses(info):
    attempts = info.get('multistart', [])
    if not isinstance(attempts, list) or not all(isinstance(row, dict) for row in attempts):
        raise ValueError('Invalid saved multistart records')
    for index, row in enumerate(attempts):
        if isinstance(row.get('index'), bool) or row.get('index') != index:
            raise ValueError('Multistart indices must be an ordered sequence from zero')
        if not isinstance(row.get('success'), bool):
            raise ValueError('Invalid multistart success value')
        if row['success'] or row.get('chi2') is not None:
            _finite(row.get('chi2'), 'Multistart chi2', nonnegative=True)
        if row.get('selected') and (not row['success'] or row.get('chi2') is None
                                   or not np.isclose(row['chi2'], info['chi2'], rtol=1e-9, atol=1e-9)):
            raise ValueError('Selected restart disagrees with saved fit')
    if attempts and sum(row.get('selected') is True for row in attempts) != 1:
        raise ValueError('Saved restarts must identify exactly one selected fit')
    profiles = info.get('profiles', {})
    if not isinstance(profiles, dict):
        raise ValueError('Invalid saved profiles')
    if profiles:
        base = _finite(profiles.get('base_chi2'), 'Profile base chi2', nonnegative=True)
        if not np.isclose(base, info['chi2'], rtol=1e-9, atol=1e-9):
            raise ValueError('Profile baseline differs from saved fit')
        for name, curve in profiles.items():
            if name in ('parameter_names', 'warnings', 'interpretation', 'base_chi2'):
                continue
            if not isinstance(curve, list):
                raise ValueError(f'Invalid profile curve {name}')
            for row in curve:
                if not isinstance(row, dict) or not isinstance(row.get('success'), bool):
                    raise ValueError(f'Invalid profile point {name}')
                _finite(row.get('target'), f'{name} target')
                if row['success']:
                    value = _finite(row.get('chi2'), f'{name} chi2', nonnegative=True)
                    delta = _finite(row.get('delta_chi2'), f'{name} delta chi2')
                    if not np.isclose(value - base, delta, rtol=1e-9, atol=1e-9):
                        raise ValueError(f'Inconsistent profile delta chi2 for {name}')
                elif row.get('chi2') is not None or row.get('delta_chi2') is not None:
                    raise ValueError(f'Failed profile point {name} must have null chi2/delta')
    bootstrap = info.get('bootstrap', {})
    if not isinstance(bootstrap, dict):
        raise ValueError('Invalid saved bootstrap summary')
    if bootstrap:
        samples = bootstrap.get('samples', [])
        if not isinstance(samples, list) or not all(isinstance(row, dict) for row in samples):
            raise ValueError('Invalid bootstrap samples')
        names = bootstrap.get('parameter_names', [])
        if (not isinstance(names, list) or not names or len(set(names)) != len(names)
                or not all(isinstance(name, str) and name in info['parameters'] for name in names)):
            raise ValueError('Invalid bootstrap parameter names')
        for index, row in enumerate(samples):
            if isinstance(row.get('index'), bool) or row.get('index') != index:
                raise ValueError('Bootstrap indices must be an ordered sequence from zero')
            if not isinstance(row.get('success'), bool):
                raise ValueError('Invalid bootstrap success value')
            if row['success'] or row.get('chi2') is not None:
                _finite(row.get('chi2'), 'Bootstrap chi2', nonnegative=True)
            if row['success']:
                values, errors = row.get('parameters'), row.get('stderr')
                if (not isinstance(values, list) or len(values) != len(names)
                        or not isinstance(errors, list) or len(errors) != len(names)):
                    raise ValueError('Bootstrap parameter/error vectors do not match names')
                for value in values:
                    _finite(value, 'Bootstrap parameter')
                for value in errors:
                    if value is not None:
                        _finite(value, 'Bootstrap stderr', nonnegative=True)
        if 'replicates' in bootstrap and bootstrap['replicates'] != len(samples):
            raise ValueError('Bootstrap replicate count differs from sample records')
        if 'successful' in bootstrap and bootstrap['successful'] != sum(row['success'] for row in samples):
            raise ValueError('Bootstrap success count differs from sample records')
        intervals = bootstrap.get('intervals', {})
        if not isinstance(intervals, dict):
            raise ValueError('Invalid bootstrap intervals')
        for name, interval in intervals.items():
            if not isinstance(interval, dict) or not isinstance(interval.get('fixed'), bool):
                raise ValueError(f'Invalid bootstrap interval {name}')
            if interval.get('n_success') != sum(row['success'] for row in samples):
                raise ValueError(f'Bootstrap interval {name} has incorrect successful count')
            if interval['n_success']:
                bounds = [_finite(interval.get(key), f'{name} interval {key}')
                          for key in ('lower', 'median', 'upper')]
                if bounds != sorted(bounds):
                    raise ValueError(f'Bootstrap interval {name} is not ordered')
            elif any(interval.get(key) is not None for key in ('lower', 'median', 'upper')):
                raise ValueError(f'Bootstrap interval {name} requires successful replicates')


def _residual_stats(rows):
    values = np.asarray([row['residual_sigma'] for row in rows])
    return {'n_points': len(rows), 'chi2': float(values @ values),
            'residual_mean': float(values.mean()),
            'residual_rms': float(np.sqrt(np.mean(values**2)))}


def _report_summary(info, rows, result_path, predictions_path):
    summary = {key: info[key] for key in
               ('method', 'chi2', 'dof', 'n_points', 'n_parameters', 'kex', 'pB',
                'derived_se', 'parameters', 'parameter_order', 'covariance_note',
                'jacobian_rank', 'scaled_condition', 'at_bounds', 'warnings',
                'multistart', 'profiles', 'bootstrap', 'provenance') if key in info}
    summary.update(schema_version=1,
                   interpretation='Saved-fit report; no optimization performed. Residuals are (observed - predicted) / supplied sigma.',
                   inputs={label: {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
                           for label, path in (('result_json', result_path), ('predictions_csv', predictions_path))},
                   overall=_residual_stats(rows), by_residue={}, by_dataset={}, by_residue_dataset=[])
    for residue in dict.fromkeys(row['residue'] for row in rows):
        selected = [row for row in rows if row['residue'] == residue]
        summary['by_residue'][residue] = _residual_stats(selected)
        for dataset in dict.fromkeys(row['dataset_index'] for row in selected):
            points = [row for row in selected if row['dataset_index'] == dataset]
            summary['by_residue_dataset'].append({'residue': residue, 'dataset_index': dataset,
                                                  **_residual_stats(points)})
    for dataset in dict.fromkeys(row['dataset_index'] for row in rows):
        summary['by_dataset'][str(dataset)] = _residual_stats([row for row in rows if row['dataset_index'] == dataset])
    return summary


def _summary_lines(summary):
    def number(value):
        return 'unavailable' if value is None else f'{value:.8g}'

    lines = ['Sideband saved-fit report', summary['interpretation'], '']
    for key, entry in summary['inputs'].items():
        lines.extend([f'{key}: {entry["path"]}', f'SHA256: {entry["sha256"]}'])
    lines.extend(['', 'Standardized residual statistics: n, chi2, mean, RMS'])
    groups = [('all', summary['overall'])]
    groups += [(f'residue {key}', value) for key, value in summary['by_residue'].items()]
    groups += [(f'dataset {key}', value) for key, value in summary['by_dataset'].items()]
    groups += [(f'{row["residue"]} / dataset {row["dataset_index"]}', row)
               for row in summary['by_residue_dataset']]
    for label, values in groups:
        lines.append(f'{label}: {values["n_points"]}, {number(values["chi2"])}, '
                     f'{number(values["residual_mean"])}, {number(values["residual_rms"])}')
    lines.extend(['', f'dof: {summary.get("dof", "unavailable")}',
                  'Parameters: value +/- local standard error (supplied absolute sigma)'])
    for name, item in summary['parameters'].items():
        state = 'varied' if item['vary'] else 'fixed assumption'
        lines.append(f'{name}: {number(item["value"])} +/- {number(item.get("stderr"))} [{state}]')
    for name in ('kex', 'pB'):
        if name in summary:
            lines.append(f'{name}: {number(summary[name])} +/- {number(summary.get("derived_se", {}).get(name))}')
    lines.append(summary.get('covariance_note', 'Local standard errors are conditional on the fitted model and fixed inputs.'))
    if 'multistart' in summary:
        lines.extend(['', 'All restart attempts:'])
        for row in summary['multistart']:
            lines.append(f'#{row.get("index")}: success={row["success"]}, selected={row.get("selected", False)}, '
                         f'chi2={number(row.get("chi2"))}, status={row.get("status")}; {row.get("message", "")}')
    for name, curve in summary.get('profiles', {}).items():
        if name in ('parameter_names', 'warnings', 'interpretation', 'base_chi2'):
            continue
        lines.extend(['', f'Likelihood profile {name}: (failed points retained; negative deltas retained)'])
        for row in curve:
            lines.append(f'{number(row["target"])}: success={row["success"]}, chi2={number(row.get("chi2"))}, '
                         f'delta={number(row.get("delta_chi2"))}, status={row.get("status")}; {row.get("message", "")}')
    bootstrap = summary.get('bootstrap', {})
    if bootstrap:
        lines.extend(['', f'Parametric bootstrap: {bootstrap.get("successful")} / {bootstrap.get("replicates")} successful, '
                      f'seed={bootstrap.get("seed")}, confidence={bootstrap.get("confidence")}',
                      bootstrap.get('interpretation', 'Intervals are conditional on the selected model and fixed inputs.')])
        for name, interval in bootstrap.get('intervals', {}).items():
            lines.append(f'{name}: [{number(interval.get("lower"))}, {number(interval.get("upper"))}], '
                         f'median={number(interval.get("median"))}, n={interval.get("n_success")}'
                         + (' [fixed assumption]' if interval.get('fixed') else ''))
        lines.append('All bootstrap replicates:')
        for row in bootstrap.get('samples', []):
            lines.append(f'#{row.get("index")}: success={row["success"]}, chi2={number(row.get("chi2"))}; {row.get("message", "")}')
    warnings = list(summary.get('warnings', [])) + bootstrap.get('warnings', []) + summary.get('profiles', {}).get('warnings', [])
    if warnings:
        lines.extend(['', 'Warnings:', *dict.fromkeys(warnings)])
    return lines


def _summary_pdf_pages(pdf, summary):
    import matplotlib.pyplot as plt

    lines = [part for line in _summary_lines(summary)
             for part in (textwrap.wrap(line, width=100, replace_whitespace=False) or [''])]
    for start in range(0, len(lines), 54):
        fig = plt.figure(figsize=(8.5, 11))
        fig.text(.07, .95, '\n'.join(lines[start:start + 54]), va='top',
                 family='monospace', fontsize=8.5, linespacing=1.5)
        pdf.savefig(fig)
        plt.close(fig)
    for name, curve in summary.get('profiles', {}).items():
        if name in ('parameter_names', 'warnings', 'interpretation', 'base_chi2'):
            continue
        fig, ax = plt.subplots(figsize=(8.5, 6))
        fig.subplots_adjust(left=.12, right=.97, bottom=.18, top=.89)
        # NaN gaps prevent the line from bridging a failed nuisance fit.
        y = [row['delta_chi2'] if row['success'] else np.nan for row in curve]
        ax.plot([row['target'] for row in curve], y, 'o-')
        failures = [row['target'] for row in curve if not row['success']]
        if failures:
            ax.plot(failures, [.03] * len(failures), 'rx', transform=ax.get_xaxis_transform(),
                    label='failed target (gap)', clip_on=False)
            ax.legend()
        ax.axhline(0, color='black', linewidth=.7)
        ax.set(xlabel=name, ylabel='Delta chi-squared', title=f'Likelihood profile: {name}')
        ax.grid(alpha=.25)
        fig.text(.5, .018, 'Within-model likelihood differences; no automatic confidence interval.',
                 ha='center', fontsize=8)
        pdf.savefig(fig)
        plt.close(fig)


def regenerate_report(result_path, out_prefix, predictions_path=None):
    """Write a reconciled JSON/text/PDF report from saved results, without fitting.

    Paths are returned as ``summary_json``, ``summary_txt`` and ``pdf``. Without
    an explicit CSV path, ``NAME_result.json`` pairs with ``NAME_predictions.csv``.
    Existing files and dangling symlinks are rejected before any output is made.
    """
    result_path = Path(result_path).expanduser().absolute()
    if predictions_path is None:
        if not result_path.name.endswith('_result.json'):
            raise ValueError('Result must end in _result.json or provide --predictions')
        predictions_path = result_path.with_name(result_path.name[:-len('_result.json')] + '_predictions.csv')
    predictions_path = Path(predictions_path).expanduser().absolute()
    prefix = Path(out_prefix).expanduser().absolute()
    paths = {name: Path(str(prefix) + suffix) for name, suffix in
             (('summary_json', '_summary.json'), ('summary_txt', '_summary.txt'), ('pdf', '.pdf'))}
    for path in paths.values():
        if path.exists() or path.is_symlink():
            raise FileExistsError(f'Report output already exists: {path}')
    info, rows = _load_report_inputs(result_path, predictions_path)
    summary = _report_summary(info, rows, result_path, predictions_path)
    text = '\n'.join(_summary_lines(summary)) + '\n'
    encoded = json.dumps(summary, indent=2, allow_nan=False) + '\n'
    buffer = BytesIO()
    fit_pdf(buffer, rows, summary)
    contents = {'summary_json': encoded.encode(), 'summary_txt': text.encode(), 'pdf': buffer.getvalue()}
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for name, path in paths.items():
        with path.open('xb') as stream:
            stream.write(contents[name])
    return {name: str(path) for name, path in paths.items()}
