"""Residual diagnostics: is the supplied sigma plausible and are the residuals structureless?

All uncertainty statements in SBONEST (local covariance, profile intervals,
bootstrap, AICc/BIC, F-tests) assume the supplied absolute sigma is right and
the model describes the data. These checks use only the standardized residuals
of a finished fit, so they add no fitting cost: the reduced chi-square and its
expected spread, a sigma scale estimated from the residuals, per-block reduced
chi-square, a Wald-Wolfowitz runs test and lag-1 autocorrelation of residuals
ordered by offset within each residue/dataset block, outlier counts, and a
covariance rescaled by chi2/dof as an explicitly labelled alternative. They flag
problems; they do not identify their cause.
"""
import math
from statistics import NormalDist

import numpy as np

_NORMAL = NormalDist()


def runs_test(values):
    """Two-sided Wald-Wolfowitz runs test on the signs of ``values`` (zeros dropped)."""
    signs = np.sign(np.asarray(values, dtype=float))
    signs = signs[signs != 0]
    n_pos, n_neg = int(np.sum(signs > 0)), int(np.sum(signs < 0))
    n = n_pos + n_neg
    if n_pos == 0 or n_neg == 0 or n < 4:
        return {'runs': int(n > 0), 'n_positive': n_pos, 'n_negative': n_neg,
                'expected_runs': None, 'z': None, 'p_value': None}
    runs = 1 + int(np.sum(signs[1:] != signs[:-1]))
    expected = 1 + 2 * n_pos * n_neg / n
    variance = 2 * n_pos * n_neg * (2 * n_pos * n_neg - n) / (n**2 * (n - 1))
    if variance <= 0:
        return {'runs': runs, 'n_positive': n_pos, 'n_negative': n_neg,
                'expected_runs': expected, 'z': None, 'p_value': None}
    z = (runs - expected) / math.sqrt(variance)
    p = 2 * (1 - _NORMAL.cdf(abs(z)))
    return {'runs': runs, 'n_positive': n_pos, 'n_negative': n_neg,
            'expected_runs': expected, 'z': float(z), 'p_value': float(p)}


def lag1_autocorrelation(values):
    """Lag-1 autocorrelation of a sequence, or None when it is undefined."""
    values = np.asarray(values, dtype=float)
    if values.size < 3:
        return None
    centered = values - values.mean()
    denominator = float(centered @ centered)
    if denominator <= 0:
        return None
    return float((centered[1:] @ centered[:-1]) / denominator)


def _block_summary(block):
    """Statistics for residuals of one residue/dataset block ordered by offset."""
    ordered = sorted(block, key=lambda row: row['offset_ppm'])
    values = np.asarray([row['residual_sigma'] for row in ordered], dtype=float)
    n = int(values.size)
    chi2 = float(values @ values)
    corr = lag1_autocorrelation(values)
    return {'n_points': n, 'chi2': chi2, 'reduced_chi2': chi2 / n if n else None,
            'residual_mean': float(values.mean()) if n else None,
            'runs_test': runs_test(values), 'lag1_autocorrelation': corr,
            'lag1_threshold': 2 / math.sqrt(n) if n else None,
            'max_abs_residual': float(np.max(np.abs(values))) if n else None}


def residual_diagnostics(rows, n_parameters, covariance=None, *, runs_alpha=0.01):
    """Diagnostics from prediction rows (see ``sb_report.prediction_rows``).

    ``covariance`` (full parameter order) is rescaled by chi2/dof when given; the
    rescaled errors assume a correct model with a uniformly mis-estimated sigma.
    """
    values = np.asarray([row['residual_sigma'] for row in rows], dtype=float)
    n = int(values.size)
    k = int(n_parameters)
    dof = n - k
    if n == 0 or dof <= 0:
        raise ValueError('Residual diagnostics need more points than parameters')
    chi2 = float(values @ values)
    reduced = chi2 / dof
    spread = math.sqrt(2.0 / dof)
    blocks = {}
    for row in rows:
        blocks.setdefault((row['residue'], row['dataset_index']), []).append(row)
    block_rows = [{'residue': residue, 'dataset_index': dataset, **_block_summary(block)}
                  for (residue, dataset), block in blocks.items()]
    by_residue = {}
    for row in rows:
        by_residue.setdefault(row['residue'], []).append(row['residual_sigma'])
    by_residue = {label: {'n_points': len(v), 'reduced_chi2': float(np.asarray(v) @ np.asarray(v)) / len(v)}
                  for label, v in by_residue.items()}
    by_dataset = {}
    for row in rows:
        by_dataset.setdefault(str(row['dataset_index']), []).append(row['residual_sigma'])
    by_dataset = {label: {'n_points': len(v), 'reduced_chi2': float(np.asarray(v) @ np.asarray(v)) / len(v)}
                  for label, v in by_dataset.items()}
    overall_runs = runs_test(np.concatenate([
        np.asarray([r['residual_sigma'] for r in sorted(block, key=lambda row: row['offset_ppm'])])
        for block in blocks.values()]))
    warnings = []
    if abs(reduced - 1.0) > 3 * spread:
        direction = 'larger' if reduced > 1 else 'smaller'
        warnings.append(f'Reduced chi2 {reduced:.3g} is more than three expected spreads ({spread:.3g}) from 1: '
                        f'the supplied sigma may be {"too small" if reduced > 1 else "too large"} or the model '
                        f'{"misses structure" if reduced > 1 else "is overfitting"}; residual-based sigma scale {math.sqrt(reduced):.3g} is {direction} than 1.')
    structured = [f'{row["residue"]}/dataset {row["dataset_index"]}' for row in block_rows
                  if row['runs_test']['p_value'] is not None and row['runs_test']['p_value'] < runs_alpha]
    if structured:
        warnings.append('Runs test finds systematic residual sign structure (p < '
                        f'{runs_alpha:g}) in: ' + ', '.join(structured) + '; inspect model adequacy.')
    correlated = [f'{row["residue"]}/dataset {row["dataset_index"]}' for row in block_rows
                  if row['lag1_autocorrelation'] is not None and abs(row['lag1_autocorrelation']) > row['lag1_threshold']]
    if correlated:
        warnings.append('Lag-1 residual autocorrelation beyond 2/sqrt(n) in: ' + ', '.join(correlated) + '.')
    outliers = int(np.sum(np.abs(values) > 3))
    expected_outliers = n * 2 * (1 - _NORMAL.cdf(3))
    if outliers > max(3, 3 * expected_outliers):
        warnings.append(f'{outliers} residuals beyond 3 sigma (about {expected_outliers:.1f} expected).')
    out = {
        'n_points': n, 'n_parameters': k, 'dof': dof, 'chi2': chi2,
        'reduced_chi2': reduced, 'reduced_chi2_expected_sd': spread,
        'sigma_scale_estimate': math.sqrt(reduced),
        'residual_mean': float(values.mean()), 'residual_sd': float(values.std(ddof=1)) if n > 1 else None,
        'outliers_beyond_3sigma': outliers, 'expected_outliers_beyond_3sigma': expected_outliers,
        'overall_runs_test': overall_runs,
        'by_residue': by_residue, 'by_dataset': by_dataset, 'blocks': block_rows,
        'warnings': warnings,
        'interpretation': ('Reduced chi2 near 1 (within a few expected spreads) is consistent with the supplied '
                           'absolute sigma; runs tests and autocorrelation detect ordered structure that a '
                           'correct model should not leave. Rescaled errors multiply the local covariance by '
                           'chi2/dof and assume a correct model with uniformly mis-estimated sigma.'),
    }
    if covariance is not None:
        cov = np.asarray(covariance, dtype=float)
        scaled = cov * reduced
        out['rescaled_stderr'] = [float(math.sqrt(v)) if np.isfinite(v) and v >= 0 else None for v in np.diag(scaled)]
        out['rescale_factor'] = reduced
    return out


def diagnostics_lines(report, parameter_names=None):
    """Plain-text lines summarizing a residual_diagnostics report."""
    def number(value, digits=4):
        return 'unavailable' if value is None else f'{value:.{digits}g}'

    lines = ['Residual diagnostics (supplied absolute sigma)',
             f'reduced chi2 {number(report["reduced_chi2"])} (expected 1 +/- {number(report["reduced_chi2_expected_sd"], 2)}), '
             f'sigma scale estimate {number(report["sigma_scale_estimate"])}, '
             f'residual mean {number(report["residual_mean"], 3)}, outliers beyond 3 sigma {report["outliers_beyond_3sigma"]} '
             f'(expected {number(report["expected_outliers_beyond_3sigma"], 2)})']
    overall = report['overall_runs_test']
    lines.append(f'overall runs test: runs {overall["runs"]}, expected {number(overall["expected_runs"])}, '
                 f'z {number(overall["z"], 3)}, p {number(overall["p_value"], 3)}')
    for row in report['blocks']:
        test = row['runs_test']
        lines.append(f'{row["residue"]} / dataset {row["dataset_index"]}: n {row["n_points"]}, reduced chi2 '
                     f'{number(row["reduced_chi2"])}, runs p {number(test["p_value"], 3)}, lag-1 r '
                     f'{number(row["lag1_autocorrelation"], 3)} (|r| > {number(row["lag1_threshold"], 2)} flags), '
                     f'max |residual| {number(row["max_abs_residual"], 3)}')
    if report.get('rescaled_stderr') is not None and parameter_names:
        lines.append(f'Rescaled standard errors (local SE x sqrt(chi2/dof) = x{number(math.sqrt(report["rescale_factor"]))}; '
                     'assumes a correct model and uniformly mis-estimated sigma):')
        lines.extend(f'  {name}: {number(value, 6)}' for name, value in zip(parameter_names, report['rescaled_stderr'])
                     if value is not None)
    lines.extend(report['warnings'])
    return lines
