#!/usr/bin/env python3
"""Seeded synthetic coverage of Sideband local standard-error intervals.

This evaluates repeated-data local-SE intervals under explicit generating truth,
the configured fixed inputs and supplied absolute sigma. It does not measure
coverage of bootstrap percentile intervals or validate experimental accuracy.
"""
import argparse
import hashlib
import json
import math
import numbers
from pathlib import Path
from statistics import NormalDist
import time

import numpy as np

from run import load_config
from sb_bootstrap import bootstrap_fit, validate_bootstrap
from sb_report import provenance
from sbfit import SidebandModel


def _finite_number(value):
    return isinstance(value, numbers.Real) and not isinstance(value, (bool, np.bool_)) and math.isfinite(value)


def _wilson(covered, count):
    if not count:
        return None
    z = NormalDist().inv_cdf(.975)
    fraction, z2 = covered / count, z * z
    denominator = 1 + z2 / count
    center = (fraction + z2 / (2 * count)) / denominator
    half = z * math.sqrt(fraction * (1 - fraction) / count + z2 / (4 * count**2)) / denominator
    return [max(0., center - half), min(1., center + half)]


def coverage_summary(samples, parameter_names, varying_parameters, truth, confidence=.95):
    """Count local normal intervals using positive finite SEs only.

    Missing/nonfinite/negative SEs and zero SEs are reported separately. Coverage
    and its binomial uncertainty use eligible successful fits as the denominator;
    boundary frequency uses all successful fits. Fixed quantities are excluded.
    """
    validate_bootstrap({'replicates': max(1, len(samples)), 'seed': 0, 'confidence': confidence})
    z = -NormalDist().inv_cdf((1 - confidence) / 2)
    good = [row for row in samples if row['success']]
    varied = set(varying_parameters)
    rates = varied & {'kab', 'kba'}
    total = truth['kab'] + truth['kba']
    if not math.isfinite(total) or total <= 0:
        raise ValueError('Generating exchange rate must be finite and positive')
    truths = {**truth, 'kex': total, 'pB': truth['kab'] / total}
    results = {}
    for name in [*parameter_names, 'kex', 'pB']:
        index = parameter_names.index(name) if name in parameter_names else None
        fixed = name not in varied if index is not None else not rates
        # A fixed zero kab also fixes pB at zero when kba alone varies.
        fixed = fixed or (name == 'pB' and 'kab' not in varied and truth['kab'] == 0)
        relevant = {name} if index is not None else rates
        at_bounds = sum(bool(set(row.get('at_bounds', [])) & relevant) for row in good)
        result = {'truth': float(truths[name]), 'fixed': fixed,
                  'successful': len(good), 'failed': len(samples) - len(good),
                  'finite_se': 0, 'zero_se': 0, 'unavailable_se': 0, 'eligible_se': 0,
                  'covered': 0, 'fraction': None, 'binomial_se': None, 'wilson_95': None,
                  'at_bounds': at_bounds,
                  'boundary_fraction': at_bounds / len(good) if good else None}
        results[name] = result
        if fixed:
            result['excluded_reason'] = 'Fixed assumption; no coverage claim'
            continue
        for row in good:
            estimate = row['parameters'][index] if index is not None else row[name]
            errors = row.get('stderr') if index is not None else row.get('derived_se')
            error = errors[index] if index is not None and errors is not None else (
                errors.get(name) if errors is not None else None)
            if not _finite_number(error) or error < 0:
                result['unavailable_se'] += 1
                continue
            result['finite_se'] += 1
            if error == 0:
                result['zero_se'] += 1
                continue
            if not _finite_number(estimate):
                raise ValueError('Successful sample has a nonfinite estimate')
            result['eligible_se'] += 1
            result['covered'] += int(abs(estimate - truths[name]) <= z * error)
        count = result['eligible_se']
        if count:
            fraction = result['covered'] / count
            result.update(fraction=fraction, binomial_se=math.sqrt(fraction * (1 - fraction) / count),
                          wilson_95=_wilson(result['covered'], count))
    warnings = []
    if len(samples) < 100:
        warnings.append('Fewer than 100 replicates: this is a demonstration with large sampling uncertainty.')
    if len(good) != len(samples):
        warnings.append('Coverage excludes failed fits; failures can bias coverage among successful fits.')
    if any(not row['fixed'] and row['eligible_se'] < len(good) for row in results.values()):
        warnings.append('Unavailable or zero local SEs are excluded from coverage; inspect their counts.')
    boundary_count = sum(bool(set(row.get('at_bounds', [])) & varied) for row in good)
    if boundary_count:
        warnings.append('Some successful fits reach bounds; symmetric local normal intervals can be unreliable.')
    return {'schema_version': 1, 'study': 'synthetic_local_se_coverage',
            'replicates': len(samples), 'successful': len(good), 'failed': len(samples) - len(good),
            'confidence': float(confidence), 'normal_z': z,
            'boundary_replicates': boundary_count,
            'boundary_fraction': boundary_count / len(good) if good else None,
            'parameters': results, 'warnings': warnings,
            'coverage_denominator': 'Successful refits with positive finite local standard errors for each quantity',
            'sampling_uncertainty': 'Binomial plug-in standard error and Wilson 95% interval for each coverage fraction',
            'interpretation': 'Synthetic local-SE interval coverage conditional on the generating model, fixed inputs and supplied absolute sigma; not bootstrap-interval coverage or experimental validation.',
            'fit_protocol': 'Every independent Gaussian synthetic dataset is refitted starting at the generating truth; configured multistart and profile scans are not run.'}


def _write_json(path, value):
    with path.open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(value, indent=2, allow_nan=False) + '\n')


def run_study(config_path, truth_path, replicates, seed, out, confidence=.95):
    """Validate truth, preserve inputs, and refit one fresh dataset per replicate."""
    settings = validate_bootstrap({'replicates': replicates, 'seed': seed, 'confidence': confidence})
    config_path, truth_path, out = Path(config_path).resolve(), Path(truth_path).resolve(), Path(out)
    if out.exists() or out.is_symlink():
        raise FileExistsError(f'Output directory already exists: {out}')
    config = load_config(config_path)
    model = SidebandModel(config, config_path.parent)
    initial = model.prepare_fit()
    document = json.loads(truth_path.read_text())
    truth = document.get('truth', document) if isinstance(document, dict) else document
    if (not isinstance(truth, dict) or set(truth) != set(model.parameter_names)
            or not all(_finite_number(value) for value in truth.values())):
        raise ValueError('Truth must map every model parameter name to a finite number, without extra names')
    truth = {name: float(truth[name]) for name in model.parameter_names}
    p = np.array([truth[name] for name in model.parameter_names])
    if np.any(p < model.lower) or np.any(p > model.upper):
        raise ValueError('Generating truth must lie inside the configured parameter bounds')
    free = set(int(index) for index in model.free)
    fixed = [i for i in range(len(p)) if i not in free]
    mismatches = [model.parameter_names[i] for i in fixed
                  if not np.isclose(p[i], initial[i], rtol=16*np.finfo(float).eps, atol=0.)]
    if mismatches:
        raise ValueError('Generating truth differs from fixed init parameters: ' + ', '.join(mismatches))
    metadata = provenance(config, config_path.parent)
    metadata['source_sha256'][Path(__file__).name] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    inputs = {name: {'path': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
              for name, path in (('config', config_path), ('truth', truth_path))}
    study = {'schema_version': 1, 'settings': settings, 'config': config, 'truth': truth,
             'config_directory': str(config_path.parent), 'inputs': inputs, 'provenance': metadata,
             'parameter_names': list(model.parameter_names),
             'varying_parameters': [model.parameter_names[i] for i in model.free],
             'fixed_parameters': {model.parameter_names[i]: float(initial[i]) for i in fixed},
             'data_design': 'Observed intensities are replaced by independent Gaussian draws around generating truth; offsets, RF, fields and supplied absolute sigma remain fixed.'}
    out.mkdir(parents=True, exist_ok=False)
    (out / 'samples').mkdir()
    for filename, path in (('input_config.json', config_path), ('input_truth.json', truth_path)):
        with (out / filename).open('xb') as stream:
            stream.write(path.read_bytes())
    _write_json(out / 'study.json', study)
    began = time.monotonic()
    result = bootstrap_fit(model, p, settings,
                           on_complete=lambda row: _write_json(out / 'samples' / f"{row['index']:06d}.json", row))
    summary = coverage_summary(result['samples'], model.parameter_names, result['varying_parameters'], truth, confidence)
    summary['elapsed_s'] = time.monotonic() - began
    _write_json(out / 'samples.json', result['samples'])
    _write_json(out / 'coverage.json', summary)
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True, type=Path)
    parser.add_argument('--truth', required=True, type=Path)
    parser.add_argument('--replicates', required=True, type=int)
    parser.add_argument('--seed', required=True, type=int)
    parser.add_argument('--confidence', default=.95, type=float)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    try:
        summary = run_study(args.config, args.truth, args.replicates, args.seed, args.out, args.confidence)
    except (ValueError, KeyError, TypeError, OSError, RuntimeError) as exc:
        parser.exit(1, f'Error: {exc}\n')
    print(json.dumps({'coverage': str((args.out / 'coverage.json').absolute()),
                      'replicates': summary['replicates'], 'successful': summary['successful'],
                      'failed': summary['failed']}, indent=2))


if __name__ == '__main__':
    main()
