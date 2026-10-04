"""Experimental design: expected uncertainty of a planned Sideband measurement.

A design file names a base configuration (decoupling, residues, RF and proton
modes, bounds, fixed parameters), the generating truth and one or more
acquisition scenarios. For each scenario noise-free synthetic inputs are
written, the model is evaluated at the truth and the local Fisher information
gives the expected standard errors, correlations, rank and condition. No
optimizer runs. Expected errors are local linear statements at the truth with
the supplied absolute sigma; they rank designs and do not validate a sample.
"""
import contextlib
import copy
import hashlib
import json
import math
import numbers
from pathlib import Path

import numpy as np

from run import load_config
from sb_analysis import identifiability
from sb_report import provenance

DATASET_KEYS = ('v1n_hz', 'v1err_hz', 'T', 'sigma', 'offsets_ppm', 'offsets_rel_ppm',
                'field_mhz', 'decoupling')


def _finite(value, label, *, positive=False, nonnegative=False):
    if (isinstance(value, bool) or not isinstance(value, numbers.Real) or not np.isfinite(value)
            or (positive and value <= 0) or (nonnegative and value < 0)):
        raise ValueError(f'{label} must be a finite' + (' positive' if positive else ' nonnegative' if nonnegative else '') + ' number')
    return float(value)


def _offsets(spec, label):
    if isinstance(spec, dict):
        if set(spec) != {'min', 'max', 'n'}:
            raise ValueError(f'{label} grid needs exactly min, max and n')
        lo, hi = _finite(spec['min'], f'{label}.min'), _finite(spec['max'], f'{label}.max')
        n = spec['n']
        if isinstance(n, bool) or not isinstance(n, int) or n < 2 or lo >= hi:
            raise ValueError(f'{label} grid needs min < max and an integer n >= 2')
        return np.linspace(lo, hi, n)
    if not isinstance(spec, list) or len(spec) < 2:
        raise ValueError(f'{label} must be a grid object or a list of at least two values')
    values = np.array([_finite(v, label) for v in spec], dtype=float)
    if len(np.unique(values)) != len(values):
        raise ValueError(f'{label} contains repeated offsets')
    return values


def load_design(path):
    path = Path(path).expanduser().resolve()
    design = json.loads(path.read_text(encoding='utf-8'), parse_constant=lambda v: (_ for _ in ()).throw(ValueError(f'Non-finite JSON value {v}')))
    if not isinstance(design, dict) or set(design) - {'config', 'truth', 'truth_result', 'scenarios', 'optimize'}:
        raise ValueError('Design allows config, truth or truth_result, scenarios and optimize')
    if not isinstance(design.get('config'), str) or not design['config']:
        raise ValueError('Design needs a base config path')
    if ('truth' in design) == ('truth_result' in design):
        raise ValueError('Supply exactly one of truth or truth_result')
    scenarios = design.get('scenarios')
    if not isinstance(scenarios, list) or not scenarios:
        raise ValueError('Design needs a nonempty scenarios list')
    names = []
    for scenario in scenarios:
        if (not isinstance(scenario, dict) or set(scenario) - {'name', 'datasets'}
                or not isinstance(scenario.get('name'), str) or not scenario['name'].strip()
                or not isinstance(scenario.get('datasets'), list) or not scenario['datasets']):
            raise ValueError('Each scenario needs a name and a nonempty datasets list')
        if not all(c.isalnum() or c in '-_' for c in scenario['name']):
            raise ValueError('Scenario names may contain letters, digits, - and _ only')
        names.append(scenario['name'])
        for i, ds in enumerate(scenario['datasets']):
            label = f'{scenario["name"]}.datasets[{i}]'
            if not isinstance(ds, dict) or set(ds) - set(DATASET_KEYS):
                raise ValueError(f'{label} allows {", ".join(DATASET_KEYS)}')
            _finite(ds.get('v1n_hz'), f'{label}.v1n_hz', positive=True)
            _finite(ds.get('v1err_hz', 0.), f'{label}.v1err_hz', nonnegative=True)
            _finite(ds.get('T'), f'{label}.T', positive=True)
            _finite(ds.get('sigma'), f'{label}.sigma', positive=True)
            if 'field_mhz' in ds:
                _finite(ds['field_mhz'], f'{label}.field_mhz', positive=True)
            if ('offsets_ppm' in ds) == ('offsets_rel_ppm' in ds):
                raise ValueError(f'{label} needs exactly one of offsets_ppm or offsets_rel_ppm')
            _offsets(ds.get('offsets_ppm', ds.get('offsets_rel_ppm')), f'{label}.offsets')
            if 'decoupling' in ds and not isinstance(ds['decoupling'], dict):
                raise ValueError(f'{label}.decoupling must be an object')
    if len(set(names)) != len(names):
        raise ValueError('Scenario names must be unique')
    optimize = design.get('optimize')
    if optimize is not None:
        allowed = {'scenario', 'budget', 'criterion', 'min_per_dataset'}
        if not isinstance(optimize, dict) or set(optimize) - allowed or 'scenario' not in optimize or 'budget' not in optimize:
            raise ValueError('optimize needs scenario and budget and allows criterion and min_per_dataset')
        if optimize['scenario'] not in names:
            raise ValueError('optimize.scenario must name one of the scenarios')
        budget = optimize['budget']
        if isinstance(budget, bool) or not isinstance(budget, int) or budget < 2:
            raise ValueError('optimize.budget must be an integer of at least 2')
        minimum = optimize.get('min_per_dataset', 2)
        if isinstance(minimum, bool) or not isinstance(minimum, int) or minimum < 0:
            raise ValueError('optimize.min_per_dataset must be a nonnegative integer')
        criterion = optimize.get('criterion', 'kex')
        if not isinstance(criterion, str) or not criterion:
            raise ValueError('optimize.criterion must be kex, pB, D or a free parameter name')
        candidate = next(s for s in scenarios if s['name'] == optimize['scenario'])
        if any('offsets_ppm' not in ds for ds in candidate['datasets']):
            raise ValueError('optimize requires absolute offsets_ppm grids in the candidate scenario (one offset list per spectrum)')
        if budget < minimum * len(candidate['datasets']):
            raise ValueError('optimize.budget is smaller than min_per_dataset times the number of datasets')
    config_path = (path.parent / design['config']).resolve()
    config = load_config(str(config_path))
    if not str(config['init'].get('Method', '')).startswith('Sideband'):
        raise ValueError('Design requires a Sideband base configuration')
    if 'truth_result' in design:
        result = json.loads((path.parent / design['truth_result']).resolve().read_text(encoding='utf-8'))
        if not isinstance(result, dict) or not isinstance(result.get('parameters'), dict):
            raise ValueError('truth_result must be a saved Sideband result JSON with parameters')
        truth = {name: item.get('value') for name, item in result['parameters'].items()}
    else:
        truth = design['truth']
    if not isinstance(truth, dict) or not truth or not all(_finite(v, f'truth.{k}') is not None for k, v in truth.items()):
        raise ValueError('truth must map parameter names to finite numbers')
    return {'path': path, 'config': config, 'config_path': config_path, 'truth': dict(truth),
            'scenarios': copy.deepcopy(scenarios), 'optimize': copy.deepcopy(optimize)}


def _write_dataset(path, model, truth_vector, index, offsets_spec, sigma, field, T, v1, v1err):
    """Write one noise-free synthetic input for all active residues; return point count."""
    names = list(model.parameter_names)
    parameters = model.seParam(truth_vector)
    lines = [f'{field:.12g}', f'{T:.12g}', f'{v1:.12g} {v1err:.12g}', '# offset(ppm) intensity error']
    count = 0
    for i, residue in enumerate(model.dataset.res):
        if not residue.active:
            continue
        peak = truth_vector[names.index(f'{residue.label}.peak_ppm')]
        dw = truth_vector[names.index(f'{residue.label}.dw_ppm')]
        r2a = truth_vector[names.index(f'{residue.label}.R2a')]
        r2b = truth_vector[names.index(f'{residue.label}.R2b')]
        grid = _offsets(offsets_spec[1], 'offsets')
        offsets = grid if offsets_spec[0] == 'offsets_ppm' else peak + grid
        es = next(es for es in residue.estSpecs if es.dataset_index == index)
        values = np.asarray(model.calc(parameters, i, offsets, es), dtype=float)
        if not np.isfinite(values).all():
            raise ValueError(f'Non-finite prediction for {residue.label} in scenario dataset {index}')
        lines.append(f'# {residue.label} R2a: {r2a:.12g} R2b: {r2b:.12g} dw: {dw:.12g}')
        lines.extend(f'{o:.15g} {v:.15g} {sigma:.12g}' for o, v in zip(offsets, values))
        count += len(offsets)
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return count


def _placeholder_dataset(path, config, index, scenario_ds, truth):
    """Write a placeholder input so the model can be built before predictions exist."""
    field = scenario_ds.get('field_mhz')
    if field is None:
        template = Path(config['datasets'][0])
        field = float(template.read_text(encoding='utf-8').splitlines()[0].split()[0])
    lines = [f'{field:.12g}', f'{scenario_ds["T"]:.12g}',
             f'{scenario_ds["v1n_hz"]:.12g} {scenario_ds.get("v1err_hz", 0.):.12g}', '# offset(ppm) intensity error']
    kind = 'offsets_ppm' if 'offsets_ppm' in scenario_ds else 'offsets_rel_ppm'
    grid = _offsets(scenario_ds[kind], 'offsets')
    for entry in config['residues']:
        if entry.get('flag') != 'on':
            continue
        label = entry['name']
        for key in ('peak_ppm', 'dw_ppm', 'R2a', 'R2b'):
            if f'{label}.{key}' not in truth:
                raise ValueError(f'truth must include {label}.{key}')
        peak = truth[f'{label}.peak_ppm']
        offsets = grid if kind == 'offsets_ppm' else peak + grid
        lines.append(f'# {label} R2a: {truth[f"{label}.R2a"]:.12g} R2b: {truth[f"{label}.R2b"]:.12g} dw: {truth[f"{label}.dw_ppm"]:.12g}')
        lines.extend(f'{o:.15g} 1 {scenario_ds["sigma"]:.12g}' for o in offsets)
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return field, (kind, scenario_ds[kind])


def evaluate_scenario(design, scenario, folder, *, workers=1):
    """Expected errors for one scenario; writes its synthetic inputs under folder."""
    from sb_parallel import WorkerPool
    from sbfit import SidebandModel, _validate_fit_config

    config = copy.deepcopy(design['config'])
    truth = design['truth']
    folder.mkdir(parents=True, exist_ok=False)
    # Synthetic inputs contain the active residues only; inactive entries would
    # otherwise be reported as missing from the data.
    config['residues'] = [entry for entry in config['residues'] if entry.get('flag') == 'on']
    if not config['residues']:
        raise ValueError('Design needs at least one active residue in the base configuration')
    datasets, specs = [], []
    for i, ds in enumerate(scenario['datasets']):
        path = folder / f'data_{i}.txt'
        field, spec = _placeholder_dataset(path, config, i, ds, truth)
        datasets.append(str(path))
        specs.append((field, spec, ds))
    config['datasets'] = datasets
    config['Project Name'] = str(folder / 'design')
    init = {key: value for key, value in config['init'].items()
            if key not in ('initial', 'multistart', 'profile', 'bootstrap', 'profile_interval')}
    config['init'] = init
    sb = config['sideband']
    overrides = sb.get('datasets')
    base_override = overrides[0] if isinstance(overrides, list) and overrides else {}
    sb['datasets'] = [{**base_override, **ds.get('decoupling', {})} for ds in scenario['datasets']]
    model = SidebandModel(config, design['config_path'].parent)
    _validate_fit_config(model, init)
    names = list(model.parameter_names)
    missing = [name for name in names if name not in truth]
    extra = [name for name in truth if name not in names]
    if missing or extra:
        raise ValueError(f'truth must map exactly the model parameters; missing {missing}, unknown {extra}')
    vector = np.array([truth[name] for name in names], dtype=float)
    model.prepare_fit(p0=vector, fitting_config=init)
    points = 0
    for i, (field, spec, ds) in enumerate(specs):
        points += _write_dataset(folder / f'data_{i}.txt', model, vector, i, spec, ds['sigma'],
                                 field, ds['T'], ds['v1n_hz'], ds.get('v1err_hz', 0.))
    # Rebuild on the noise-free predictions so the saved inputs and the evaluation agree.
    model = SidebandModel(config, design['config_path'].parent)
    model.prepare_fit(p0=vector, fitting_config=init)
    pool_context = (WorkerPool(workers, config, design['config_path'].parent) if workers > 1
                    else contextlib.nullcontext())
    with pool_context as pool:
        report = identifiability(model, vector, pool=pool)
    saturation = sum(len(es.offset) * es.T for r in model.dataset.res if r.active for es in r.estSpecs)
    (folder / 'design_config.json').write_text(json.dumps(config, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    return {
        'name': scenario['name'], 'n_points': points, 'n_datasets': len(datasets),
        'total_saturation_s': float(saturation),
        'datasets': [{'path': path, 'sha256': hashlib.sha256(Path(path).read_bytes()).hexdigest(),
                      **{k: v for k, v in ds.items() if k != 'decoupling'}}
                     for path, ds in zip(datasets, scenario['datasets'])],
        'identifiability': report,
    }


def run_design(design_path, out, *, workers=1):
    from sb_parallel import validate_workers

    workers = validate_workers(workers)
    out = Path(out).expanduser().absolute()
    if out.exists() or out.is_symlink():
        raise FileExistsError(f'Design output already exists: {out}')
    design = load_design(design_path)
    out.mkdir(parents=True)
    results = [evaluate_scenario(design, scenario, out / 'scenarios' / scenario['name'], workers=workers)
               for scenario in design['scenarios']]
    optimization = None
    if design.get('optimize'):
        optimization, extra = optimize_design(design, out, results, workers=workers)
        results.extend(extra)
    summary = {
        'schema_version': 1, 'design_path': str(design['path']),
        'design_sha256': hashlib.sha256(design['path'].read_bytes()).hexdigest(),
        'config_path': str(design['config_path']), 'truth': design['truth'],
        'confidence_note': 'Expected standard errors are local linear values at the truth with the supplied absolute sigma; no optimizer ran.',
        'scenarios': results,
        'optimization': optimization,
        'provenance': provenance(design['config'], design['config_path'].parent),
    }
    text = '\n'.join(design_lines(summary)) + '\n'
    (out / 'design.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    (out / 'design.txt').write_text(text, encoding='utf-8')
    design_pdf(out / 'design.pdf', summary)
    return {'design_json': str(out / 'design.json'), 'design_txt': str(out / 'design.txt'),
            'design_pdf': str(out / 'design.pdf')}


def design_lines(summary):
    def number(value):
        return 'unavailable' if value is None else f'{value:.6g}'

    lines = ['Sideband experimental design: expected local errors at the truth',
             summary['confidence_note'], f'Base configuration: {summary["config_path"]}', '']
    lines.append('scenario | points | datasets | sum(points*T) s | rank/free | condition | SE kex | SE pB | weak | strong pairs')
    for row in summary['scenarios']:
        rep = row['identifiability']
        lines.append(' | '.join([
            row['name'], str(row['n_points']), str(row['n_datasets']), number(row['total_saturation_s']),
            f'{rep["jacobian_rank"]}/{rep["n_free"]}', number(rep['scaled_condition']),
            number(rep['derived_se']['kex']), number(rep['derived_se']['pB']),
            ','.join(rep['weak_parameters']) or '-', str(len(rep['strong_correlations']))]))
    for row in summary['scenarios']:
        rep = row['identifiability']
        lines.extend(['', f'Scenario {row["name"]}: expected SE (relative) per free parameter'])
        for name in rep['free_parameters']:
            se, rel = rep['expected_se'][name], rep['relative_se'][name]
            shown = 'n/a' if rel is None else '>1000%' if rel > 10 else f'{100*rel:.2f}%'
            lines.append(f'{name}: {number(se)} ({shown})')
        if rep['strong_correlations']:
            lines.append('Strong correlations: ' + '; '.join(
                f'{"/".join(c["parameters"])} {c["correlation"]:+.3f}' for c in rep['strong_correlations'][:8]))
        for ds in row['datasets']:
            spec = ds.get('offsets_ppm', ds.get('offsets_rel_ppm'))
            grid = f'{spec["min"]}..{spec["max"]} ({spec["n"]})' if isinstance(spec, dict) else f'{len(spec)} explicit values'
            lines.append(f'  dataset: v1n {ds["v1n_hz"]} Hz, T {ds["T"]} s, sigma {ds["sigma"]}, offsets {grid}')
    opt = summary.get('optimization')
    if opt:
        lines.extend(['', f'Optimized design: criterion {opt["criterion"]}, budget {opt["budget"]} of '
                      f'{opt["n_candidate_measurements"]} candidate spectrum rows (min {opt["min_per_dataset"]} per dataset)',
                      f'criterion value: candidate {number(opt["candidate_criterion"])} -> optimized {number(opt["optimized_criterion"])}',
                      opt['interpretation']])
        for dataset, offsets in opt['selected_offsets_ppm'].items():
            lines.append(f'  dataset {dataset}: ' + ', '.join(f'{o:.4g}' for o in offsets))
    return lines


def design_pdf(path, summary):
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    from sb_report import _paginate
    import textwrap

    names = [row['name'] for row in summary['scenarios']]
    with PdfPages(path) as pdf:
        fig, axes = plt.subplots(1, 2, figsize=(11, 5))
        for ax, key in zip(axes, ('kex', 'pB')):
            values = [row['identifiability']['derived_se'][key] for row in summary['scenarios']]
            heights = [v if v is not None else 0. for v in values]
            bars = ax.bar(range(len(names)), heights, color='#4477aa')
            for bar, v in zip(bars, values):
                if v is None:
                    ax.text(bar.get_x() + bar.get_width() / 2, 0, 'n/a', ha='center', va='bottom', fontsize=8)
            ax.set_xticks(range(len(names)))
            ax.set_xticklabels(names, rotation=30, ha='right', fontsize=8)
            ax.set_ylabel(f'expected SE of {key}')
            ax.set_title(f'Expected local SE of {key} at the truth')
            ax.grid(alpha=.25, axis='y')
        fig.text(.5, .01, summary['confidence_note'], ha='center', fontsize=7)
        fig.tight_layout(rect=(0, .04, 1, 1))
        pdf.savefig(fig)
        plt.close(fig)
        opt = summary.get('optimization')
        if opt:
            fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
            path = opt['criterion_path']
            axes[0].plot([row['n_measurements'] for row in path], [row['criterion'] for row in path], '-')
            axes[0].set(xlabel='spectrum rows kept', ylabel=f'criterion ({opt["criterion"]})',
                        title='Backward elimination path')
            axes[0].invert_xaxis()
            axes[0].grid(alpha=.25)
            for dataset, offsets in opt['selected_offsets_ppm'].items():
                axes[1].plot(offsets, [int(dataset)] * len(offsets), '|', markersize=14, label=f'dataset {dataset}')
            axes[1].set(xlabel='saturation offset (ppm)', ylabel='dataset index', title='Selected offsets')
            axes[1].grid(alpha=.25)
            axes[1].legend(fontsize=8)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
        lines = [part for line in design_lines(summary)
                 for part in (textwrap.wrap(line, width=110, replace_whitespace=False) or [''])]
        for page in _paginate(lines, 56):
            fig = plt.figure(figsize=(8.5, 11))
            fig.text(.05, .95, '\n'.join(page), va='top', family='monospace', fontsize=7.5, linespacing=1.5)
            pdf.savefig(fig)
            plt.close(fig)


def _criterion_gradient(model, vector, criterion):
    """Gradient of the criterion quantity with respect to the free parameters, or None for D."""
    names = list(model.parameter_names)
    free = [int(i) for i in model.free]
    if criterion == 'D':
        return None
    g = np.zeros(len(names))
    if criterion in ('kex', 'pB'):
        if 'kbc' in names:
            raise ValueError('kex/pB criteria are defined for the two-state model only')
        a, b = names.index('kab'), names.index('kba')
        if a not in free or b not in free:
            raise ValueError('kex/pB criteria need kab and kba to vary')
        total = vector[a] + vector[b]
        if criterion == 'kex':
            g[a] = g[b] = 1.0
        else:
            g[a], g[b] = vector[b] / total**2, -vector[a] / total**2
    else:
        if criterion not in names:
            raise ValueError(f'Unknown criterion parameter {criterion}')
        if names.index(criterion) not in free:
            raise ValueError(f'Criterion parameter {criterion} is fixed')
        g[names.index(criterion)] = 1.0
    return g[free]


def _measurement_groups(model):
    """Rows of the residual vector grouped by (dataset, offset index): one spectrum row each."""
    groups = {}
    position = 0
    for i, es, offsets, _, _ in model._prepare_data():
        for k in range(len(offsets)):
            groups.setdefault((es.dataset_index, k), []).append(position + k)
        position += len(offsets)
    return groups


def select_measurements(jacobian, groups, budget, gradient, *, min_per_dataset=2):
    """Backward elimination of spectrum rows that least worsen the criterion.

    ``groups`` maps (dataset, offset index) to the Jacobian rows of that spectrum
    row (one per residue). The criterion is the variance of ``gradient``-weighted
    parameters, or the D-criterion (log det of the information) when gradient is
    None. Removal uses exact Woodbury downdates and refuses to make the
    information singular. Returns the kept keys and the criterion path.
    """
    J = np.asarray(jacobian, dtype=float)
    keys = list(groups)
    rows = [np.asarray(groups[key], dtype=int) for key in keys]
    active = {key: True for key in keys}
    per_dataset = {}
    for key in keys:
        per_dataset[key[0]] = per_dataset.get(key[0], 0) + 1
    M = J.T @ J
    try:
        Minv = np.linalg.inv(M)
    except np.linalg.LinAlgError as exc:
        raise ValueError('Candidate design has singular information; add points or fix parameters') from exc
    sign, logdet = np.linalg.slogdet(M)
    if sign <= 0:
        raise ValueError('Candidate design information is not positive definite')

    def value():
        return float(gradient @ Minv @ gradient) if gradient is not None else float(-logdet)

    path = [{'n_measurements': len(keys), 'criterion': value()}]
    count = len(keys)
    while count > budget:
        best = None
        for index, key in enumerate(keys):
            if not active[key] or per_dataset[key[0]] <= min_per_dataset:
                continue
            Jg = J[rows[index]]
            U = Jg @ Minv
            H = U @ Jg.T
            I_H = np.eye(len(Jg)) - H
            det = np.linalg.det(I_H)
            if det <= 1e-12:
                continue
            if gradient is None:
                delta = -math.log(det)
            else:
                a = U @ gradient
                delta = float(a @ np.linalg.solve(I_H, a))
            if best is None or delta < best[0]:
                best = (delta, index, Jg, U, I_H, det)
        if best is None:
            raise ValueError('No removable measurement keeps the information nonsingular and the per-dataset minimum')
        delta, index, Jg, U, I_H, det = best
        Minv = Minv + U.T @ np.linalg.solve(I_H, U)
        logdet += math.log(det)
        active[keys[index]] = False
        per_dataset[keys[index][0]] -= 1
        count -= 1
        path.append({'n_measurements': count, 'criterion': value()})
    kept = [key for key in keys if active[key]]
    return kept, path


def optimize_design(design, out, results, *, workers=1):
    """Select the measurement budget from the candidate scenario and evaluate it."""
    from sb_analysis import local_jacobian
    from sb_parallel import WorkerPool
    from sbfit import SidebandModel

    settings = design['optimize']
    name = settings['scenario']
    candidate_row = next(row for row in results if row['name'] == name)
    folder = out / 'scenarios' / name
    config = json.loads((folder / 'design_config.json').read_text(encoding='utf-8'))
    model = SidebandModel(config, design['config_path'].parent)
    vector = np.array([design['truth'][n] for n in model.parameter_names], dtype=float)
    model.prepare_fit(p0=vector, fitting_config=config['init'])
    model._fit_data = model._prepare_data()
    criterion = settings.get('criterion', 'kex')
    gradient = _criterion_gradient(model, vector, criterion)
    pool_context = (WorkerPool(workers, config, design['config_path'].parent) if workers > 1
                    else contextlib.nullcontext())
    with pool_context as pool:
        jacobian = local_jacobian(model, vector, pool=pool)
    groups = _measurement_groups(model)
    kept, path = select_measurements(jacobian, groups, settings['budget'], gradient,
                                     min_per_dataset=settings.get('min_per_dataset', 2))
    candidate = next(s for s in design['scenarios'] if s['name'] == name)
    grids = [_offsets(ds['offsets_ppm'], 'offsets') for ds in candidate['datasets']]
    selected = {}
    for dataset, k in kept:
        selected.setdefault(dataset, []).append(float(grids[dataset][k]))
    optimized = {'name': f'{name}_optimized_{settings["budget"]}',
                 'datasets': [dict(ds, offsets_ppm=sorted(selected.get(i, []))) for i, ds in enumerate(candidate['datasets'])]}
    for ds in optimized['datasets']:
        if len(ds['offsets_ppm']) < 2:
            ds['offsets_ppm'] = sorted(set(ds['offsets_ppm']) | set(grids[optimized['datasets'].index(ds)][:2].tolist()))
    # Uniform comparison: the same budget spread evenly over the candidate grids.
    total = sum(len(g) for g in grids)
    shares = [max(2, round(settings['budget'] * len(g) / total)) for g in grids]
    uniform = {'name': f'{name}_uniform_{settings["budget"]}',
               'datasets': [dict(ds, offsets_ppm=np.linspace(g[0], g[-1], share).tolist())
                            for ds, g, share in zip(candidate['datasets'], grids, shares)]}
    evaluated = [evaluate_scenario(design, scenario, out / 'scenarios' / scenario['name'], workers=workers)
                 for scenario in (optimized, uniform)]
    return {
        'criterion': criterion, 'budget': settings['budget'], 'candidate_scenario': name,
        'min_per_dataset': settings.get('min_per_dataset', 2),
        'n_candidate_measurements': len(groups), 'selected_offsets_ppm': {str(k): v for k, v in sorted(selected.items())},
        'criterion_path': path,
        'optimized_scenario': optimized['name'], 'uniform_scenario': uniform['name'],
        'candidate_criterion': path[0]['criterion'], 'optimized_criterion': path[-1]['criterion'],
        'interpretation': ('Backward elimination of whole spectrum rows (all residues at one offset) that least '
                           'increase the criterion, using exact information downdates at the truth. The selection '
                           'is local to the truth and the model; compare the optimized and uniform scenarios.'),
        'candidate_expected_se': candidate_row['identifiability']['derived_se'],
    }, evaluated
