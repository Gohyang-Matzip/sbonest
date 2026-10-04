"""Compare a shared-exchange (global) Sideband fit with per-residue individual fits.

The global model shares kab/kba across the selected residues; the individual
model fits every residue alone with its own exchange rates and is therefore the
larger, nesting model. Both are fitted with the configured settings through the
normal run_config path, so each sub-fit has its own checkpoint, result JSON and
predictions. Statistics use the supplied absolute sigma: chi2, AICc, BIC and a
nested F-test. They compare models conditional on that sigma and on the shared
RF/proton settings; they do not establish which description is physically true.
"""
import copy
import json
import math
from pathlib import Path


def _residue_config(config, labels, keep):
    """Configuration with only ``keep`` active; per-residue settings of others dropped."""
    cfg = copy.deepcopy(config)
    cfg['residues'] = [{'name': label, 'flag': 'on' if label in keep else 'off'} for label in labels]
    init = cfg['init']
    for key in ('initial', 'bounds'):
        if key in init:
            init[key] = {name: value for name, value in init[key].items()
                         if '.' not in name or name.rsplit('.', 1)[0] in keep}
    if 'vary' in init:
        init['vary'] = [name for name in init['vary'] if '.' not in name or name.rsplit('.', 1)[0] in keep]
    for key in ('multistart', 'profile', 'bootstrap', 'profile_interval'):
        init.pop(key, None)
    if 'multistart' in config['init']:
        starts = [{name: value for name, value in start.items()
                   if '.' not in name or name.rsplit('.', 1)[0] in keep}
                  for start in config['init']['multistart'].get('starts', [])]
        init['multistart'] = {**config['init']['multistart'], 'starts': starts}
    return cfg


def information_criteria(chi2, k, n):
    """Gaussian AICc/BIC up to the common constant, with the supplied absolute sigma."""
    aic = chi2 + 2 * k
    aicc = aic + (2 * k * (k + 1) / (n - k - 1) if n - k - 1 > 0 else math.inf)
    return {'aic': aic, 'aicc': aicc, 'bic': chi2 + k * math.log(n)}


def nested_f_test(chi2_small, k_small, chi2_large, k_large, n):
    """F-test of the smaller (shared) model against the larger (individual) model."""
    from scipy.stats import f as f_dist

    extra, dof = k_large - k_small, n - k_large
    if extra <= 0 or dof <= 0 or chi2_large <= 0:
        return {'f': None, 'p_value': None, 'extra_parameters': extra, 'dof': dof,
                'note': 'F-test undefined: models are not nested with spare degrees of freedom'}
    f = max(0., (chi2_small - chi2_large) / extra) / (chi2_large / dof)
    return {'f': float(f), 'p_value': float(f_dist.sf(f, extra, dof)), 'extra_parameters': extra, 'dof': dof,
            'note': 'Nested F-test assuming correct absolute sigma and Gaussian residuals'}


def _load_result(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def run_comparison(config, config_dir, out, *, no_pdf=True, workers=1):
    """Fit the global and the individual models; write comparison JSON/TXT/PDF."""
    out = Path(out).expanduser().absolute()
    if out.exists() or out.is_symlink():
        raise FileExistsError(f'Comparison output already exists: {out}')
    labels = [entry['name'] for entry in config['residues'] if entry.get('flag') == 'on']
    if len(labels) < 2:
        raise ValueError('Model comparison needs at least two active residues')
    out.mkdir(parents=True)
    fits = {}
    global_cfg = _residue_config(config, [e['name'] for e in config['residues']], set(labels))
    global_cfg['Project Name'] = str(out / 'global' / 'fit')
    fits['global'] = _fit(global_cfg, config_dir, no_pdf, workers)
    individual = {}
    for label in labels:
        cfg = _residue_config(config, [e['name'] for e in config['residues']], {label})
        cfg['Project Name'] = str(out / 'individual' / label / 'fit')
        individual[label] = _fit(cfg, config_dir, no_pdf, workers)
    summary = summarize(fits['global'], individual, labels)
    summary['outputs'] = {'global': fits['global']['result_json'],
                          'individual': {label: row['result_json'] for label, row in individual.items()}}
    (out / 'comparison.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    (out / 'comparison.txt').write_text('\n'.join(comparison_lines(summary)) + '\n', encoding='utf-8')
    comparison_pdf(out / 'comparison.pdf', summary)
    return {'comparison_json': str(out / 'comparison.json'), 'comparison_txt': str(out / 'comparison.txt'),
            'comparison_pdf': str(out / 'comparison.pdf')}


def _fit(cfg, config_dir, no_pdf, workers):
    from sbfit import run_config

    row = {'project': cfg['Project Name'], 'result_json': cfg['Project Name'] + '_result.json',
           'success': False, 'message': ''}
    try:
        run_config(cfg, config_dir, no_pdf, workers=workers)
        row['success'] = True
    except (ValueError, RuntimeError, OSError) as exc:
        row['message'] = f'{type(exc).__name__}: {exc}'
    if Path(row['result_json']).exists():
        row['result'] = _load_result(row['result_json'])
        row['success'] = bool(row['result'].get('success')) and row['success']
    return row


def summarize(global_fit, individual, labels):
    summary = {'schema_version': 1, 'residues': labels, 'global': None, 'individual': {},
               'comparison': None, 'warnings': [],
               'interpretation': ('Shared kab/kba (global) versus per-residue exchange (individual, the nesting '
                                  'model). chi2/AICc/BIC/F use the supplied absolute sigma and the configured '
                                  'fixed inputs; they rank descriptions under those assumptions only.')}
    if not global_fit['success']:
        summary['warnings'].append(f'Global fit failed: {global_fit["message"] or global_fit.get("result", {}).get("message", "")}')
    else:
        r = global_fit['result']
        summary['global'] = _entry(r)
    total_chi2, total_k, n_total, complete = 0., 0, 0, True
    for label in labels:
        row = individual[label]
        if not row['success']:
            complete = False
            summary['individual'][label] = {'success': False, 'message': row['message'] or row.get('result', {}).get('message', '')}
            summary['warnings'].append(f'Individual fit for {label} failed: {summary["individual"][label]["message"]}')
            continue
        r = row['result']
        summary['individual'][label] = {'success': True, **_entry(r)}
        total_chi2 += r['chi2']
        total_k += r['n_parameters']
        n_total += r['n_points']
    if summary['global'] is not None and complete:
        g = summary['global']
        if n_total != g['n_points']:
            summary['warnings'].append('Individual fits do not cover the same points as the global fit')
        n = g['n_points']
        summary['comparison'] = {
            'n_points': n,
            'global': {'chi2': g['chi2'], 'k': g['n_parameters'], **information_criteria(g['chi2'], g['n_parameters'], n)},
            'individual_sum': {'chi2': total_chi2, 'k': total_k, **information_criteria(total_chi2, total_k, n)},
            'delta_chi2': g['chi2'] - total_chi2,
            'delta_aicc_individual_minus_global': information_criteria(total_chi2, total_k, n)['aicc'] - information_criteria(g['chi2'], g['n_parameters'], n)['aicc'],
            'delta_bic_individual_minus_global': information_criteria(total_chi2, total_k, n)['bic'] - information_criteria(g['chi2'], g['n_parameters'], n)['bic'],
            'f_test': nested_f_test(g['chi2'], g['n_parameters'], total_chi2, total_k, n),
        }
        summary['comparison']['preferred_by_aicc'] = ('global' if summary['comparison']['delta_aicc_individual_minus_global'] > 0 else 'individual')
        if summary['comparison']['delta_chi2'] < -1e-8 * max(1., total_chi2):
            summary['warnings'].append('The global fit has a lower chi2 than the sum of individual fits; at least one individual fit is not at its minimum.')
    return summary


def _entry(result):
    return {'chi2': result['chi2'], 'n_points': result['n_points'], 'n_parameters': result['n_parameters'],
            'dof': result['dof'], 'kex': result['kex'], 'pB': result['pB'],
            'kex_se': result.get('derived_se', {}).get('kex'), 'pB_se': result.get('derived_se', {}).get('pB'),
            'jacobian_rank': result.get('jacobian_rank'), 'at_bounds': result.get('at_bounds', []),
            'warnings': result.get('warnings', [])}


def comparison_lines(summary):
    def number(value, digits=8):
        return 'unavailable' if value is None else f'{value:.{digits}g}'

    lines = ['Sideband model comparison: shared exchange (global) versus per-residue exchange (individual)',
             summary['interpretation'], '']
    g = summary['global']
    if g is not None:
        lines.append(f'global: chi2 {number(g["chi2"])}, k {g["n_parameters"]}, n {g["n_points"]}, '
                     f'kex {number(g["kex"])} +/- {number(g["kex_se"], 4)}, pB {number(g["pB"])} +/- {number(g["pB_se"], 4)}')
    for label, row in summary['individual'].items():
        if row.get('success'):
            lines.append(f'{label}: chi2 {number(row["chi2"])}, k {row["n_parameters"]}, n {row["n_points"]}, '
                         f'kex {number(row["kex"])} +/- {number(row["kex_se"], 4)}, pB {number(row["pB"])} +/- {number(row["pB_se"], 4)}'
                         + (f'; at bounds: {", ".join(row["at_bounds"])}' if row['at_bounds'] else ''))
        else:
            lines.append(f'{label}: failed; {row.get("message", "")}')
    c = summary['comparison']
    if c is not None:
        lines.extend(['', 'Comparison (individual sum minus global unless noted):',
                      f'delta chi2 (global - individual): {number(c["delta_chi2"])}',
                      f'AICc global {number(c["global"]["aicc"])}, individual {number(c["individual_sum"]["aicc"])}, delta {number(c["delta_aicc_individual_minus_global"])} -> preferred: {c["preferred_by_aicc"]}',
                      f'BIC global {number(c["global"]["bic"])}, individual {number(c["individual_sum"]["bic"])}, delta {number(c["delta_bic_individual_minus_global"])}',
                      f'F-test: F {number(c["f_test"]["f"], 5)}, p {number(c["f_test"]["p_value"], 4)}, extra parameters {c["f_test"]["extra_parameters"]}, dof {c["f_test"]["dof"]}; {c["f_test"]["note"]}'])
    if summary['warnings']:
        lines.extend(['', 'Warnings:', *summary['warnings']])
    return lines


def comparison_pdf(path, summary):
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    with PdfPages(path) as pdf:
        fig, axes = plt.subplots(1, 2, figsize=(11, 5))
        labels = [label for label, row in summary['individual'].items() if row.get('success')]
        for ax, key in zip(axes, ('kex', 'pB')):
            values = [summary['individual'][label][key] for label in labels]
            errors = [summary['individual'][label][f'{key}_se'] or 0. for label in labels]
            ax.errorbar(range(len(labels)), values, yerr=errors, fmt='o', capsize=3, label='individual')
            if summary['global'] is not None:
                g = summary['global']
                ax.axhline(g[key], color='black', linewidth=1, label='global')
                if g[f'{key}_se']:
                    ax.axhspan(g[key] - g[f'{key}_se'], g[key] + g[f'{key}_se'], color='gray', alpha=.2)
            ax.set_xticks(range(len(labels)))
            ax.set_xticklabels(labels)
            ax.set_ylabel(key)
            ax.set_title(f'{key}: per-residue versus shared')
            ax.grid(alpha=.25)
            ax.legend(fontsize=8)
        fig.text(.5, .01, 'Error bars: local standard errors with supplied absolute sigma.', ha='center', fontsize=7)
        fig.tight_layout(rect=(0, .04, 1, 1))
        pdf.savefig(fig)
        plt.close(fig)
        fig = plt.figure(figsize=(8.5, 11))
        fig.text(.05, .95, '\n'.join(comparison_lines(summary)), va='top', family='monospace', fontsize=7.5, linespacing=1.5)
        pdf.savefig(fig)
        plt.close(fig)


THREE_STATE_METHODS = ("Sideband_3st_Linear", "Sideband_3st_Triangle")


def three_state_config(config, method, *, h_ppm_c=None, starts=None):
    """Derive a three-state configuration from a fitted two-state configuration.

    State C gets ``h_ppm_c`` (per residue; default: the state-B proton shift) and
    explicit multistart starts for the extra rates and shift built from the
    configured initial kab/kba; ``starts`` may override the (kbc, kcb, dwC ratio)
    triples. The two-state analyses are not carried over.
    """
    if method not in THREE_STATE_METHODS:
        raise ValueError(f"method must be one of {THREE_STATE_METHODS}")
    cfg = copy.deepcopy(config)
    init = cfg['init']
    init['Method'] = method
    for key in ('multistart', 'profile', 'bootstrap', 'profile_interval', 'kex', 'pB'):
        init.pop(key, None)
    initial = init.setdefault('initial', {})
    kab, kba = float(initial.get('kab', 15.)), float(initial.get('kba', 285.))
    kex = kab + kba
    initial.update(kab=kab, kba=kba, kbc=kex / 3., kcb=kex)
    if method == 'Sideband_3st_Triangle':
        initial.update(kca=0.1 * kex, kac=0.01 * kex)
    labels = [entry['name'] for entry in cfg['residues'] if entry.get('flag') == 'on']
    residues = cfg['sideband']['residues']
    for label in labels:
        shift = residues[label]
        shift['h_ppm_c'] = float((h_ppm_c or {}).get(label, shift['h_ppm_b']))
        dw = float(initial.get(f'{label}.dw_ppm', 1.0))
        initial.setdefault(f'{label}.dwC_ppm', -dw)
        initial.setdefault(f'{label}.R2c', float(initial.get(f'{label}.R2b', 20.)))
    if 'vary' in init:
        extra = ['kbc', 'kcb'] + (['kca', 'kac'] if method == 'Sideband_3st_Triangle' else [])
        init['vary'] = list(init['vary']) + extra + [f'{label}.{key}' for label in labels for key in ('dwC_ppm', 'R2c')]
    for name in ('kbc', 'kcb', 'kca', 'kac'):
        init.get('bounds', {}).pop(name, None)
    triples = starts or [(kex / 3., kex, -1.), (kex, kex / 3., -1.), (kex / 3., kex, 2.), (kex, kex / 3., 2.),
                         (kex, kex, 0.5)]
    explicit = []
    for kbc, kcb, ratio in triples:
        start = {'kbc': float(kbc), 'kcb': float(kcb)}
        for label in labels:
            start[f'{label}.dwC_ppm'] = float(ratio) * float(initial.get(f'{label}.dw_ppm', 1.0))
        explicit.append(start)
    init['multistart'] = {'starts': explicit}
    return cfg


def compare_models(config, config_dir, out, *, models=('Sideband', 'Sideband_3st_Linear'), h_ppm_c=None,
                   no_pdf=True, workers=1):
    """Fit the same data with several Sideband models and compare AICc/BIC."""
    out = Path(out).expanduser().absolute()
    if out.exists() or out.is_symlink():
        raise FileExistsError(f'Comparison output already exists: {out}')
    if len(models) < 2 or len(set(models)) != len(models):
        raise ValueError('compare --models needs at least two distinct models')
    for method in models:
        if method != 'Sideband' and method not in THREE_STATE_METHODS:
            raise ValueError(f'Unknown model {method}')
    if config['init'].get('Method', 'Sideband') != 'Sideband':
        raise ValueError('compare --models starts from a two-state (Sideband) configuration')
    out.mkdir(parents=True)
    fits = {}
    for method in models:
        cfg = copy.deepcopy(config) if method == 'Sideband' else three_state_config(config, method, h_ppm_c=h_ppm_c)
        for key in ('profile', 'bootstrap', 'profile_interval'):
            cfg['init'].pop(key, None)
        cfg['Project Name'] = str(out / method / 'fit')
        (out / method).mkdir()
        (out / method / 'config.json').write_text(json.dumps(cfg, indent=2) + '\n', encoding='utf-8')
        fits[method] = _fit(cfg, config_dir, no_pdf, workers)
    summary = {'schema_version': 1, 'models': list(models), 'fits': {}, 'comparison': None, 'warnings': [],
               'h_ppm_c': h_ppm_c or 'state-B proton shift reused for state C',
               'interpretation': ('The same data fitted with two- and three-state Sideband models. AICc/BIC use '
                                  'the supplied absolute sigma; three-state starts come from a small multistart '
                                  'around the two-state exchange rates. A lower AICc does not prove an extra '
                                  'state: check the populations, the shift of state C, boundary flags and the '
                                  'residual diagnostics.')}
    n = None
    for method, row in fits.items():
        if not row['success']:
            summary['fits'][method] = {'success': False, 'message': row['message'] or row.get('result', {}).get('message', '')}
            summary['warnings'].append(f'{method} fit failed: {summary["fits"][method]["message"]}')
            continue
        r = row['result']
        n = r['n_points']
        entry = {'success': True, 'chi2': r['chi2'], 'n_points': n, 'n_parameters': r['n_parameters'], 'dof': r['dof'],
                 **information_criteria(r['chi2'], r['n_parameters'], n), 'kex': r['kex'], 'pB': r['pB'],
                 'exchange': r.get('exchange'), 'at_bounds': r.get('at_bounds', []), 'warnings': r.get('warnings', []),
                 'reduced_chi2': r['residual_diagnostics']['reduced_chi2'] if 'residual_diagnostics' in r else r['chi2'] / r['dof'],
                 'multistart_success': [a['success'] for a in r.get('multistart', [])] or None,
                 'result_json': row['result_json']}
        summary['fits'][method] = entry
    done = {m: e for m, e in summary['fits'].items() if e.get('success')}
    if len(done) >= 2:
        ranked = sorted(done, key=lambda m: done[m]['aicc'])
        reference = done['Sideband'] if 'Sideband' in done else done[ranked[0]]
        summary['comparison'] = {
            'preferred_by_aicc': ranked[0], 'preferred_by_bic': min(done, key=lambda m: done[m]['bic']),
            'delta_aicc_vs_two_state': {m: e['aicc'] - reference['aicc'] for m, e in done.items()},
            'delta_bic_vs_two_state': {m: e['bic'] - reference['bic'] for m, e in done.items()},
            'delta_chi2_vs_two_state': {m: reference['chi2'] - e['chi2'] for m, e in done.items()},
            'extra_parameters_vs_two_state': {m: e['n_parameters'] - reference['n_parameters'] for m, e in done.items()},
            'note': ('Three-state models reduce to two states only at the boundary of their parameter space, '
                     'so no F-test is reported; AICc/BIC differences and the fitted populations are the evidence.'),
        }
        for m, e in done.items():
            if m != 'Sideband' and e['exchange'] and min(e['exchange']['populations'].values()) < 0.005:
                summary['warnings'].append(f'{m}: a state population is below 0.5%; the extra state is not supported by the data.')
            if m != 'Sideband' and e['at_bounds']:
                absorbing = [name for name in e['at_bounds'] if name in ('kba', 'kcb')]
                summary['warnings'].append(f'{m}: parameters at bounds ({", ".join(e["at_bounds"])}); inspect identifiability.'
                                           + (' A return rate at its lower bound makes a state absorbing, so the '
                                              'reported populations are meaningless and the model is degenerate.' if absorbing else ''))
    (out / 'comparison.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    (out / 'comparison.txt').write_text('\n'.join(models_lines(summary)) + '\n', encoding='utf-8')
    return {'comparison_json': str(out / 'comparison.json'), 'comparison_txt': str(out / 'comparison.txt')}


def models_lines(summary):
    def number(value, digits=8):
        return 'unavailable' if value is None else f'{value:.{digits}g}'

    lines = ['Sideband model comparison: two-state versus three-state exchange', summary['interpretation'],
             f'h_ppm_c: {summary["h_ppm_c"]}', '']
    for method, entry in summary['fits'].items():
        if not entry.get('success'):
            lines.append(f'{method}: failed; {entry.get("message", "")}')
            continue
        populations = entry['exchange']['populations'] if entry.get('exchange') else {}
        lines.append(f'{method}: chi2 {number(entry["chi2"])}, k {entry["n_parameters"]}, n {entry["n_points"]}, '
                     f'reduced chi2 {number(entry["reduced_chi2"], 5)}, AICc {number(entry["aicc"])}, BIC {number(entry["bic"])}; '
                     + 'populations ' + ', '.join(f'{s} {v:.4f}' for s, v in populations.items())
                     + (f'; at bounds: {", ".join(entry["at_bounds"])}' if entry['at_bounds'] else ''))
        if entry.get('exchange'):
            lines.append('  rates: ' + ', '.join(f'{k} {v:.5g}' for k, v in entry['exchange']['rates'].items()))
    c = summary['comparison']
    if c:
        lines.extend(['', f'Preferred by AICc: {c["preferred_by_aicc"]}; by BIC: {c["preferred_by_bic"]}',
                      'delta AICc vs two-state: ' + ', '.join(f'{m} {v:+.3f}' for m, v in c['delta_aicc_vs_two_state'].items()),
                      'delta chi2 (two-state minus model): ' + ', '.join(f'{m} {v:+.4g}' for m, v in c['delta_chi2_vs_two_state'].items()),
                      c['note']])
    if summary['warnings']:
        lines.extend(['', 'Warnings:', *summary['warnings']])
    return lines
