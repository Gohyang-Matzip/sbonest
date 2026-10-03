"""Matched nitrogen-only comparison; preserve all existing study artifacts.

Run: .venv/bin/python manuscript/sideband_30ppm/compare_conventional.py
Primary comparison keeps signed longitudinal deviation magnetization. The
actual ONEST worker, which clamps negative predictions, is a sensitivity check.
"""
import copy
import hashlib
import json
import os
from pathlib import Path
import sys
import time

for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(key, '1')
os.environ.setdefault('MPLBACKEND', 'Agg')
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from scipy.linalg import expm
import matplotlib.pyplot as plt
from estmodel import b1_weights, matrix_calc_worker_chunk
from fit import generate_initial_parameters
from sbfit import SidebandModel

SOURCE = HERE / 'results' / '1200'
OUT = HERE / 'results' / 'conventional_comparison_20261003'
FIG = HERE / 'figures' / 'conventional_comparison_20261003'
TRUTH = np.array([15., 285., 1.08, 120., 3., 1.5, 12., 15.])


def nitrogen_signed(args):
    """ONEST worker's homogeneous 6x6 block, without its final zero clamp.

    Its seventh coordinate is initially zero and remains zero, so the affine
    recovery column is dormant. Dropping that coordinate preserves its exact
    signed propagation; the J=0 control below verifies the observable mapping.
    """
    kab, kba, peak, dw, r1, r2a, r2b, offsets, duration, nu, nuerr, field = args
    p_b = kab / (kab + kba)
    p_a = 1 - p_b
    wx, weights = b1_weights(nu, nuerr)
    wa = (peak - np.asarray(offsets))[:, None] * field * 2 * np.pi
    wb = wa + dw * field * 2 * np.pi
    a = np.zeros((len(wa), len(wx), 6, 6))
    a[..., 0, 0] = a[..., 1, 1] = -kab - r2a
    a[..., 3, 3] = a[..., 4, 4] = -kba - r2b
    a[..., 2, 2], a[..., 5, 5] = -kab-r1, -kba-r1
    a[..., 0, 1], a[..., 1, 0] = -wa, wa
    a[..., 3, 4], a[..., 4, 3] = -wb, wb
    a[..., 1, 2] = a[..., 4, 5] = -wx
    a[..., 2, 1] = a[..., 5, 4] = wx
    a[..., 0, 3] = a[..., 1, 4] = a[..., 2, 5] = kba
    a[..., 3, 0] = a[..., 4, 1] = a[..., 5, 2] = kab
    end = expm(a*duration) @ np.array([0., 0., p_a, 0., 0., p_b])
    return np.sum(weights * end[..., 2] / p_a, axis=1)


class NitrogenModel(SidebandModel):
    """Change only the forward model; inherit data, RF, bounds and optimizer."""
    clipped = False

    def calc(self, parameters, residue, offsets, spectrum):
        nu, nuerr = self.rf_values(parameters, spectrum)
        args = (parameters['kab'], parameters['kba'],
                parameters['dGs'][residue], parameters['dws'][residue],
                parameters['r1s'][residue], parameters['r2as'][residue],
                parameters['r2bs'][residue], np.atleast_1d(offsets),
                spectrum.T, nu, nuerr, spectrum.field)
        forward = matrix_calc_worker_chunk if self.clipped else nitrogen_signed
        values = np.asarray(forward(args))
        return values if np.ndim(offsets) else float(values[0])


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def vector(record, model):
    return np.array([record['parameters'][name]['value'] for name in model.parameter_names])


def metrics(p):
    names = ('kab_s-1', 'kba_s-1', 'v1n_scale', 'peak_ppm', 'dw_ppm',
             'R1N_s-1', 'R2NA_s-1', 'R2NB_s-1')
    values = dict(zip(names, map(float, p)))
    values.update(kex_s_inv=float(p[0]+p[1]), pB=float(p[0]/(p[0]+p[1])))
    return values


def biases(p):
    truth, fitted = metrics(TRUTH), metrics(p)
    return {key: {'absolute': fitted[key]-truth[key],
                  'relative_percent': 100*(fitted[key]/truth[key]-1)} for key in truth}


def prediction(model, p, full_data, dense=None):
    parameters = model.seParam(p)
    return np.array([model.calc(parameters, i, offsets if dense is None else dense, es)
                     for i, es, offsets, _, _ in full_data])


def residual_metrics(values, keep):
    def stats(x):
        x = np.asarray(x).ravel()
        return {'n': int(x.size), 'mean_z': float(x.mean()),
                'rms_z': float(np.sqrt(np.mean(x*x))),
                'max_abs_z': float(np.abs(x).max()), 'chi2': float(x@x)}
    return {'full': stats(values), 'retained': stats(values[:, keep]),
            'masked_sideband_regions': stats(values[:, ~keep]),
            'per_rf': [dict(nominal_v1n_hz=nu, full=stats(row),
                            retained=stats(row[keep]),
                            masked_sideband_regions=stats(row[~keep]))
                       for nu, row in zip((25, 50, 100), values)]}


def run_fit(config, clipped, dataset_label, exact_mask=None):
    model = NitrogenModel(config, SOURCE)
    model.clipped = clipped
    if exact_mask is not None:
        for _, es, *_ in model._prepare_data():
            for attr in ('offset', 'int', 'intstd'):
                setattr(es, attr, np.asarray(getattr(es, attr))[exact_mask].tolist())
    model._initializing = True
    try:
        standard = model._expand_initial(generate_initial_parameters(model, config['init']))
    finally:
        model._initializing = False
    for name, val in config['init'].get('initial', {}).items():
        standard[model.parameter_names.index(name)] = val
    fit_config = copy.deepcopy(config['init'])
    fit_config['initial'] = {}
    starts = {'standard_grid': standard, 'truth': TRUTH,
              'perturbed': np.array([24., 376., .95, 120.08, 2.7, 1.2, 18., 25.])}
    attempts, successes = [], []
    for label, initial in starts.items():
        began = time.monotonic()
        try:
            p, cov = model.fit(initial, fitting_config=fit_config)
            record = model.diagnostics(p, cov)
            record.pop('config')
            record['method'] = 'N-only legacy clipped' if clipped else 'N-only signed'
            record.update(start=label, initial=initial.tolist(), success=True,
                          parameter_vector=p.tolist(), metrics=metrics(p), biases=biases(p),
                          reduced_chi2=record['chi2']/record['dof'],
                          elapsed_s=time.monotonic()-began)
            attempts.append(record)
            successes.append((record['chi2'], p.copy(), label))
            print(dataset_label, record['method'], label,
                  f"chi2/dof={record['reduced_chi2']:.6g}",
                  f"kex={record['kex']:.6g} pB={record['pB']:.8g}", flush=True)
        except Exception as exc:
            attempts.append(dict(start=label, initial=initial.tolist(), success=False,
                                 error=f'{type(exc).__name__}: {exc}',
                                 elapsed_s=time.monotonic()-began))
    if not successes:
        raise RuntimeError(f'No successful fit for {dataset_label}')
    _, best, best_label = min(successes, key=lambda x: x[0])
    result = {'best_start': best_label, 'attempts': attempts,
              'best': next(row for row in attempts if row['start'] == best_label)}
    return model, best, result


def main():
    OUT.mkdir(exist_ok=False)
    FIG.mkdir(exist_ok=False)
    config = json.loads((SOURCE/'full.json').read_text())
    original = json.loads((HERE/'results/summary.json').read_text())
    nh = SidebandModel(config, SOURCE)
    full = nh._prepare_data()
    ppm, observed, sigma = (np.array([row[i] for row in full]) for i in (2, 3, 4))
    assert all(np.array_equal(row, ppm[0]) for row in ppm)
    keep = ~((abs((ppm[0]-120)*full[0][1].field) >= 1250) &
             (abs((ppm[0]-120)*full[0][1].field) <= 1850))
    assert keep.sum() == 101 and len(keep) == 147
    for index, nu in enumerate((25, 50, 100)):
        masked = np.loadtxt(SOURCE/f'masked_{nu}.txt', skiprows=5)
        np.testing.assert_array_equal(masked, np.c_[ppm[index,keep], observed[index,keep], sigma[index,keep]])
    zero_config = copy.deepcopy(config)
    zero_config['sideband']['decoupling']['J_hz'] = 0
    zero = SidebandModel(zero_config, SOURCE)
    nitrogen = NitrogenModel(config, SOURCE)
    controls = []
    for control_p in (TRUTH, np.array([24., 376., .95, 120.08, 2.7, 1.2, 18., 25.])):
        signed = prediction(nitrogen, control_p, full)
        nh_zero = prediction(zero, control_p, full)
        nitrogen.clipped = True
        clipped = prediction(nitrogen, control_p, full)
        nitrogen.clipped = False
        np.testing.assert_allclose(signed, nh_zero, atol=1e-10, rtol=0)
        np.testing.assert_allclose(clipped, np.maximum(signed, 0), atol=1e-10, rtol=0)
        controls.append({'parameters': control_p.tolist(),
                         'signed_J0_max_abs_difference': float(abs(signed-nh_zero).max()),
                         'legacy_vs_clipped_signed_max_abs_difference': float(abs(clipped-np.maximum(signed,0)).max()),
                         'negative_signed_points': int((signed < 0).sum()),
                         'zero_clamp_max_change': float(abs(clipped-signed).max())})
    report = {'scope': 'Synthetic 1.2 GHz; one noise realization and noiseless profiles; eight free parameters.',
              'primary_model': 'Signed nitrogen-only two-state Bloch-McConnell propagation; dormant identity coordinate removed from ONEST worker and legacy zero clamp disabled.',
              'legacy_check': 'Actual estmodel.matrix_calc_worker_chunk, with its final zero clamp; observations never clipped.',
              'optimizer': 'Same SidebandModel.fit optimizer, bounds and RF-scale handling for both nitrogen-only variants.',
              'bounds': {'kab': [0,None], 'kba': [1e-8,None], 'v1n_scale': [.8,1.2],
                         'peak_ppm': [None,None], 'dw_ppm': [None,None],
                         'R1N': [0,None], 'R2NA': [0,None], 'R2NB': [0,None]},
              'truth': metrics(TRUTH), 'J0_controls': controls, 'fits': {},
              'interpretation': 'Noisy deviations from truth are realization-specific errors. Noiseless misspecified-fit shifts estimate the pseudo-true displacement for this grid; neither establishes experimental bias. Standard errors under a misspecified model are not validated confidence intervals.'}
    arrays = {'ppm': ppm[0], 'observed': observed, 'sigma': sigma, 'keep': keep,
              'dense_ppm': np.linspace(105,135,2401)}
    for kind in ('full','masked','exact'):
        existing = copy.deepcopy(original['fields']['1200']['fits'][kind])
        p = vector(existing, nh)
        key = 'NH_'+kind
        predicted = prediction(nh,p,full)
        values = (predicted-observed)/sigma
        existing.update(metrics=metrics(p), biases=biases(p), reduced_chi2=existing['chi2']/existing['dof'],
                        residuals_against_noisy_full=residual_metrics(values,keep),
                        provenance='Existing summary.json result; not refitted in this comparison.')
        existing.pop('config')
        report['fits'][key] = existing
        arrays[key+'_prediction'] = predicted
        if kind == 'full':
            arrays[key+'_dense'] = prediction(nh,p,full,arrays['dense_ppm'])
    for clipped in (False, True):
        for noiseless, masked in ((False,False),(False,True),(True,False),(True,True)):
            cfg = copy.deepcopy(config)
            cfg['datasets'] = [f"{'exact' if noiseless else 'masked' if masked else 'full'}_{nu}.txt" for nu in (25,50,100)]
            # Exact masked input is an in-memory subset; original data files stay intact.
            label = ('exact_' if noiseless else 'noisy_') + ('masked' if masked else 'full')
            model, p, record = run_fit(cfg, clipped, label,
                                       exact_mask=keep if noiseless and masked else None)
            key = ('legacy_' if clipped else 'N_')+label
            predicted = prediction(model,p,full)
            record['best']['residuals_against_noisy_full'] = residual_metrics((predicted-observed)/sigma,keep)
            report['fits'][key] = record
            arrays[key+'_prediction'] = predicted
            if not clipped and not noiseless and not masked:
                arrays[key+'_dense'] = prediction(model,p,full,arrays['dense_ppm'])
            write_json(OUT/'summary.json',report)
    for key in ('noisy_full','noisy_masked','exact_full','exact_masked'):
        signed, clipped = (report['fits'][prefix+key]['best'] for prefix in ('N_','legacy_'))
        report.setdefault('clipping_sensitivity',{})[key] = {
            'chi2_signed': signed['chi2'], 'chi2_legacy': clipped['chi2'],
            'parameter_metric_difference_legacy_minus_signed': {name: clipped['metrics'][name]-signed['metrics'][name] for name in signed['metrics']}}
    source_files = [Path(__file__),ROOT/'estmodel.py',ROOT/'sbfit.py',ROOT/'sideband.py',
                    HERE/'results/summary.json', *sorted(SOURCE.glob('*.txt')),SOURCE/'full.json']
    report['source_sha256'] = {str(path.relative_to(ROOT)):hashlib.sha256(path.read_bytes()).hexdigest() for path in source_files}
    write_json(OUT/'summary.json',report)
    np.savez(OUT/'predictions.npz',**arrays)
    make_figure(arrays,report)
    print('PASS: signed J=0 parity, matched source masks, multistart fits, legacy clipping sensitivity',flush=True)


def make_figure(a, report):
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,
                         'axes.spines.top':False,'axes.spines.right':False})
    fig, axes = plt.subplots(2,3,figsize=(11.0,5.7),sharex=True,sharey='row',
                             layout='constrained',gridspec_kw={'height_ratios':[1.35,1]})
    colors = {'NH_full':'#0072B2','N_noisy_full':'#D55E00'}
    names = {'NH_full':'Coupled NH fit','N_noisy_full':'Nitrogen-only fit'}
    nfield = 1200*.101329118
    for col,nu in enumerate((25,50,100)):
        upper,lower=axes[:,col]
        for ax in (upper,lower):
            for lo,hi in ((-1850,-1250),(1250,1850)):
                ax.axvspan(max(105,120+lo/nfield),min(135,120+hi/nfield),color='.92',zorder=0)
            ax.set_xlim(105,135)
            ax.set_xticks([105,115,125,135])
        upper.plot(a['ppm'],a['observed'][col],'.',color='.35',ms=3,label='Synthetic observations',zorder=3)
        for key in colors:
            upper.plot(a['dense_ppm'],a[key+'_dense'][col],color=colors[key],lw=1.2,label=names[key])
            z=(a[key+'_prediction'][col]-a['observed'][col])/a['sigma'][col]
            lower.plot(a['ppm'],z,'o-',color=colors[key],ms=2.2,lw=.65)
        upper.set_title(f"{'ABC'[col]}   Nominal ν₁N = {nu} Hz",loc='left',fontweight='bold')
        lower.set_title(f"{'DEF'[col]}   Standardized residuals",loc='left',fontweight='bold')
        lower.axhline(0,color='.25',lw=.6)
        lower.set_xlabel('¹⁵N irradiation position (ppm)')
    axes[0,0].set_ylabel('I / I₀')
    axes[1,0].set_ylabel('(Prediction − observation) / σ')
    axes[0,0].set_ylim(-.02,.59)
    handles, labels = axes[0,0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, fontsize=9,
               loc='outside lower center', ncol=3)
    fig.suptitle('1.2 GHz · Same 441 observations · Eight fitted parameters per model',fontsize=11)
    for suffix in ('png','pdf'):
        fig.savefig(FIG/f'figure5_conventional_comparison.{suffix}',dpi=300,bbox_inches='tight')
    plt.close(fig)
    caption = ('Figure 5. Nitrogen-only and coupled NH fits to the same synthetic 1.2 GHz data. '
               'A–C, all 147 offsets at each nominal nitrogen RF amplitude are fitted jointly, with eight '
               'free parameters including a shared RF scale. Points are the original noisy observations; '
               'curves are forward-model predictions from the full-data fits. D–F, standardized residuals '
               '(prediction minus observation, divided by the supplied absolute intensity error σ=0.001). '
               'Shading marks the prespecified |offset from 120 ppm|=1250–1850 Hz mask, clipped to the '
               'acquisition window; no shaded points were excluded from these plotted fits. The primary '
               'nitrogen-only model uses the conventional two-state Bloch–McConnell matrix with its legacy '
               'zero clamp disabled to match signed longitudinal-deviation detection. J=0 parity and an '
               'actual clipped-ONEST sensitivity comparison are provided in the source JSON. These are '
               'simulations from one parameter set and one noise realization, not experimental validation.')
    write_json(FIG/'figure5_conventional_comparison.json',{
        'claim':'A nitrogen-only model cannot reproduce the coupled-spin features in this synthetic complete profile.',
        'panel_plan':{'A-C':'Same observations and both full-data fits at 25/50/100 Hz.',
                      'D-F':'Standardized residuals on the original acquired grid.'},
        'caption':caption,'generator':'manuscript/sideband_30ppm/compare_conventional.py',
        'sources':['results/conventional_comparison_20261003/summary.json',
                   'results/conventional_comparison_20261003/predictions.npz'],
        'evidence':'Synthetic data only; identical observations, error weights and free-parameter count.'})


if __name__ == '__main__':
    main()
