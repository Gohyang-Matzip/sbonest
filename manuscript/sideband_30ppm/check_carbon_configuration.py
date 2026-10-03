"""Check a 13C-1H configuration with the existing, unchanged SBONEST engine.

Run from the repository with .venv/bin/python. Synthetic pair-model check only.
"""
import hashlib
import json
import os
from pathlib import Path
import sys

for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(key, '1')
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from sbfit import SidebandModel
from sideband import composite_segments, profile
from test_sideband import reference


def main():
    out = HERE / 'results/carbon_configuration_20261003'
    out.mkdir(exist_ok=False)
    cfg = json.loads((HERE / 'results/1200/full.json').read_text())
    cfg['Project Name'] = 'carbon_pair_configuration_check'
    cfg['init']['initial'] = {}
    cfg['init']['max_nfev'] = 120
    dec = cfg['sideband']['decoupling']
    dec.update(h_carrier_ppm=4.7, J_hz=140.0)
    cfg['sideband']['residues']['A1'] = dict(h_ppm_a=3.8, h_ppm_b=4.1)
    carbon_mhz = 301.8  # Explicit illustrative 13C frequency at nominal 1.2 GHz.
    truth = np.array([15., 285., 1.08, 55., 1.2, 1.2, 15., 20.])
    initial = np.array([18., 330., 1.0, 55.02, 1.1, 1.3, 13., 23.])
    ppm = np.linspace(45, 65, 121)
    segments = composite_segments(70e-6, 'RR')
    args = dict(T=.4, kab=15., kba=285., dw=1.2*carbon_mhz, r1=1.2,
                r2a=15., r2b=20., ha=(3.8-4.7)*1200, hb=(4.1-4.7)*1200,
                J=140., r1h=2., r2h=25.)
    cfg['datasets'] = []
    differences = []
    for rf in (25., 50., 100.):
        y = reference(segments, (ppm-55)*carbon_mhz, nu=rf*1.08, **args)
        predicted = profile(segments, (ppm-55)*carbon_mhz, nu=rf*1.08, **args)
        np.testing.assert_allclose(predicted, y, atol=2e-10, rtol=0)
        differences.append(float(abs(predicted-y).max()))
        name = f'carbon_{rf:g}.txt'
        with (out/name).open('x') as f:
            f.write(f'{carbon_mhz}\n0.4\n{rf} 0\n# offset intensity error\n'
                    '# A1 R2a: 13 R2b: 23 dw: 1.1\n')
            np.savetxt(f, np.c_[ppm, y, np.full(len(ppm), .001)], fmt='%.15g')
        cfg['datasets'].append(name)
    (out/'config.json').write_text(json.dumps(cfg, indent=2)+'\n')
    model = SidebandModel(cfg, out)
    for _, es, offsets, intensity, _ in model._prepare_data():
        np.testing.assert_allclose(model.calc(model.seParam(truth), 0, offsets, es),
                                   intensity, atol=2e-10, rtol=0)
    p, cov = model.fit(initial)
    report = model.diagnostics(p, cov)
    np.testing.assert_allclose(p, truth, atol=1e-5, rtol=1e-6)
    assert model.rank == 8 and report['chi2']/report['dof'] < 1e-10
    report.update(scope='Noiseless synthetic isolated 13C-1H pair; no experimental or multi-spin validation.',
                  carbon_larmor_mhz=carbon_mhz, truth=truth.tolist(), initial=initial.tolist(),
                  fit=p.tolist(), max_reference_difference=max(differences),
                  relative_parameter_errors_percent=(100*(p/truth-1)).tolist(),
                  source_sha256={str(f.relative_to(ROOT)):hashlib.sha256(f.read_bytes()).hexdigest()
                                 for f in (Path(__file__), ROOT/'sideband.py', ROOT/'sbfit.py', ROOT/'test_sideband.py')})
    (out/'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k:report[k] for k in ('scope','max_reference_difference','fit','chi2','dof')}, indent=2))


if __name__ == '__main__':
    main()
