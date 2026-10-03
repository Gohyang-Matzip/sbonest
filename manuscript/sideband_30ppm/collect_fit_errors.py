"""Recover full covariance at the saved optima without changing any fits.

Run with the repository .venv Python. Retain absolute sigma=0.001 weights.
"""
from copy import deepcopy
import hashlib
import json
import os
from pathlib import Path
import sys

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(name, '1')
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
import numpy as np
from fit import _block_jacobian
from sbfit import SidebandModel
from compare_conventional import NitrogenModel

OUT = HERE/'results/fit_errors_20261003'


def main():
    OUT.mkdir(exist_ok=False)
    original = json.loads((HERE/'results/summary.json').read_text())
    comparison = json.loads((HERE/'results/conventional_comparison_20261003/summary.json').read_text())
    carbon = json.loads((HERE/'results/carbon_configuration_20261003/summary.json').read_text())
    sources = [HERE/'results/summary.json', HERE/'results/conventional_comparison_20261003/summary.json',
               HERE/'results/carbon_configuration_20261003/summary.json']
    records = []
    for kind, record in original['fields']['1200']['fits'].items():
        records.append(('NH_'+kind, SidebandModel, record, record['config'], HERE/'results/1200', False))
    for kind in ('N_noisy_full', 'N_noisy_masked', 'N_exact_full', 'N_exact_masked'):
        record = comparison['fits'][kind]['best']
        config = deepcopy(original['fields']['1200']['fits']['full']['config'])
        dataset = 'exact' if 'exact' in kind else 'masked' if 'masked' in kind else 'full'
        config['datasets'] = [f'{dataset}_{nu}.txt' for nu in (25, 50, 100)]
        records.append((kind, NitrogenModel, record, config, HERE/'results/1200', kind=='N_exact_masked'))
    records.append(('CH_exact', SidebandModel, carbon, carbon['config'], HERE/'results/carbon_configuration_20261003', False))
    report = {'method': 'Unscaled local 1 SE from weighted Jacobian at each saved optimum; sigma=0.001.',
              'derived_method': 'sqrt(g.T @ covariance @ g), retaining cross-covariances.',
              'interpretation': 'Formal within-model SE. Misspecified-model SE excludes model discrepancy; noiseless SE uses assigned hypothetical sigma, not replicate variability.',
              'fits': {}}
    for key, cls, record, cfg, folder, mask in records:
        model = cls(cfg, folder)
        model._fit_data = model._prepare_data()
        if mask:
            model._fit_data = [(i, es, x[keep], y[keep], e[keep])
                for i, es, x, y, e in model._fit_data
                for keep in [~((abs((x-120)*es.field)>=1250)&(abs((x-120)*es.field)<=1850))]]
        names = model.parameter_names
        p = np.array([record['parameters'][n]['value'] for n in names])
        size = sum(len(row[2]) for row in model._fit_data)
        jac = _block_jacobian(model.errFunc, [size], 2+len(model.rf_names), len(model.local_keys), relative_step=1e-6)(p)
        _, singular, vt = np.linalg.svd(jac, full_matrices=False)
        covariance = (vt.T/singular**2)@vt
        se = np.sqrt(np.diag(covariance))
        old_se = np.array([record['parameters'][n]['stderr'] for n in names])
        np.testing.assert_allclose(se, old_se, rtol=2e-6, atol=1e-10)
        residual = model.errFunc(p)
        np.testing.assert_allclose(residual@residual, record['chi2'], rtol=1e-6, atol=1e-10)
        rank = np.linalg.matrix_rank(jac/np.linalg.norm(jac, axis=0), tol=1e-8)
        assert rank==len(p) and size==record['n_points']
        values = deepcopy(record['parameters'])
        kab, kba = p[:2]
        gradients = np.zeros((3,len(p)))
        gradients[0,:2] = 1
        gradients[1,:2] = 100*np.array([kba,-kab])/(kab+kba)**2
        ia, idw = names.index('A1.peak_ppm'), names.index('A1.dw_ppm')
        gradients[2,[ia,idw]] = 1
        derived = np.array([kab+kba,100*kab/(kab+kba),p[ia]+p[idw]])
        derived_se = np.sqrt(np.diag(gradients@covariance@gradients.T))
        for name,v,e in zip(('kex','pB_percent','delta_B'),derived,derived_se):
            values[name] = dict(value=float(v),stderr=float(e),vary=False,derived=True)
        if 'derived_se' in record:
            np.testing.assert_allclose(derived_se[:2], [record['derived_se']['kex'],100*record['derived_se']['pB']], rtol=2e-6)
        report['fits'][key] = dict(parameters=values, parameter_order=names,
            covariance=covariance.tolist(), n_points=size, dof=record['dof'], chi2=record['chi2'],
            rank=int(rank), max_relative_stderr_difference=float(np.max(abs(se/old_se-1))))
        sources.extend(folder/name for name in cfg['datasets'])
        print(key, 'SE parity passed; kex, pB%, deltaB SE:', derived_se, flush=True)
    sources += [Path(__file__), HERE/'compare_conventional.py', ROOT/'fit.py', ROOT/'sbfit.py', ROOT/'sideband.py']
    report['source_sha256'] = {str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    (OUT/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print('PASS: all saved SE and objectives reproduced; full covariance propagated; no fits changed.')


if __name__=='__main__':
    main()
