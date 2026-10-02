"""Three-peak extension: independent R1H/R2H for each peak, shared exchange/RF.

Run: .venv/bin/python check_peakwise_proton_relaxation.py --out results/peakwise_H_fit
Research check only; production run.py configuration remains unchanged.
"""

import argparse
import copy
import hashlib
import json
from pathlib import Path

from check_proton_relaxation import SOURCE, jacobian, residual
import numpy as np
from scipy.optimize import least_squares
from fit import _block_jacobian
from run import load_config
from sbfit import SidebandModel
from sideband import composite_segments, profile
from test_sideband import reference

PEAKS = [
    ("A1", [120.0, 3.0, 1.5, 12.0, 15.0, 2.0, 25.0], [6.2, 6.5]),
    ("G2", [118.5, 2.2, 1.3, 10.0, 20.0, 1.2, 18.0], [7.2, 7.5]),
    ("S3", [121.5, -2.5, 1.7, 16.0, 25.0, 3.0, 35.0], [9.8, 9.6]),
]
LOCAL = ["peak_ppm", "dw_ppm", "R1", "R2a", "R2b", "R1H", "R2H"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    out = parser.parse_args().out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    rng = np.random.default_rng(20261003)
    cfg = load_config(SOURCE / "full.json")
    cfg["sideband"]["residues"] = {
        name: dict(zip(("h_ppm_a", "h_ppm_b"), h)) for name, _, h in PEAKS
    }
    cfg["residues"] = [{"name": name, "flag": "on"} for name, _, _ in PEAKS]
    field, sigma = 121.5949416, 0.001
    ppm = np.linspace(105.0, 135.0, 147)
    paths = {kind: [] for kind in ("exact", "noisy")}
    seg = composite_segments(70e-6, "RR")
    for rf in (25.0, 50.0, 100.0):
        lines = {
            kind: [f"{field}\n0.4\n{rf} 0\n# offset intensity error\n"]
            for kind in paths
        }
        for name, local, h in PEAKS:
            shift, dw, r1, r2a, r2b, r1h, r2h = local
            args = dict(
                T=0.4,
                nu=rf * 1.08,
                kab=15.0,
                kba=285.0,
                dw=dw * field,
                r1=r1,
                r2a=r2a,
                r2b=r2b,
                ha=(h[0] - 8.5) * 1200,
                hb=(h[1] - 8.5) * 1200,
                r1h=r1h,
                r2h=r2h,
            )
            y = reference(seg, (ppm - shift) * field, **args)
            np.testing.assert_allclose(
                y, profile(seg, (ppm - shift) * field, **args), atol=2e-10
            )
            observed = y + rng.normal(0, sigma, len(y))
            for kind in paths:
                lines[kind].append(f"# {name} R2a: {r2a} R2b: {r2b} dw: {dw}\n")
                vals = y if kind == "exact" else observed
                lines[kind].extend(
                    f"{x:.14g} {v:.14g} {sigma}\n" for x, v in zip(ppm, vals)
                )
        for kind in paths:
            path = out / f"{kind}_{rf:g}.txt"
            path.write_text("".join(lines[kind]))
            paths[kind].append(str(path))
    names = ["kab", "kba", "v1n_scale"] + [
        f"{peak}.{key}" for peak, _, _ in PEAKS for key in LOCAL
    ]
    truth = np.r_[15.0, 285.0, 1.08, np.array([local for _, local, _ in PEAKS]).ravel()]
    lower = np.r_[
        0.0, 1e-8, 0.8, np.tile([-np.inf, -np.inf, 0.0, 0.0, 0.0, 0.0, 0.0], 3)
    ]
    upper = np.r_[np.inf, np.inf, 1.2, np.tile([np.inf] * 5 + [200.0, 500.0], 3)]
    h_indices = [3 + 7 * i + k for i in range(3) for k in (5, 6)]
    fits = {}
    for kind, label, vary_h in (
        ("exact", "exact_peakwise", True),
        ("noisy", "fixed_H", False),
        ("noisy", "fit_R2H_peakwise", "R2H"),
        ("noisy", "fit_both_peakwise", True),
    ):
        models = []
        for peak, _, _ in PEAKS:
            c = copy.deepcopy(cfg)
            c["datasets"] = paths[kind]
            c["residues"] = [
                {"name": name, "flag": "on" if name == peak else "off"}
                for name, _, _ in PEAKS
            ]
            model = SidebandModel(c)
            model._fit_data = model._prepare_data()
            models.append(model)

        def full_fun(p):
            return np.concatenate(
                [
                    residual(model, np.r_[p[:3], p[3 + 7 * i : 10 + 7 * i]])
                    for i, model in enumerate(models)
                ]
            )

        np.testing.assert_allclose(
            full_fun(truth), 0, atol=1e-6
        ) if kind == "exact" else None
        grouped = _block_jacobian(full_fun, [441] * 3, 3, 7, relative_step=1e-5)
        free = np.array(
            [
                i
                for i in range(len(truth))
                if i not in h_indices
                or vary_h is True
                or (vary_h == "R2H" and names[i].endswith(".R2H"))
            ]
        )
        start = truth.copy()
        start[:3] = [18.0, 310.0, 1.03]
        for i in range(3):
            start[3 + 7 * i] += 0.03
            start[4 + 7 * i] *= 0.97
            for k in (5, 6):
                if 3 + 7 * i + k in free:
                    start[3 + 7 * i + k] *= 1.5

        def expand(q):
            p = start.copy()
            p[free] = q
            return p

        result = least_squares(
            lambda q: full_fun(expand(q)),
            start[free],
            jac=lambda q: grouped(expand(q))[:, free],
            bounds=(lower[free], upper[free]),
            x_scale="jac",
            ftol=1e-9,
            xtol=1e-9,
            gtol=1e-9,
            max_nfev=200,
        )
        p = expand(result.x)
        fun = lambda q: full_fun(expand(q))
        central = jacobian(fun, result.x, lower[free], upper[free])
        grouped_final = grouped(p)[:, free]
        relative = np.linalg.norm(central - grouped_final, axis=0) / np.linalg.norm(
            central, axis=0
        )
        assert relative.max() < 0.02, relative
        norms = np.linalg.norm(central, axis=0)
        sv = np.linalg.svd(central / norms, compute_uv=False)
        rank = int(np.count_nonzero(sv > sv[0] * 1e-8))
        assert rank == len(free)
        _, s, vt = np.linalg.svd(central, full_matrices=False)
        cov = np.zeros((len(p), len(p)))
        cov[np.ix_(free, free)] = (vt.T / s**2) @ vt
        se = np.sqrt(np.diag(cov))
        chi2 = float(result.fun @ result.fun)
        correlations = {
            peak: float(cov[8 + 7 * i, 9 + 7 * i] / se[8 + 7 * i] / se[9 + 7 * i])
            for i, (peak, _, _) in enumerate(PEAKS)
            if vary_h is True
        }
        data = {
            "success": bool(result.success),
            "message": result.message,
            "chi2": chi2,
            "n_parameters": len(free),
            "n_points": 1323,
            "dof": 1323 - len(free),
            "kex": float(p[0] + p[1]),
            "pB": float(p[0] / (p[0] + p[1])),
            "v1n_scale": float(p[2]),
            "parameters": {
                n: {
                    "value": float(p[i]),
                    "stderr": float(se[i]),
                    "vary": bool(i in free),
                }
                for i, n in enumerate(names)
            },
            "jacobian_rank": rank,
            "scaled_condition": float(sv[0] / sv[-1]),
            "at_bounds": [names[free[i]] for i in np.flatnonzero(result.active_mask)],
            "peakwise_R1H_R2H_correlation": correlations,
            "grouped_vs_central_jacobian_max_relative_column_change": float(
                relative.max()
            ),
            "nfev": result.nfev,
        }
        (out / f"{label}.json").write_text(
            json.dumps(data, indent=2, allow_nan=False) + "\n"
        )
        fits[label] = data
        print(
            label,
            "chi2",
            round(chi2, 5),
            "free",
            len(free),
            "H",
            np.round(p[h_indices], 4),
            "H_se",
            np.round(se[h_indices], 4),
            flush=True,
        )
        assert result.success, result.message
        if kind == "exact":
            np.testing.assert_allclose(p, truth, atol=0.002, rtol=1e-5)
            assert chi2 < 1e-8
    (out / "summary.json").write_text(
        json.dumps(
            {
                "seed": 20261003,
                "scope": "H rates independent per peak; within each peak, shared across A/B states and three RF datasets",
                "truth": dict(zip(names, truth.tolist())),
                "config": cfg,
                "input_files": paths,
                "script_sha256": hashlib.sha256(
                    Path(__file__).read_bytes()
                ).hexdigest(),
                "fits": fits,
            },
            indent=2,
        )
        + "\n"
    )
    print(
        "PASS: independent generator, per-peak 24-parameter exact recovery, grouped/central Jacobian parity",
        flush=True,
    )


if __name__ == "__main__":
    main()
