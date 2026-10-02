"""Test proton-relaxation identifiability with the existing 30 ppm NH model.

Run: .venv/bin/python check_proton_relaxation.py --out results/proton_relaxation_fit
Uses existing synthetic data; does not add parameters to the production CLI.
R1H/R2H are shared across states and RF files for the single residue A1.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import time

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[key] = "1"
os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np
import scipy
from scipy.optimize import least_squares

from run import load_config
from sbfit import SidebandModel
from test_sideband import reference

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "manuscript/sideband_30ppm/results/1200"
TRUTH = np.array([15.0, 285.0, 1.08, 120.0, 3.0, 1.5, 12.0, 15.0, 2.0, 25.0])
LOWER = np.array([0.0, 1e-8, 0.8, -np.inf, -np.inf, 0.0, 0.0, 0.0, 0.0, 0.0])
UPPER = np.array(
    [np.inf, np.inf, 1.2, np.inf, np.inf, np.inf, np.inf, np.inf, 200.0, 500.0]
)


def model_for(kind):
    model = SidebandModel(load_config(SOURCE / f"{kind}.json"))
    model._fit_data = model._prepare_data()
    return model


def residual(model, parameters):
    for decoupling in model.decoupling:
        decoupling["R1H"], decoupling["R2H"] = parameters[-2:]
    return model.errFunc(parameters[:-2])


def jacobian(fun, values, lower, upper, step_scale=1e-5):
    """Dense central differences; a finite one-sided step at a physical bound."""
    columns = []
    for i, value in enumerate(values):
        step = step_scale * max(abs(value), 1.0)
        plus, minus = values.copy(), values.copy()
        plus[i] = min(value + step, upper[i])
        minus[i] = max(value - step, lower[i])
        columns.append((fun(plus) - fun(minus)) / (plus[i] - minus[i]))
    return np.column_stack(columns)


def fit_case(model, start, free, label, out):
    began = time.monotonic()
    free = np.asarray(free, dtype=int)
    names = model.parameter_names + ["R1H", "R2H"]

    def expand(values):
        p = np.array(start, dtype=float)
        p[free] = values
        return p

    def fun(values):
        r = residual(model, expand(values))
        assert np.isfinite(r).all()
        return r

    def jac(values):
        return jacobian(fun, values, LOWER[free], UPPER[free])

    result = least_squares(
        fun,
        start[free],
        jac=jac,
        bounds=(LOWER[free], UPPER[free]),
        x_scale="jac",
        ftol=1e-9,
        xtol=1e-9,
        gtol=1e-9,
        max_nfev=200,
    )
    p = expand(result.x)
    chi2 = float(result.fun @ result.fun)
    j = jac(result.x)
    norms = np.linalg.norm(j, axis=0)
    singular = np.linalg.svd(j / np.where(norms > 0, norms, 1.0), compute_uv=False)
    rank = int(np.count_nonzero(singular > singular[0] * 1e-8))
    covariance = np.zeros((len(p), len(p)))
    if rank == len(free):
        _, s, vt = np.linalg.svd(j, full_matrices=False)
        covariance[np.ix_(free, free)] = (vt.T / s**2) @ vt
    else:
        covariance[np.ix_(free, free)] = np.nan
    se = np.sqrt(np.diag(covariance))
    corr = []
    for i in free:
        for k in free:
            if k > i and (i >= 8 or k >= 8) and se[i] > 0 and se[k] > 0:
                corr.append(
                    {
                        "parameters": [names[i], names[k]],
                        "correlation": float(covariance[i, k] / se[i] / se[k]),
                    }
                )
    corr.sort(key=lambda pair: -abs(pair["correlation"]))
    kex = p[0] + p[1]
    gradient = np.zeros((3, len(p)))
    gradient[0, :2] = 1
    gradient[1, :2] = [p[1] / kex**2, -p[0] / kex**2]
    gradient[2, 2] = 1
    derived_se = np.sqrt(np.diag(gradient @ covariance @ gradient.T))
    # A tighter difference step checks that local covariance is not a step artifact.
    j_finer = jacobian(fun, result.x, LOWER[free], UPPER[free], step_scale=5e-6)
    relative_column_change = np.linalg.norm(j_finer - j, axis=0) / np.maximum(
        norms, 1e-30
    )
    n = len(result.fun)
    k = len(free)
    data = {
        "label": label,
        "success": bool(result.success),
        "message": result.message,
        "parameters": {
            name: {
                "value": float(p[i]),
                "stderr": float(se[i]) if np.isfinite(se[i]) else None,
                "vary": bool(i in free),
            }
            for i, name in enumerate(names)
        },
        "kex": float(kex),
        "pB": float(p[0] / kex),
        "v1n_scale": float(p[2]),
        "derived_stderr": {
            name: float(v) if np.isfinite(v) else None
            for name, v in zip(("kex", "pB", "v1n_scale"), derived_se)
        },
        "chi2": chi2,
        "dof": n - k,
        "reduced_chi2": chi2 / (n - k),
        "n_points": n,
        "n_parameters": k,
        "jacobian_rank": rank,
        "scaled_condition": float(singular[0] / singular[-1]),
        "at_bounds": [names[free[i]] for i in np.flatnonzero(result.active_mask)],
        "proton_correlations": corr,
        "jacobian_step_check_max_relative_column_change": float(
            relative_column_change.max()
        ),
        "covariance_note": "Local absolute-sigma covariance; not rescaled by chi2/dof. Bound-limited errors are not confidence intervals.",
        "AICc_known_sigma_up_to_common_constant": chi2
        + 2 * k
        + 2 * k * (k + 1) / (n - k - 1),
        "nfev": result.nfev,
        "elapsed_s": time.monotonic() - began,
    }
    (out / f"{label}.json").write_text(
        json.dumps(data, indent=2, allow_nan=False) + "\n"
    )
    print(
        label,
        "success",
        result.success,
        "chi2",
        round(chi2, 6),
        "H",
        np.round(p[-2:], 5),
        "H_se",
        np.round(se[-2:], 5),
        "seconds",
        round(data["elapsed_s"], 1),
        flush=True,
    )
    assert result.success, (label, result.message)
    assert relative_column_change.max() < 0.02, (label, relative_column_change)
    return p, data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    noisy = model_for("full")
    exact = model_for("exact")
    names = noisy.parameter_names + ["R1H", "R2H"]
    assert len(names) == 10
    # The loaded independent-generator data must match the reused propagator.
    assert np.max(np.abs(residual(exact, TRUTH))) < 1e-6
    # Check changed proton rates against the independent complex density matrix.
    test = TRUTH.copy()
    test[-2:] = [7.0, 40.0]
    residual(noisy, test)
    es = noisy.dataset.res[0].estSpecs[0]
    ppm = np.array([106.5, 119.8, 123.0, 132.5])
    expected = reference(
        noisy.segments[0],
        (ppm - 120) * es.field,
        T=es.T,
        nu=es.v1 * 1.08,
        kab=15.0,
        kba=285.0,
        dw=3 * es.field,
        r1=1.5,
        r2a=12.0,
        r2b=15.0,
        ha=-2.3 * 1200,
        hb=-2 * 1200,
        r1h=7.0,
        r2h=40.0,
    )
    actual = noisy.calc(noisy.seParam(test[:-2]), 0, ppm, es)
    np.testing.assert_allclose(actual, expected, atol=2e-10, rtol=2e-10)
    inputs = [SOURCE / name for name in ("full.json", "exact.json")]
    inputs += [Path(p) for p in noisy.config["datasets"] + exact.config["datasets"]]
    inputs += [
        ROOT / name
        for name in (
            "sbfit.py",
            "sideband.py",
            "test_sideband.py",
            "check_proton_relaxation.py",
        )
    ]
    metadata = {
        "truth": dict(zip(names, TRUTH.tolist())),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "source_sha256": {
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in inputs
        },
        "conditions": {
            "field_H_MHz": 1200,
            "field_N_MHz": es.field,
            "T_s": 0.4,
            "p90_s": 70e-6,
            "cycle": "RR",
            "window_ppm": [105, 135],
            "nominal_RF_Hz": [25, 50, 100],
            "n_points": 441,
            "sigma": 0.001,
            "proton_relaxation_scope": "One shared R1H/R2H across both states and all RF datasets for A1",
        },
        "proton_bounds_s_inverse": {"R1H": [0, 200], "R2H": [0, 500]},
        "config": noisy.config,
        "checks": [
            "Independent synthetic generator parity",
            "Changed H rates match independent density-matrix propagator",
        ],
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    start = np.array([18.0, 320.0, 1.0, 120.05, 2.7, 1.2, 13.0, 18.0, 8.0, 40.0])
    fits = {}
    p_exact, fits["exact_both"] = fit_case(exact, start, range(10), "exact_both", out)
    np.testing.assert_allclose(p_exact, TRUTH, atol=2e-3, rtol=1e-5)
    assert fits["exact_both"]["chi2"] < 1e-8
    for label, free in (
        ("fixed_H", list(range(8))),
        ("fit_R1H", list(range(8)) + [8]),
        ("fit_R2H", list(range(8)) + [9]),
        ("fit_both", list(range(10))),
    ):
        p0 = start.copy()
        for i in (8, 9):
            if i not in free:
                p0[i] = TRUTH[i]
        p, fits[label] = fit_case(noisy, p0, free, label, out)
        if label == "fixed_H":
            baseline = json.loads((SOURCE / "full_result.json").read_text())
            assert abs(fits[label]["chi2"] - baseline["chi2"]) < 1e-4
        if label == "fit_both":
            best = p
    # Alternate starts test whether apparent degeneracy is only initialization.
    for i, h in enumerate(([0.2, 5.0], [60.0, 100.0])):
        p0 = start.copy()
        p0[-2:] = h
        p, fits[f"both_start_{i}"] = fit_case(
            noisy, p0, range(10), f"both_start_{i}", out
        )
        if fits[f"both_start_{i}"]["chi2"] < fits["fit_both"]["chi2"]:
            best = p
    p0 = best.copy()
    p0[2] = TRUTH[2]
    _, fits["both_calibrated_v1n"] = fit_case(
        noisy, p0, [i for i in range(10) if i != 2], "both_calibrated_v1n", out
    )
    (out / "fits.json").write_text(json.dumps(fits, indent=2, allow_nan=False) + "\n")
    # Profile each H rate: reoptimize all remaining N/RF/H parameters at each value.
    profiles = {}
    for index, grid in (
        (8, [0.0, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 80.0, 160.0]),
        (9, [5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 60.0]),
    ):
        rows = []
        free = [i for i in range(10) if i != index]
        for value in grid:
            p0 = best.copy()
            p0[index] = value
            _, data = fit_case(
                noisy, p0, free, f"profile_{names[index]}_{value:g}", out
            )
            rows.append(data)
        profiles[names[index]] = rows
    (out / "profiles.json").write_text(
        json.dumps(profiles, indent=2, allow_nan=False) + "\n"
    )
    print(
        "PASS: parity, exact 10-parameter recovery, original fit reproduction, convergence, finite-difference stability",
        flush=True,
    )


if __name__ == "__main__":
    main()
