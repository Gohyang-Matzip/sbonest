"""Refit N/exchange/RF parameters with random, peak-specific fixed H rates.

Reuses the independent-generator three-peak data from check_peakwise_proton_relaxation.py.
Run: .venv/bin/python check_random_proton_rates.py --out results/random_H_rates
"""

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import copy
import csv
import hashlib
import json
from pathlib import Path
import time

from check_proton_relaxation import ROOT, jacobian
import numpy as np
from scipy.optimize import least_squares
from fit import _block_jacobian
from sbfit import SidebandModel

SOURCE = ROOT / "results/peakwise_H_fit/summary.json"
PEAKS = ("A1", "G2", "S3")
LOCAL = ("peak_ppm", "dw_ppm", "R1", "R2a", "R2b")
NAMES = ["kab", "kba", "v1n_scale"] + [f"{p}.{k}" for p in PEAKS for k in LOCAL]
LOWER = np.r_[0.0, 1e-8, 0.8, np.tile([-np.inf, -np.inf, 0.0, 0.0, 0.0], 3)]
UPPER = np.r_[np.inf, np.inf, 1.2, np.tile([np.inf] * 5, 3)]
MODELS = {}


def initialize():
    source = json.loads(SOURCE.read_text())
    for kind in ("noisy", "exact"):
        MODELS[kind] = []
        for peak in PEAKS:
            cfg = copy.deepcopy(source["config"])
            cfg["datasets"] = source["input_files"][kind]
            cfg["residues"] = [
                {"name": p, "flag": "on" if p == peak else "off"} for p in PEAKS
            ]
            model = SidebandModel(cfg)
            model._fit_data = model._prepare_data()
            MODELS[kind].append(model)


def fit_job(job):
    began = time.monotonic()
    models = MODELS[job["kind"]]
    h = np.asarray(job["H_rates"])
    for model, rates in zip(models, h):
        for d in model.decoupling:
            d["R1H"], d["R2H"] = rates

    def fun(p):
        return np.concatenate(
            [
                m.errFunc(np.r_[p[:3], p[3 + 5 * i : 8 + 5 * i]])
                for i, m in enumerate(models)
            ]
        )

    grouped = _block_jacobian(fun, [441] * 3, 3, 5, relative_step=1e-6)
    result = least_squares(
        fun,
        job["start"],
        jac=grouped,
        bounds=(LOWER, UPPER),
        x_scale="jac",
        ftol=1e-9,
        xtol=1e-9,
        gtol=1e-9,
        max_nfev=200,
    )
    p = result.x
    data = {k: v for k, v in job.items() if k != "start"}
    data.update(
        success=bool(result.success),
        message=result.message,
        nfev=result.nfev,
        chi2=float(result.fun @ result.fun),
        dof=len(result.fun) - len(p),
        parameters=dict(zip(NAMES, map(float, p))),
        kex=float(p[0] + p[1]),
        pB=float(p[0] / (p[0] + p[1])),
        v1n_hz=(np.array([25.0, 50.0, 100.0]) * p[2]).tolist(),
        at_bounds=[NAMES[i] for i in np.flatnonzero(result.active_mask)],
        elapsed_s=time.monotonic() - began,
    )
    if job.get("covariance"):
        j = jacobian(fun, p, LOWER, UPPER)
        norms = np.linalg.norm(j, axis=0)
        scaled = np.linalg.svd(j / norms, compute_uv=False)
        rank = int(np.count_nonzero(scaled > scaled[0] * 1e-8))
        assert rank == len(p)
        _, s, vt = np.linalg.svd(j, full_matrices=False)
        cov = (vt.T / s**2) @ vt
        data["stderr"] = dict(zip(NAMES, map(float, np.sqrt(np.diag(cov)))))
        data["stderr"]["kex"] = float(np.sqrt(cov[0, 0] + cov[1, 1] + 2 * cov[0, 1]))
        g = np.array([p[1], -p[0]]) / (p[0] + p[1]) ** 2
        data["stderr"]["pB"] = float(np.sqrt(g @ cov[:2, :2] @ g))
        data["jacobian_rank"] = rank
        error = np.linalg.norm(j - grouped(p), axis=0) / norms
        data["jacobian_check_max_relative_column_change"] = float(error.max())
        assert error.max() < 0.005
    return data


def write_result(out, result):
    (out / (result["id"] + ".json")).write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )
    print(
        result["id"],
        "success",
        result["success"],
        "chi2",
        round(result["chi2"], 3),
        "kex",
        round(result["kex"], 4),
        "RF",
        round(result["parameters"]["v1n_scale"], 6),
        "s",
        round(result["elapsed_s"], 1),
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--draws", type=int, default=16)
    parser.add_argument("--workers", type=int, default=3)
    args = parser.parse_args()
    if args.draws < 1 or args.workers < 1:
        parser.error("draws and workers must be positive")
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    source = json.loads(SOURCE.read_text())
    truth = source["truth"]
    h_true = np.array([[truth[f"{p}.R1H"], truth[f"{p}.R2H"]] for p in PEAKS])
    true_p = np.array([truth[n] for n in NAMES])
    initialize()
    start = true_p.copy()
    start[:3] = [18.0, 310.0, 1.03]
    results = []
    for kind in ("noisy", "exact"):
        job = dict(
            id=f"baseline_{kind}",
            arm="baseline",
            draw=-1,
            kind=kind,
            H_rates=h_true.tolist(),
            start=start.tolist(),
            covariance=True,
        )
        result = fit_job(job)
        write_result(out, result)
        assert result["success"]
        results.append(result)
        if kind == "noisy":
            baseline = result
            assert abs(result["chi2"] - source["fits"]["fixed_H"]["chi2"]) < 0.001
        else:
            np.testing.assert_allclose(
                [result["parameters"][n] for n in NAMES], true_p, atol=1e-5
            )
            assert result["chi2"] < 1e-8
    seed = 20261004
    rng = np.random.default_rng(seed)
    near = rng.uniform(0.8, 1.2, (args.draws, 3, 2))
    broad = np.exp(rng.uniform(np.log(0.5), np.log(2.0), (args.draws, 3, 2)))
    jobs = []
    for arm in ("near_both", "broad_R1H", "broad_R2H", "broad_both"):
        for draw in range(args.draws):
            factors = (near if arm == "near_both" else broad)[draw].copy()
            if arm == "broad_R1H":
                factors[:, 1] = 1
            if arm == "broad_R2H":
                factors[:, 0] = 1
            jobs.append(
                dict(
                    id=f"{arm}_{draw:02d}",
                    arm=arm,
                    draw=draw,
                    kind="noisy",
                    H_rates=(h_true * factors).tolist(),
                    factors=factors.tolist(),
                    start=[baseline["parameters"][n] for n in NAMES],
                )
            )
    source_files = [
        SOURCE,
        Path(__file__),
        ROOT / "check_proton_relaxation.py",
        ROOT / "sbfit.py",
        ROOT / "sideband.py",
        ROOT / "fit.py",
    ]
    source_files += [Path(p) for paths in source["input_files"].values() for p in paths]
    metadata = dict(
        seed=seed,
        draws_per_arm=args.draws,
        workers=args.workers,
        truth=truth,
        arms={
            "near_both": "Independent uniform factors [0.8,1.2] for each peak and H rate",
            "broad_R1H": "Independent log-uniform R1H factors [0.5,2]; R2H true",
            "broad_R2H": "Independent log-uniform R2H factors [0.5,2]; R1H true",
            "broad_both": "Same paired broad factors, both H rates randomized",
        },
        meaning="Random H rates are fixed during refitting, not merely random initial guesses. The observed noise realization is unchanged.",
        scope="Peak-specific H rates, shared within a peak across RF files and A/B states. Shared kab/kba/v1n_scale. 18 free N/exchange/RF parameters.",
        source_sha256={
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in source_files
        },
        jobs=jobs,
    )
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    with ProcessPoolExecutor(max_workers=args.workers, initializer=initialize) as pool:
        for future in as_completed([pool.submit(fit_job, job) for job in jobs]):
            result = future.result()
            write_result(out, result)
            results.append(result)
    failed = [r["id"] for r in results if not r["success"]]
    if failed:
        raise RuntimeError(f"Fits failed; retain outputs and investigate: {failed}")
    # Select the largest parameter shifts and largest residual penalty for checks.
    targets = ["kex", "v1n_scale"] + [
        f"{p}.{k}" for p in PEAKS for k in ("dw_ppm", "R1", "R2a", "R2b")
    ]

    def value(r, name):
        return r[name] if name == "kex" else r["parameters"][name]

    def impact(r):
        return max(
            abs(value(r, n) - value(baseline, n))
            / abs(300.0 if n == "kex" else truth[n])
            for n in targets
        )

    randomized = [r for r in results if r["arm"] != "baseline"]
    selected = sorted(randomized, key=impact, reverse=True)[:2]
    worst_chi2 = max(randomized, key=lambda r: r["chi2"])
    if worst_chi2["id"] not in [r["id"] for r in selected]:
        selected.append(worst_chi2)
    for original in selected:
        # Noise-free refits isolate model-assumption bias from the fixed noise draw.
        p0 = [original["parameters"][n] for n in NAMES]
        for kind, label in (("exact", "exact_check"), ("noisy", "restart_check")):
            alternate = np.array(p0)
            if kind == "noisy":
                alternate[:3] *= [1.2, 0.85, 0.98]
                for i in range(3):
                    alternate[3 + 5 * i] += 0.02
                    alternate[4 + 5 * i] *= 1.03
                    alternate[5 + 5 * i : 8 + 5 * i] *= 1.1
            job = dict(
                id=f"{label}_{original['id']}",
                arm=label,
                draw=original["draw"],
                kind=kind,
                H_rates=original["H_rates"],
                start=alternate.tolist(),
                source_fit=original["id"],
            )
            result = fit_job(job)
            write_result(out, result)
            assert result["success"]
            if label == "restart_check":
                assert abs(result["chi2"] - original["chi2"]) < 0.002
            results.append(result)
    results.sort(key=lambda r: r["id"])
    (out / "results.json").write_text(
        json.dumps(results, indent=2, allow_nan=False) + "\n"
    )
    rows = []
    for r in results:
        row = {
            k: r[k]
            for k in (
                "id",
                "arm",
                "kind",
                "draw",
                "success",
                "chi2",
                "dof",
                "kex",
                "pB",
            )
        }
        row.update(r["parameters"])
        for peak, h in zip(PEAKS, r["H_rates"]):
            row.update({peak + ".R1H": h[0], peak + ".R2H": h[1]})
        rows.append(row)
    with (out / "results.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(
        "PASS:",
        len(results),
        "fits; unchanged-noise controls, exact recovery, central Jacobian, selected exact/restart checks",
        flush=True,
    )


if __name__ == "__main__":
    main()
