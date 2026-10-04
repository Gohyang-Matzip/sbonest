"""Compare two nitrogen RF amplitudes against the existing three-RF experiment.

Run: .venv/bin/python scripts/check_two_rf.py --out results/two_RF_check --workers 3
Uses exactly the existing offsets, observations, noise, and 64 random H-rate sets.
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]  # repository root: flat modules live there
sys.path.insert(0, str(ROOT))

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import copy
import csv
import hashlib
import json
import time

from check_random_proton_rates import ROOT, SOURCE, NAMES, PEAKS, LOWER, UPPER
from check_proton_relaxation import jacobian
import numpy as np
from scipy.optimize import least_squares
from fit import _block_jacobian
from sbfit import SidebandModel

PRIOR = ROOT / "results/random_H_rates"
PAIRS = {"25_50": [0, 1], "25_100": [0, 2], "50_100": [1, 2]}
CACHE = {}


def models_for(indices, kind):
    key = (tuple(indices), kind)
    if key not in CACHE:
        source = json.loads(SOURCE.read_text())
        models = []
        for peak in PEAKS:
            cfg = copy.deepcopy(source["config"])
            cfg["datasets"] = [source["input_files"][kind][i] for i in indices]
            cfg["residues"] = [
                {"name": p, "flag": "on" if p == peak else "off"} for p in PEAKS
            ]
            m = SidebandModel(cfg)
            m._fit_data = m._prepare_data()
            assert sum(len(row[2]) for row in m._fit_data) == 147 * len(indices)
            models.append(m)
        CACHE[key] = models
    return CACHE[key]


def fit_job(job):
    began = time.monotonic()
    models = models_for(job["indices"], job["kind"])
    for m, h in zip(models, job["H_rates"]):
        for d in m.decoupling:
            d["R1H"], d["R2H"] = h

    def fun(p):
        return np.concatenate(
            [
                m.errFunc(np.r_[p[:3], p[3 + 5 * i : 8 + 5 * i]])
                for i, m in enumerate(models)
            ]
        )

    sizes = [sum(len(row[2]) for row in m._fit_data) for m in models]
    grouped = _block_jacobian(fun, sizes, 3, 5, relative_step=1e-6)
    result = least_squares(
        fun,
        job["start"],
        jac=grouped,
        bounds=(LOWER, UPPER),
        x_scale="jac",
        ftol=1e-9,
        xtol=1e-9,
        gtol=1e-9,
        max_nfev=300,
    )
    p = result.x
    out = {k: v for k, v in job.items() if k != "start"}
    out.update(
        success=bool(result.success),
        message=result.message,
        nfev=result.nfev,
        chi2=float(result.fun @ result.fun),
        n_points=len(result.fun),
        n_parameters=len(p),
        dof=len(result.fun) - len(p),
        parameters=dict(zip(NAMES, map(float, p))),
        kex=float(p[0] + p[1]),
        pB=float(p[0] / (p[0] + p[1])),
        nominal_RF_Hz=models[0].dataset.v1s,
        v1n_hz=(np.array(models[0].dataset.v1s) * p[2]).tolist(),
        at_bounds=[NAMES[i] for i in np.flatnonzero(result.active_mask)],
        elapsed_s=time.monotonic() - began,
    )
    if job.get("covariance"):
        j = jacobian(fun, p, LOWER, UPPER)
        norms = np.linalg.norm(j, axis=0)
        singular = np.linalg.svd(j / norms, compute_uv=False)
        rank = int(np.count_nonzero(singular > singular[0] * 1e-8))
        assert rank == len(p), (job["id"], rank)
        _, s, vt = np.linalg.svd(j, full_matrices=False)
        cov = (vt.T / s**2) @ vt
        out["stderr"] = dict(zip(NAMES, map(float, np.sqrt(np.diag(cov)))))
        out["stderr"]["kex"] = float(np.sqrt(cov[0, 0] + cov[1, 1] + 2 * cov[0, 1]))
        g = np.array([p[1], -p[0]]) / (p[0] + p[1]) ** 2
        out["stderr"]["pB"] = float(np.sqrt(g @ cov[:2, :2] @ g))
        out["jacobian_rank"] = rank
        out["scaled_condition"] = float(singular[0] / singular[-1])
        err = np.linalg.norm(j - grouped(p), axis=0) / norms
        out["jacobian_check_max_relative_column_change"] = float(err.max())
        assert err.max() < 0.005, (job["id"], err)
    return out


def save(out, r):
    (out / (r["id"] + ".json")).write_text(
        json.dumps(r, indent=2, allow_nan=False) + "\n"
    )
    print(
        r["id"],
        "success",
        r["success"],
        "chi2",
        round(r["chi2"], 3),
        "kex",
        round(r["kex"], 4),
        "RF",
        round(r["parameters"]["v1n_scale"], 6),
        "s",
        round(r["elapsed_s"], 1),
        flush=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=3)
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("workers must be positive")
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    meta = json.loads((PRIOR / "metadata.json").read_text())
    prior = json.loads((PRIOR / "results.json").read_text())
    prior_by_id = {r["id"]: r for r in prior}
    # Require the preceding experiment's data and source to remain unchanged.
    for name, sha in meta["source_sha256"].items():
        assert hashlib.sha256((ROOT / name).read_bytes()).hexdigest() == sha, name
    truth = meta["truth"]
    ptrue = np.array([truth[n] for n in NAMES])
    htrue = [[truth[f"{p}.R1H"], truth[f"{p}.R2H"]] for p in PEAKS]
    start = ptrue.copy()
    start[:3] = [18.0, 310.0, 1.03]
    results = []
    control = fit_job(
        dict(
            id="control_three_RF",
            pair="25_50_100",
            indices=[0, 1, 2],
            kind="noisy",
            arm="control",
            draw=-1,
            H_rates=htrue,
            start=start.tolist(),
            covariance=True,
        )
    )
    save(out, control)
    results.append(control)
    assert control["success"]
    assert abs(control["chi2"] - prior_by_id["baseline_noisy"]["chi2"]) < 1e-6
    np.testing.assert_allclose(
        [control["parameters"][n] for n in NAMES],
        [prior_by_id["baseline_noisy"]["parameters"][n] for n in NAMES],
        atol=1e-5,
        rtol=1e-7,
    )
    baselines = {}
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        jobs = [
            dict(
                id=f"{pair}_baseline_{kind}",
                pair=pair,
                indices=indices,
                kind=kind,
                arm="baseline",
                draw=-1,
                H_rates=htrue,
                start=start.tolist(),
                covariance=True,
            )
            for pair, indices in PAIRS.items()
            for kind in ("noisy", "exact")
        ]
        for future in as_completed([pool.submit(fit_job, j) for j in jobs]):
            r = future.result()
            save(out, r)
            results.append(r)
            assert r["success"]
            assert r["n_points"] == 882 and r["dof"] == 864
            if r["kind"] == "exact":
                np.testing.assert_allclose(
                    [r["parameters"][n] for n in NAMES], ptrue, atol=1e-5
                )
                assert r["chi2"] < 1e-8
            else:
                baselines[r["pair"]] = r
        jobs = []
        # Interleave pairs so all designs progress together on the same H draws.
        for original in meta["jobs"]:
            for pair, indices in PAIRS.items():
                j = copy.deepcopy(original)
                j.update(
                    id=f"{pair}_{original['id']}",
                    pair=pair,
                    indices=indices,
                    source_fit=original["id"],
                    start=[baselines[pair]["parameters"][n] for n in NAMES],
                )
                jobs.append(j)
        source_files = [
            Path(__file__),
            ROOT / "scripts/check_random_proton_rates.py",
            PRIOR / "metadata.json",
            PRIOR / "results.json",
        ]
        metadata = dict(
            interpretation="Two nitrogen RF amplitudes at one 1.2 GHz B0, not two B0 fields",
            pairs=PAIRS,
            nominal_RF_Hz=[25.0, 50.0, 100.0],
            truth=truth,
            random_H_source="results/random_H_rates/metadata.json",
            draws_per_arm=16,
            arms=meta["arms"],
            workers=args.workers,
            jobs=jobs,
            acquisition="Identical observations/noise on retained files; 882 rather than 1323 points",
            equal_time="Assuming noise proportional to 1/sqrt(scans), 1.5x scans on two RF profiles scales SE by sqrt(2/3). Model bias is unchanged under uniform weighting; no new acquisition simulated.",
            source_sha256={
                str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                for p in source_files
            },
            upstream_source_sha256=meta["source_sha256"],
        )
        (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
        for future in as_completed([pool.submit(fit_job, j) for j in jobs]):
            r = future.result()
            save(out, r)
            results.append(r)
    failures = [r["id"] for r in results if not r["success"]]
    if failures:
        raise RuntimeError(f"Fits failed; outputs preserved: {failures}")
    # Recheck the largest shifts and the largest residual in each pair.
    targets = ["kex", "v1n_scale"] + [
        f"{p}.{k}" for p in PEAKS for k in ("dw_ppm", "R1", "R2a", "R2b")
    ]

    def value(r, n):
        return r[n] if n == "kex" else r["parameters"][n]

    for pair, indices in PAIRS.items():
        rows = [r for r in results if r["pair"] == pair and r["arm"] in meta["arms"]]
        b = baselines[pair]

        def impact(r):
            return max(
                abs(value(r, n) - value(b, n)) / abs(300.0 if n == "kex" else truth[n])
                for n in targets
            )

        chosen = {
            r["id"]: r
            for r in (max(rows, key=impact), max(rows, key=lambda r: r["chi2"]))
        }
        for original in chosen.values():
            for kind, arm in (("exact", "exact_check"), ("noisy", "restart_check")):
                p0 = np.array([original["parameters"][n] for n in NAMES])
                if kind == "noisy":
                    p0[:3] *= [1.2, 0.85, 0.98]
                    for i in range(3):
                        p0[3 + 5 * i] += 0.02
                        p0[4 + 5 * i] *= 1.03
                        p0[5 + 5 * i : 8 + 5 * i] *= 1.1
                p0 = np.clip(p0, LOWER, UPPER)
                r = fit_job(
                    dict(
                        id=f"{arm}_{original['id']}",
                        pair=pair,
                        indices=indices,
                        kind=kind,
                        arm=arm,
                        draw=original["draw"],
                        H_rates=original["H_rates"],
                        source_fit=original["id"],
                        start=p0.tolist(),
                    )
                )
                save(out, r)
                results.append(r)
                assert r["success"]
                if kind == "noisy":
                    assert abs(r["chi2"] - original["chi2"]) < 0.002, (
                        r["id"],
                        r["chi2"],
                        original["chi2"],
                    )
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
                "pair",
                "arm",
                "kind",
                "draw",
                "success",
                "n_points",
                "chi2",
                "dof",
                "kex",
                "pB",
            )
        }
        row.update(r["parameters"])
        for p, h in zip(PEAKS, r["H_rates"]):
            row.update({p + ".R1H": h[0], p + ".R2H": h[1]})
        rows.append(row)
    with (out / "results.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(
        "PASS:",
        len(results),
        "fits; three-RF parity; exact recovery; all paired H draws; restart checks",
        flush=True,
    )


if __name__ == "__main__":
    main()
