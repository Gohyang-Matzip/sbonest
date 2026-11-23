"""Repeat the 600/800 MHz SBONEST versus ONEST Matrix comparison over many noise seeds.

Run: .venv/bin/python scripts/compare_field_seeds.py \
       --source results/field_comparison_600_800_20261003_02/inputs --out results/field_seeds_01 --seeds 50
The archived study (scripts/compare_field_models.py) has one noise draw per field. This
study adds independent Gaussian noise (sigma 0.001) to the same noise-free profiles for
every seed and fits both models with the same three starts, so bias and the calibration
of the local standard errors can be separated from one draw's luck. Synthetic data only.
"""

# ruff: noqa: E402 -- compare_field_models sets BLAS limits and import paths first.
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import argparse
import contextlib
import copy
import csv
import hashlib
import io
import json
import shutil
import tempfile

import compare_field_models as base
import numpy as np

SIGMA = .001
BASE_SEED = 20261123
MODELS = ("SBONEST", "ONEST_Matrix")
QUANTITIES = ("kex", "pB", "v1n_scale")


def noisy_inputs(source, folder, field, seed):
    """Write the three RF datasets of one field with fresh noise; return their configuration."""
    config = json.loads((source / "full.json").read_text())
    rng = np.random.default_rng(np.random.SeedSequence([BASE_SEED, field, seed]))
    names = []
    for rf in (25, 50, 100):
        lines = (source / f"exact_{rf}.txt").read_text().splitlines()
        out = lines[:5]
        for line in lines[5:]:
            offset, value, _ = line.split()
            out.append(f"{offset} {float(value) + SIGMA * rng.normal():.15g} {SIGMA}")
        names.append(f"noisy_{rf}.txt")
        (folder / names[-1]).write_text("\n".join(out) + "\n")
    config["datasets"] = names
    return config


def one_seed(source, field, seed):
    """Best-of-three-starts fits of both models to one noise draw."""
    folder = Path(tempfile.mkdtemp(prefix=f"field-seeds-{field}-{seed}-"))
    try:
        config = noisy_inputs(source, folder, field, seed)
        row = {"field_mhz": field, "seed": seed}
        for label in MODELS:
            with contextlib.redirect_stdout(io.StringIO()):
                _, p, result = base.fit_starts(copy.deepcopy(config), folder, label)
            best = result["best"]
            values = {"kex": best["kex"], "pB": best["pB"], "v1n_scale": float(p[2])}
            errors = {**best["derived_se"], "v1n_scale": best["parameters"]["v1n_scale"]["stderr"]}
            row[label] = {"reduced_chi2": best["reduced_chi2"], "multistart_agreement": result["multistart_agreement"],
                          **{name: {"value": float(values[name]), "se": float(errors[name])} for name in QUANTITIES}}
        return row
    finally:
        shutil.rmtree(folder)


def summarize(rows, truth):
    """Bias, scatter and standard-error calibration of every model and quantity, plus paired differences."""
    out = {}
    count = len(rows)
    for label in MODELS:
        out[label] = {"mean_reduced_chi2": float(np.mean([row[label]["reduced_chi2"] for row in rows])),
                      "starts_agree": int(sum(row[label]["multistart_agreement"] for row in rows))}
        for name in QUANTITIES:
            values = np.array([row[label][name]["value"] for row in rows])
            errors = np.array([row[label][name]["se"] for row in rows])
            sd = float(values.std(ddof=1))
            bias = float(values.mean() - truth[name])
            out[label][name] = {
                "mean": float(values.mean()), "bias": bias, "bias_se": sd / np.sqrt(count),
                "bias_z": bias / (sd / np.sqrt(count)), "empirical_sd": sd, "mean_local_se": float(errors.mean()),
                "coverage_1se": float(np.mean(np.abs(values - truth[name]) <= errors)),
                "coverage_2se": float(np.mean(np.abs(values - truth[name]) <= 2 * errors))}
    out["paired_ONEST_minus_SBONEST"] = {}
    for name in QUANTITIES:
        delta = np.array([row["ONEST_Matrix"][name]["value"] - row["SBONEST"][name]["value"] for row in rows])
        out["paired_ONEST_minus_SBONEST"][name] = {"mean": float(delta.mean()),
                                                   "se": float(delta.std(ddof=1) / np.sqrt(count))}
    return out


def report_lines(report):
    """Plain-text table of the per-field summaries."""
    lines = ["600/800 MHz SBONEST versus ONEST Matrix over independent noise seeds (synthetic data)",
             report["scope"], report["uncertainty"], ""]
    truth = report["truth"]
    for field, entry in report["fields"].items():
        lines.append(f"{field} MHz, {report['seeds']} seeds")
        lines.append("model         quantity    truth      bias +/- SE(bias)       z   empirical SD  mean local SE  "
                     "within 1 SE  within 2 SE")
        for label in MODELS:
            for name in QUANTITIES:
                item = entry[label][name]
                lines.append(f"{label:<13} {name:<10} {truth[name]:<8.6g} {item['bias']:+.5g} +/- {item['bias_se']:.3g}"
                             f"   {item['bias_z']:+6.2f}   {item['empirical_sd']:.5g}   {item['mean_local_se']:.5g}"
                             f"   {item['coverage_1se']:.2f}   {item['coverage_2se']:.2f}")
            lines.append(f"{label:<13} mean chi2/dof {entry[label]['mean_reduced_chi2']:.5g}; three starts agree in "
                         f"{entry[label]['starts_agree']} of {report['seeds']} seeds")
        for name, item in entry["paired_ONEST_minus_SBONEST"].items():
            lines.append(f"paired ONEST - SBONEST {name}: {item['mean']:+.5g} +/- {item['se']:.3g}")
        lines.append("")
    return lines


def main():
    """Command line of the multi-seed field comparison."""
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", type=Path, required=True, help="Folder holding 600/ and 800/ with full.json and exact_*.txt")
    parser.add_argument("--out", type=Path, required=True, help="New output directory")
    parser.add_argument("--seeds", type=int, default=50)
    parser.add_argument("--workers", type=int, default=1, help="Processes; seeds are independent")
    args = parser.parse_args()
    if args.seeds < 2 or args.workers < 1:
        parser.error("--seeds must be at least 2 and --workers positive")
    args.out.mkdir(parents=True, exist_ok=False)
    truth = base.metrics(base.TRUTH)
    truth = {"kex": truth["kex_s_inv"], "pB": truth["pB"], "v1n_scale": truth["v1n_scale"]}
    report = {"scope": f"Noise-free 30 ppm profiles of the archived benchmark plus Gaussian noise, sigma {SIGMA}; "
                       f"draw k of a field uses SeedSequence([{BASE_SEED}, field_MHz, k]).",
              "uncertainty": "Bias SE is the empirical SD over seeds divided by sqrt(seeds). Local SE is the "
                             "absolute-sigma covariance error of each fit; for a calibrated error about 68% of fits "
                             "lie within 1 SE of the truth and 95% within 2. Fixed R1H=2, R2H=25 s^-1 as in the "
                             "archived study. ONEST_Matrix is the RF-matched forward model with its zero clamp.",
              "truth": truth, "seeds": args.seeds, "sigma": SIGMA, "fields": {}}
    tasks = [(args.source / str(field), field, seed) for field in (600, 800) for seed in range(args.seeds)]
    if args.workers > 1:
        import multiprocessing

        with multiprocessing.get_context("spawn").Pool(args.workers) as pool:
            rows = pool.starmap(one_seed, tasks)
    else:
        rows = [one_seed(*task) for task in tasks]
    for field in (600, 800):
        report["fields"][str(field)] = summarize([row for row in rows if row["field_mhz"] == field], truth)
    sources = [Path(__file__), Path(base.__file__), base.STUDY / "compare_conventional.py",
               *(base.ROOT / name for name in ("sbfit.py", "sideband.py", "estmodel.py", "est_data.py", "fit.py")),
               *(args.source / str(field) / name for field in (600, 800)
                 for name in ("full.json", "exact_25.txt", "exact_50.txt", "exact_100.txt"))]
    report["source_sha256"] = {str(path.resolve().relative_to(base.ROOT)) if path.resolve().is_relative_to(base.ROOT)
                               else str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in sources}
    base.save_json(args.out / "summary.json", report)
    with (args.out / "seeds.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["field_MHz", "seed", "model", "reduced_chi2", "starts_agree",
                         *(f"{name}{suffix}" for name in QUANTITIES for suffix in ("", "_se"))])
        for row in rows:
            for label in MODELS:
                item = row[label]
                writer.writerow([row["field_mhz"], row["seed"], label, item["reduced_chi2"], item["multistart_agreement"],
                                 *(item[name][key] for name in QUANTITIES for key in ("value", "se"))])
    (args.out / "REPORT.txt").write_text("\n".join(report_lines(report)) + "\n")
    print("\n".join(report_lines(report)))


if __name__ == "__main__":
    main()
