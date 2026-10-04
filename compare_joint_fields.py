#!/usr/bin/env python3
"""Joint 600 + 800 MHz Sideband fits with field groups, on the archived synthetic benchmark.

Run: .venv/bin/python compare_joint_fields.py \
      --source results/field_comparison_600_800_20261003_02/inputs \
      --out results/field_comparison_joint_NEW [--workers N]

Four joint models are fitted to the six "full" datasets (25/50/100 Hz at each
field) through the normal checkpointed run path and compared with the archived
separate-field SBONEST fits and the generating truth:

- fixed_H_shared_N: R1H = 2, R2H = 25 s^-1 fixed (the original benchmark
  assumption) with one set of nitrogen relaxation rates;
- auto_H_unbounded: field-group proton rates fitted without upper bounds, which
  documents how far unidentifiable proton rates drift at 600/800 MHz;
- auto_H_shared_N: field-group proton rates fitted inside physically generous
  bounds (R1H <= 50, R2H <= 500 s^-1 by default), shared nitrogen relaxation;
- auto_H_per_field_N: the bounded proton rates plus field-group nitrogen rates.

The synthetic truth has identical relaxation at both fields and the fixed-H
values are the true ones, so the extra parameters must not bias kex/pB and AICc
should not prefer them strongly. The result quantifies that, nothing more.
"""
# ruff: noqa: E402 -- Set BLAS limits before importing NumPy.
import argparse
import copy
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import time

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(name, "1")
os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np

ROOT = Path(__file__).resolve().parent
TRUTH = {"kab": 15., "kba": 285., "v1n_scale": 1.08, "peak_ppm": 120., "dw_ppm": 3., "R1": 1.5, "R2a": 12., "R2b": 15.,
         "R1H": 2., "R2H": 25.}
RF_LEVELS = (25, 50, 100)
VARIANTS = {
    "fixed_H_shared_N": {"proton": "fixed", "nitrogen": "shared", "proton_bounds": False},
    "auto_H_unbounded": {"proton": "fit", "nitrogen": "shared", "proton_bounds": False},
    "auto_H_shared_N": {"proton": "fit", "nitrogen": "shared", "proton_bounds": True},
    "auto_H_per_field_N": {"proton": "fit", "nitrogen": "per_field", "proton_bounds": True},
}


def joint_config(folder, variant, out_prefix, proton_bounds=(50., 500.)):
    base = json.loads((folder / "600" / "full.json").read_text())
    config = copy.deepcopy(base)
    config["Project Name"] = str(out_prefix)
    config["datasets"] = [f"{field}/full_{rf}.txt" for field in (600, 800) for rf in RF_LEVELS]
    sb = config["sideband"]
    decoupling = dict(sb["decoupling"])
    if VARIANTS[variant]["proton"] == "fit":
        decoupling.pop("R1H", None)
        decoupling.pop("R2H", None)
        sb["proton_relaxation"] = {"mode": "fit"}
        if VARIANTS[variant]["proton_bounds"]:
            bounds = config["init"].setdefault("bounds", {})
            for entry in config["residues"]:
                if entry.get("flag") == "on":
                    # Ungrouped names address every field group.
                    bounds[f"{entry['name']}.R1H"] = [0., float(proton_bounds[0])]
                    bounds[f"{entry['name']}.R2H"] = [0., float(proton_bounds[1])]
    decoupling.pop("h_larmor_mhz")
    sb["decoupling"] = decoupling
    sb["datasets"] = [{"h_larmor_mhz": float(field)} for field in (600, 800) for _ in RF_LEVELS]
    sb["nitrogen_relaxation"] = {"mode": VARIANTS[variant]["nitrogen"]}
    return config


def information_criteria(chi2, k, n):
    aic = chi2 + 2 * k
    return {"aic": aic, "aicc": aic + 2 * k * (k + 1) / (n - k - 1), "bic": chi2 + k * np.log(n)}


def main():
    from sbfit import run_config

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", type=Path, required=True, help="inputs directory with 600/ and 800/")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--proton-bounds", type=float, nargs=2, default=(50., 500.), metavar=("R1H_MAX", "R2H_MAX"),
                        help="Upper bounds for fitted proton rates in the bounded variants (s^-1)")
    parser.add_argument("--archived-summary", type=Path,
                        help="summary.json of the separate-field study for comparison (default: next to --source)")
    args = parser.parse_args()
    if args.workers < 1:
        parser.error("--workers must be a positive integer")
    args.out.mkdir(parents=True, exist_ok=False)
    inputs = args.out / "inputs"
    sources = [Path(__file__)]
    for field in (600, 800):
        (inputs / str(field)).mkdir(parents=True)
        for name in ["full.json"] + [f"full_{rf}.txt" for rf in RF_LEVELS]:
            shutil.copy2(args.source / str(field) / name, inputs / str(field) / name)
            sources.append(inputs / str(field) / name)
    archived_path = args.archived_summary or (args.source.parent / "summary.json")
    archived = json.loads(archived_path.read_text()) if archived_path.exists() else None
    report = {"scope": "Joint 600 + 800 MHz fits of the archived synthetic 30 ppm benchmark (one noise draw per field, sigma 0.001).",
              "truth": TRUTH, "proton_bounds": list(args.proton_bounds), "variants": {}, "separate_field_reference": None,
              "comparison": ("All variants share kab/kba, the RF scale and the nitrogen shifts across fields and are "
                             "fitted jointly through run_config with checkpoints; AICc/BIC use the supplied absolute sigma."),
              "note": ("The synthetic truth has the same relaxation at both fields and the fixed proton rates are the "
                       "generating values, so field-specific parameters cannot improve the description; the comparison "
                       "measures the cost of the extra freedom, not a benefit.")}
    if archived:
        report["separate_field_reference"] = {
            field: {"kex": entry["fits"]["SBONEST_full"]["best"]["kex"],
                    "kex_se": entry["fits"]["SBONEST_full"]["best"]["derived_se"]["kex"],
                    "pB": entry["fits"]["SBONEST_full"]["best"]["pB"],
                    "pB_se": entry["fits"]["SBONEST_full"]["best"]["derived_se"]["pB"],
                    "reduced_chi2": entry["fits"]["SBONEST_full"]["best"]["reduced_chi2"]}
            for field, entry in archived["fields"].items()}
        sources.append(archived_path)
    for variant in VARIANTS:
        prefix = args.out / variant / "fit"
        prefix.parent.mkdir()
        config = joint_config(inputs, variant, prefix, tuple(args.proton_bounds))
        (prefix.parent / "config.json").write_text(json.dumps(config, indent=2) + "\n")
        began = time.monotonic()
        run_config(copy.deepcopy(config), inputs, True, workers=args.workers)
        result = json.loads(Path(str(prefix) + "_result.json").read_text())
        values = {name: item["value"] for name, item in result["parameters"].items()}
        errors = {name: item["stderr"] for name, item in result["parameters"].items()}
        criteria = information_criteria(result["chi2"], result["n_parameters"], result["n_points"])
        entry = {"settings": VARIANTS[variant], "n_points": result["n_points"], "n_parameters": result["n_parameters"],
                 "chi2": result["chi2"], "reduced_chi2": result["chi2"] / result["dof"], **criteria,
                 "kex": result["kex"], "kex_se": result["derived_se"]["kex"], "pB": result["pB"],
                 "pB_se": result["derived_se"]["pB"], "v1n_scale": values["v1n_scale"],
                 "kex_bias_truth_se": (result["kex"] - 300.) / result["derived_se"]["kex"] if result["derived_se"]["kex"] else None,
                 "parameters": values, "stderr": errors, "field_groups": result["field_groups"],
                 "warnings": result["warnings"], "jacobian_rank": result["jacobian_rank"],
                 "reduced_chi2_diagnostic": result["residual_diagnostics"]["reduced_chi2"],
                 "elapsed_s": time.monotonic() - began, "result_json": str(Path(str(prefix) + "_result.json"))}
        report["variants"][variant] = entry
        print(f"{variant}: k={entry['n_parameters']} chi2/dof={entry['reduced_chi2']:.4f} "
              f"kex={entry['kex']:.3f}+/-{entry['kex_se']:.3f} pB={entry['pB']:.5f} AICc={entry['aicc']:.2f} "
              f"({entry['elapsed_s']:.0f} s)", flush=True)
        (args.out / "summary.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    best = min(report["variants"], key=lambda name: report["variants"][name]["aicc"])
    report["preferred_by_aicc"] = best
    report["delta_aicc_vs_fixed_H_shared_N"] = {name: entry["aicc"] - report["variants"]["fixed_H_shared_N"]["aicc"]
                                                for name, entry in report["variants"].items()}
    report["source_sha256"] = {os.path.relpath(p, ROOT) if p.is_relative_to(ROOT) else str(p):
                               hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    (args.out / "summary.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    with (args.out / "parameters.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["variant", "parameter", "value", "local_SE", "reduced_chi2", "aicc"])
        for name, entry in report["variants"].items():
            for parameter, value in entry["parameters"].items():
                writer.writerow([name, parameter, value, entry["stderr"][parameter], entry["reduced_chi2"], entry["aicc"]])
            writer.writerow([name, "kex", entry["kex"], entry["kex_se"], entry["reduced_chi2"], entry["aicc"]])
            writer.writerow([name, "pB", entry["pB"], entry["pB_se"], entry["reduced_chi2"], entry["aicc"]])
    print(f"Preferred by AICc: {best}; summary at {args.out / 'summary.json'}")


if __name__ == "__main__":
    main()
