"""Compare SBONEST and RF-matched ONEST Matrix on archived 600/800 MHz data.

Run: .venv/bin/python compare_field_models.py --out results/field_comparison_600_800
Use a new output directory; inputs and existing results are never overwritten.
"""

# ruff: noqa: E402 -- Set BLAS limits before importing NumPy.
import argparse
import copy
import csv
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import shutil
import sys
import time

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(name, "1")
os.environ.setdefault("MPLBACKEND", "Agg")

ROOT = Path(__file__).resolve().parent
STUDY = ROOT / "manuscript/sideband_30ppm"
sys.path.insert(0, str(STUDY))

import matplotlib.pyplot as plt
import numpy as np

from compare_conventional import NitrogenModel, TRUTH, metrics, prediction
from fit import generate_initial_parameters
from sbfit import SidebandModel
from test_sideband import reference

MODELS = ("SBONEST", "ONEST_Matrix", "N_signed_control")


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def fit_starts(config, folder, label):
    model = (SidebandModel if label == "SBONEST" else NitrogenModel)(config, folder)
    if label != "SBONEST":
        model.clipped = label == "ONEST_Matrix"
    model._initializing = True
    try:
        standard = model._expand_initial(generate_initial_parameters(model, config["init"]))
    finally:
        model._initializing = False
    for key, value in config["init"].get("initial", {}).items():
        standard[model.parameter_names.index(key)] = value
    fit_config = copy.deepcopy(config["init"])
    fit_config["initial"] = {}
    starts = {"standard_grid": standard, "truth": TRUTH,
              "perturbed": np.array([24., 376., .95, 120.08, 2.7, 1.2, 18., 25.])}
    attempts = []
    for name, initial in starts.items():
        began = time.monotonic()
        p, cov = model.fit(initial, fitting_config=fit_config)
        record = model.diagnostics(p, cov)
        record.pop("config")
        assert record["jacobian_rank"] == record["n_parameters"] == 8
        assert np.isfinite(cov).all()
        gradients = np.zeros((3, len(p)))
        gradients[0, :2] = 1
        gradients[1, :2] = [p[1], -p[0]] / (p[0] + p[1])**2
        gradients[2, [3, 4]] = 1
        derived = dict(zip(("kex", "pB", "delta_B"),
                           np.sqrt(np.diag(gradients @ cov @ gradients.T)).tolist()))
        record.update(method=label, start=name, initial=initial.tolist(),
                      parameter_vector=p.tolist(), parameter_order=model.parameter_names,
                      covariance=cov.tolist(), derived_se=derived, metrics=metrics(p),
                      reduced_chi2=record["chi2"] / record["dof"],
                      elapsed_s=time.monotonic() - began)
        attempts.append(record)
        print(folder.name, Path(config["datasets"][0]).stem, label, name,
              f"chi2/dof={record['reduced_chi2']:.7g} kex={record['kex']:.7g}", flush=True)
    best = min(attempts, key=lambda row: row["chi2"])
    spread = max(row["chi2"] for row in attempts) - best["chi2"]
    agreement = spread < max(1e-6, best["chi2"] * 1e-6)
    if not agreement:
        print(f"NOTE: {label} has distinct attained minima; chi2 spread={spread:.7g}", flush=True)
    return model, np.array(best["parameter_vector"]), {
        "best": best, "attempts": attempts, "multistart_chi2_spread": spread,
        "multistart_agreement": agreement}


def verify_inputs(config, folder):
    nh = SidebandModel(config, folder)
    full = nh._prepare_data()
    errors = []
    for _, es, offsets, observed, sigma in full:
        exact = np.loadtxt(folder / f"exact_{es.v1:g}.txt", skiprows=5)
        masked = np.loadtxt(folder / f"masked_{es.v1:g}.txt", skiprows=5)
        np.testing.assert_array_equal(masked, np.c_[offsets, observed, sigma])
        np.testing.assert_array_equal(exact[:, 0], offsets)
        np.testing.assert_allclose(sigma, .001, rtol=0, atol=0)
        hfield = config["sideband"]["decoupling"]["h_larmor_mhz"]
        independent = reference(nh.segments[es.dataset_index], (offsets - 120) * es.field,
                                T=es.T, nu=es.v1 * 1.08, kab=15, kba=285,
                                dw=3 * es.field, r1=1.5, r2a=12, r2b=15,
                                ha=-2.3 * hfield, hb=-2 * hfield)
        np.testing.assert_allclose(exact[:, 1], independent, atol=2e-10, rtol=0)
        errors.append(float(np.max(abs(exact[:, 1] - independent))))
    zero_config = copy.deepcopy(config)
    zero_config["sideband"]["decoupling"]["J_hz"] = 0
    zero = SidebandModel(zero_config, folder)
    nitrogen = NitrogenModel(config, folder)
    controls = []
    for p in (TRUTH, np.array([24., 376., .95, 120.08, 2.7, 1.2, 18., 25.])):
        nitrogen.clipped = False
        signed = prediction(nitrogen, p, full)
        nh_zero = prediction(zero, p, full)
        nitrogen.clipped = True
        clipped = prediction(nitrogen, p, full)
        np.testing.assert_allclose(signed, nh_zero, atol=1e-10, rtol=0)
        np.testing.assert_allclose(clipped, np.maximum(signed, 0), atol=1e-10, rtol=0)
        controls.append({"parameters": p.tolist(),
                         "J0_max_abs_difference": float(abs(signed - nh_zero).max()),
                         "ONEST_clamp_parity": float(abs(clipped - np.maximum(signed, 0)).max())})
    return full, {"independent_generator_max_errors": errors, "J0_controls": controls,
                  "mask_excluded_points": 0}


def make_figure(field, arrays, out):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 3, figsize=(11, 6), sharex=True, sharey="row",
                             layout="constrained", gridspec_kw={"height_ratios": [1.5, 1]})
    for column, rf in enumerate((25, 50, 100)):
        top, bottom = axes[:, column]
        top.plot(arrays["ppm"], arrays["observed"][column], ".", color=".45",
                 ms=3, label="Synthetic observations")
        for label, color, style in (("SBONEST", "#0072B2", "-"),
                                     ("ONEST_Matrix", "#D55E00", "--")):
            top.plot(arrays["dense_ppm"], arrays[label + "_dense"][column],
                     color=color, ls=style, lw=1.3, label=label.replace("_", " "))
            z = (arrays[label + "_prediction"][column] - arrays["observed"][column]) / .001
            bottom.plot(arrays["ppm"], z, marker="o", ms=2, lw=.8, color=color, ls=style)
        top.set_title(f"{'ABC'[column]}   Nitrogen RF = {rf} Hz", loc="left", weight="bold")
        bottom.set_title(f"{'DEF'[column]}   Residuals", loc="left", weight="bold")
        bottom.axhline(0, color=".3", lw=.7)
        bottom.axhspan(-2, 2, color=".93", zorder=0)
        bottom.set_xlabel("¹⁵N irradiation position (ppm)")
        bottom.set_xlim(105, 135)
        bottom.set_xticks([105, 115, 125, 135])
    axes[0, 0].set_ylabel("Normalized intensity (I/I₀)")
    axes[1, 0].set_ylabel("(Prediction − observation) / σ")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=3, frameon=False)
    fig.suptitle(f"{field} MHz · Same {arrays['observed'].size} observations · 8 parameters per model")
    for extension in ("png", "pdf"):
        fig.savefig(out / f"comparison_{field}.{extension}", dpi=240, bbox_inches="tight")
    plt.close(fig)
    save_json(out / f"comparison_{field}_caption.json", {
        "claim": "Compare model predictions and residuals on identical synthetic observations.",
        "panel_plan": {"A–C": "25/50/100 Hz profiles and both fits",
                       "D–F": "Standardized residuals on the acquired grid; gray band is ±2 sigma"},
        "caption": "One synthetic noise realization, sigma=0.001. Both models fit a shared RF scale. "
                   "SBONEST fixes R1H=2 and R2H=25 s^-1. ONEST uses its actual clipped Matrix worker "
                   "with the same optimizer and bounds. All acquired points are retained.",
        "data": f"predictions_{field}.npz", "fits": "summary.json",
        "generator": "compare_field_models.py"})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=STUDY / "results")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    report = {"scope": "600 and 800 MHz, existing synthetic 30 ppm benchmark; one noise draw per field.",
              "truth": metrics(TRUTH), "sigma": .001, "fields": {},
              "comparison": "Same data, weights, eight free parameters, bounds and SidebandModel.fit optimizer. "
                            "ONEST_Matrix uses estmodel.matrix_calc_worker_chunk including its zero clamp. "
                            "N_signed_control disables only that clamp. RF scale is fitted for all models; "
                            "this is a matched forward-model comparison, not the unmodified ONEST CLI.",
              "uncertainty": "Absolute-sigma local 1 SE, retaining covariance for derived parameters; "
                             "no reduced-chi-square scaling. Misspecified-model SE excludes model discrepancy.",
              "optimization": "Best objective among three specified starts; all attempts retained. "
                              "Distinct attained minima are reported, not treated as successful global convergence.",
              "proton_relaxation": "Fixed at R1H=2, R2H=25 s^-1, matching the original manuscript benchmark; "
                                   "not the newer automatic proton-relaxation workflow.",
              "versions": {name: importlib.metadata.version(name)
                           for name in ("numpy", "scipy", "matplotlib", "optimalcontrol-nmr")}}
    sources = [Path(__file__), STUDY / "compare_conventional.py", ROOT / "sbfit.py",
               ROOT / "sideband.py", ROOT / "estmodel.py", ROOT / "est_data.py",
               ROOT / "fit.py", ROOT / "test_sideband.py"]
    for field in (600, 800):
        folder = args.out / "inputs" / str(field)
        folder.mkdir(parents=True)
        for name in ["full.json"] + [f"{kind}_{rf}.txt" for kind in ("full", "exact", "masked")
                                      for rf in (25, 50, 100)]:
            source = args.source / str(field) / name
            shutil.copy2(source, folder / name)
            sources.append(folder / name)
        config = json.loads((folder / "full.json").read_text())
        full, checks = verify_inputs(config, folder)
        entry = {"checks": checks, "config": config, "fits": {}}
        report["fields"][str(field)] = entry
        arrays = {"ppm": full[0][2], "observed": np.array([row[3] for row in full]),
                  "dense_ppm": np.linspace(105, 135, 1201)}
        for kind in ("full", "exact"):
            cfg = copy.deepcopy(config)
            cfg["datasets"] = [f"{kind}_{rf}.txt" for rf in (25, 50, 100)]
            for label in MODELS:
                model, p, result = fit_starts(cfg, folder, label)
                entry["fits"][f"{label}_{kind}"] = result
                if label == "SBONEST" and kind == "exact":
                    np.testing.assert_allclose(p, TRUTH, atol=1e-4, rtol=1e-6)
                    assert result["best"]["chi2"] < 1e-8
                if kind == "full":
                    arrays[label + "_prediction"] = prediction(model, p, full)
                    arrays[label + "_dense"] = prediction(model, p, full, arrays["dense_ppm"])
                save_json(args.out / "summary.json", report)
        np.savez(args.out / f"predictions_{field}.npz", **arrays)
        make_figure(field, arrays, args.out)
    report["source_sha256"] = {os.path.relpath(p, ROOT): hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in sources}
    save_json(args.out / "summary.json", report)
    with (args.out / "parameters.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["field_MHz", "fit", "parameter", "value", "local_SE", "reduced_chi2"])
        for field, entry in report["fields"].items():
            for label, result in entry["fits"].items():
                best = result["best"]
                for name, record in best["parameters"].items():
                    writer.writerow([field, label, name, record["value"], record["stderr"], best["reduced_chi2"]])
                for name in ("kex", "pB"):
                    writer.writerow([field, label, name, best[name], best["derived_se"][name], best["reduced_chi2"]])
    print("PASS: input/generator parity, J=0 controls, exact NH recovery, full-rank covariance. "
          "Multistart agreement is reported separately for each fit.")


if __name__ == "__main__":
    main()
