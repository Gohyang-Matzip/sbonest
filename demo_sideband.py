#!/usr/bin/env python3
"""Reproduce full fits on independent synthetic NH data; preserve every output.

python demo_sideband.py --out example/sideband_demo
python demo_sideband.py --offset-step 5 --out example/sideband_uniform_5Hz
The output directory must be new. Generates ONEST inputs and runs the real CLI.
"""

# ruff: noqa: E402 -- Set BLAS limits before importing NumPy.
import argparse
import copy
import json
import os
import subprocess
import sys
import time
from pathlib import Path

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(name, "1")
os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np

from sbfit import SidebandModel
from sideband import composite_segments
from test_sideband import reference


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="example/sideband_demo")
    parser.add_argument(
        "--offset-step",
        type=float,
        help="Uniform sampling step in Hz over -2400..2400; must divide 4800 exactly",
    )
    args = parser.parse_args()
    if args.offset_step is not None:
        if not np.isfinite(args.offset_step) or args.offset_step <= 0:
            parser.error("--offset-step must be finite and positive")
        intervals = 4800.0 / args.offset_step
        if intervals < 1 or not np.isclose(
            intervals, round(intervals), rtol=0, atol=1e-8
        ):
            parser.error("--offset-step must divide the 4800 Hz span exactly")
    root = Path(__file__).resolve().parent
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=False)
    sf_n = 1200.0 * 0.101329118
    truth = dict(
        kab=15.0,
        kba=285.0,
        dw=3 * sf_n,
        r1=1.5,
        r2a=12.0,
        r2b=15.0,
        ha=-2760.0,
        hb=-2400.0,
    )
    true_params = {
        "kab": 15.0,
        "kba": 285.0,
        "v1n_scale": 1.08,
        "A1.peak_ppm": 120.0,
        "A1.dw_ppm": 3.0,
        "A1.R1": 1.5,
        "A1.R2a": 12.0,
        "A1.R2b": 15.0,
    }
    offsets = np.unique(
        np.r_[
            np.arange(-2400, 2401, 100),
            np.arange(-100, 101, 8),
            np.arange(270, 451, 12),
            -1720 + np.arange(-42, 43, 6),
            1720 + np.arange(-42, 43, 6),
            -1340 + np.arange(-24, 25, 8),
            1340 + np.arange(-24, 25, 8),
        ]
    )
    if args.offset_step is not None:
        offsets = np.linspace(-2400.0, 2400.0, int(round(intervals)) + 1)
        np.testing.assert_allclose(
            np.diff(offsets), args.offset_step, rtol=0, atol=1e-10
        )
    segments = composite_segments(70e-6, "RR")
    rng = np.random.default_rng(20261002)
    datasets = {key: [] for key in ("noisefree", "noisy", "core")}
    generated = {}
    for nu in (25.0, 50.0, 100.0):
        y = reference(segments, offsets, T=0.4, nu=nu * 1.08, **truth)
        noisy = y + rng.normal(0, 0.001, len(y))
        generated[nu] = noisy
        for tag in datasets:
            pick = (
                np.abs(offsets) < 600 if tag == "core" else np.ones(len(y), dtype=bool)
            )
            intensity = y if tag == "noisefree" else noisy
            path = out / f"{tag}_{nu:g}Hz.txt"
            lines = [
                f"{sf_n:.12g}",
                "0.4",
                f"{nu:g} 0",
                "# offset(ppm) intensity error",
                "# A1 R2a: 13 R2b: 18 dw: 2.7",
            ]
            lines += [
                f"{120 + x / sf_n:.12g} {v:.12g} 0.001"
                for x, v in zip(offsets[pick], intensity[pick])
            ]
            path.write_text("\n".join(lines) + "\n")
            datasets[tag].append(path.name)
    template = {
        "Project Name": "noisefree_scale",
        "datasets": datasets["noisefree"],
        "residues": [{"name": "A1", "flag": "on"}],
        "init": {
            "Method": "Sideband",
            "kex": {"min": 200.0, "max": 400.0, "nsteps": 3},
            "pB": {"min": 0.03, "max": 0.07, "nsteps": 3},
            "initial": {"A1.peak_ppm": 120.05, "A1.R1": 1.2},
        },
        "sideband": {
            "decoupling": {
                "h_larmor_mhz": 1200.0,
                "h_carrier_ppm": 8.5,
                "p90_s": 0.000070,
                "cycle": "RR",
                "b1_scale": 1.0,
                "J_hz": 92.0,
                "R1H": 2.0,
                "R2H": 25.0,
            },
            "residues": {"A1": {"h_ppm_a": 6.2, "h_ppm_b": 6.5}},
            "v1n": {"mode": "scale", "initial": 1.0, "bounds": [0.8, 1.2]},
        },
    }
    cases = [
        ("noisefree_scale", "noisefree", "scale"),
        ("noisy_scale", "noisy", "scale"),
        ("noisy_fixed", "noisy", "fixed"),
        ("noisy_core", "core", "scale"),
        ("noisy_per_dataset", "noisy", "per_dataset"),
    ]
    summaries = {}
    for name, data, mode in cases:
        cfg = copy.deepcopy(template)
        cfg["Project Name"], cfg["datasets"] = name, datasets[data]
        if mode != "scale":
            cfg["sideband"]["v1n"] = {"mode": mode}
        path = out / f"{name}.json"
        path.write_text(json.dumps(cfg, indent=2) + "\n")
        start = time.perf_counter()
        with (out / f"{name}.log").open("w") as log:
            result = subprocess.run(
                [sys.executable, str(root / "run.py"), str(path), "--no-pdf"],
                cwd=out,
                stdout=log,
                stderr=subprocess.STDOUT,
            )
        if result.returncode:
            raise RuntimeError(f"{name} failed; see {out / (name + '.log')}")
        summaries[name] = json.loads((out / f"{name}_result.json").read_text())
        summaries[name]["runtime_s"] = time.perf_counter() - start
        fit = summaries[name]
        print(
            f"{name}: chi2/dof={fit['chi2'] / fit['dof']:.5g}, kex={fit['kex']:.5g}, "
            f"pB={fit['pB']:.5g}, {fit['runtime_s']:.1f} s",
            flush=True,
        )
    exact = summaries["noisefree_scale"]
    assert exact["chi2"] < 1e-8
    for key, true in true_params.items():
        assert abs(exact["parameters"][key]["value"] - true) < 1e-3, (
            key,
            exact["parameters"][key],
        )
    # Covariance is checked against independent noise; it is not a proof for real data.
    noisy = summaries["noisy_scale"]
    assert abs(noisy["parameters"]["v1n_scale"]["value"] - 1.08) < 0.02
    assert 0.5 < noisy["chi2"] / noisy["dof"] < 1.5
    (out / "validation.json").write_text(
        json.dumps(
            {
                "generator": (
                    "Independent complex density matrix (test_sideband.reference)"
                ),
                "seed": 20261002,
                "noise_sigma": 0.001,
                "offset_grid": {
                    "kind": "uniform" if args.offset_step is not None else "adaptive",
                    "step_hz": args.offset_step,
                    "min_hz": float(offsets[0]),
                    "max_hz": float(offsets[-1]),
                    "points_per_spectrum": len(offsets),
                },
                "truth": true_params,
                "fits": summaries,
            },
            indent=2,
            allow_nan=False,
        )
        + "\n"
    )
    plot(out, offsets, generated, template, summaries)
    print(f"PASS: full-parameter recovery and actual CLI; outputs in {out}")


def plot(out, offsets, observed, template, results):
    import matplotlib.pyplot as plt

    models = {}
    for tag in ("noisy_scale", "noisy_fixed"):
        cfg = json.loads((out / f"{tag}.json").read_text())
        m = SidebandModel(cfg, out)
        p = np.array(
            [results[tag]["parameters"][n]["value"] for n in m.parameter_names]
        )
        models[tag] = m, m.seParam(p)
    fig, axes = plt.subplots(3, 3, figsize=(12, 8), constrained_layout=True)
    windows = [(-150, 500), (-1770, -1670), (1670, 1770)]
    for row, nu in enumerate(observed):
        for col, (lo, hi) in enumerate(windows):
            ax = axes[row, col]
            use = (offsets >= lo) & (offsets <= hi)
            ax.errorbar(
                offsets[use],
                observed[nu][use],
                yerr=0.001,
                fmt="o",
                ms=3,
                color="black",
                alpha=0.65,
                label="Synthetic data (σ=0.001)",
            )
            x = np.linspace(lo, hi, 250)
            for tag, color, label in [
                ("noisy_scale", "#0072B2", "Fit v1n"),
                ("noisy_fixed", "#D55E00", "Fix v1n to nominal"),
            ]:
                m, P = models[tag]
                es = m.dataset.res[0].estSpecs[row]
                ax.plot(
                    x,
                    m.calc(P, 0, 120 + x / es.field, es),
                    color=color,
                    lw=1.4,
                    label=label,
                )
            ax.set_xlim(lo, hi)
            ax.set_ylabel(f"{nu:g} Hz nominal\nI / I₀")
            if row == 0:
                ax.set_title(
                    ["Main and minor dips", "Negative sideband", "Positive sideband"][
                        col
                    ]
                )
            if row == 2:
                ax.set_xlabel("¹⁵N RF offset from site A (Hz)")
            ax.spines[["top", "right"]].set_visible(False)
    axes[0, 0].legend(fontsize=8)
    fig.suptitle(
        "SBONEST · independent synthetic validation · true v1n = 1.08 × nominal"
    )
    fig.savefig(out / "sideband_fit.png", dpi=170)
    fig.savefig(out / "sideband_fit.pdf")
    plt.close(fig)


if __name__ == "__main__":
    main()
