"""Plot saved demo fits over the complete sampled offset range; no refitting."""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from sbfit import SidebandModel


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "folder",
        nargs="?",
        type=Path,
        default=Path(__file__).resolve().parent / "example/sideband_demo",
    )
    folder = parser.parse_args().folder.resolve()
    prefix = folder / "sideband_full_profile"
    if any(prefix.with_suffix(ext).exists() for ext in (".png", ".pdf", ".npz")):
        raise FileExistsError("Full-profile outputs exist; archive before regenerating")
    models = {}
    for name in ("noisy_scale", "noisy_fixed"):
        config = json.loads((folder / f"{name}.json").read_text())
        saved = json.loads((folder / f"{name}_result.json").read_text())
        model = SidebandModel(config, folder)
        p = np.array(
            [saved["parameters"][key]["value"] for key in model.parameter_names]
        )
        residuals = model.errFunc(p)
        np.testing.assert_allclose(residuals @ residuals, saved["chi2"], rtol=1e-8)
        models[name] = model, model.seParam(p)

    model = models["noisy_scale"][0]
    spectra = model.dataset.res[0].estSpecs
    validation = json.loads((folder / "validation.json").read_text())
    reference_ppm = validation["truth"]["A1.peak_ppm"]
    grid = validation.get("offset_grid", {})
    sampling = (
        f" · uniform {grid['step_hz']:g} Hz sampling"
        if grid.get("kind") == "uniform"
        else ""
    )
    observed_offsets = [
        (np.asarray(es.offset) - reference_ppm) * es.field for es in spectra
    ]
    lo = round(min(x.min() for x in observed_offsets))
    hi = round(max(x.max() for x in observed_offsets))
    x = np.arange(lo, hi + 1.0, 1.0)
    fig, axes = plt.subplots(
        len(spectra), 1, figsize=(12, 8), sharex=True, sharey=True, layout="constrained"
    )
    curves = {name: [] for name in models}
    for row, (ax, es, observed_x) in enumerate(zip(axes, spectra, observed_offsets)):
        for name, color, style, label in (
            ("noisy_scale", "#0072B2", "-", "Fit ν₁N"),
            ("noisy_fixed", "#D55E00", "--", "Fix ν₁N to nominal"),
        ):
            m, params = models[name]
            spectrum = m.dataset.res[0].estSpecs[row]
            y = m.calc(params, 0, reference_ppm + x / spectrum.field, spectrum)
            assert y.shape == x.shape and np.isfinite(y).all()
            curves[name].append(y)
            ax.plot(x, y, color=color, ls=style, lw=1.3, label=label, zorder=3)
        ax.errorbar(
            observed_x,
            es.int,
            yerr=es.intstd,
            fmt="o",
            ms=1.6 if grid.get("kind") == "uniform" else 2.6,
            color="#222222",
            alpha=0.4 if grid.get("kind") == "uniform" else 0.65,
            elinewidth=0.65,
            label="Synthetic data (σ = 0.001)",
        )
        ax.set_ylabel("I / I₀", fontsize=12)
        ax.text(
            0.015,
            0.09,
            f"Nominal ν₁N = {es.v1:g} Hz",
            transform=ax.transAxes,
            fontsize=11,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.9),
        )
        ax.set_xlim(lo, hi)
        ax.set_ylim(-0.025, 0.60)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="x", color=".9", lw=0.6)
        ax.tick_params(labelsize=10)
    axes[0].legend(loc="lower right", fontsize=9, frameon=False)
    axes[-1].set_xlabel("¹⁵N RF offset from site A (Hz)", fontsize=12)
    axes[-1].set_xticks(np.arange(-2000, 2001, 500))
    fig.suptitle(
        "Full CEST profiles · 90x–240y–90x (RR)\n"
        f"Independent synthetic data · p90 = 70 μs · T = 0.4 s{sampling}",
        fontsize=14,
    )
    fig.savefig(prefix.with_suffix(".png"), dpi=180)
    fig.savefig(prefix.with_suffix(".pdf"))
    plt.close(fig)
    np.savez(
        prefix.with_suffix(".npz"),
        offsets_hz=x,
        fitted=np.array(curves["noisy_scale"]),
        fixed=np.array(curves["noisy_fixed"]),
        observed_offsets_hz=np.array(observed_offsets),
        observed=np.array([es.int for es in spectra]),
        std=np.array([es.intstd for es in spectra]),
        nominal_v1n_hz=np.array([es.v1 for es in spectra]),
        reference_ppm=reference_ppm,
    )
    print(f"Saved {prefix}.png/.pdf/.npz; {lo:g} to {hi:g} Hz; saved-fit chi2 verified")


if __name__ == "__main__":
    main()
