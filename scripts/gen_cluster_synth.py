"""Generate synthetic ONEST CEST datasets with known cluster structure (ground truth K=1..4)."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]  # repository root: flat modules live there
sys.path.insert(0, str(ROOT))

import numpy as np

sys.path.insert(0, "/Users/donghanlee/work/projects/ONEST")
from estmodel import est_model

OUT = "."
B0, T, V1, PB, R1, R2A, R2B, SIG = 80.12, 0.4, 25.0, 0.05, 1.0, 10.0, 20.0, 0.01
DGS = [110.0, 112.0, 114.0, 116.0, 118.0, 120.0, 122.0, 124.0]
DWS = [-5.0, 3.0, -3.5, 4.0, -2.5, 3.5, -4.0, 2.5]

SCENARIOS = {  # name -> (seed, kex per residue B1..B8)
    "k1": (1, [300] * 8),
    "k2": (2, [150] * 4 + [600] * 4),
    "k3": (3, [100] * 3 + [350] * 3 + [1000] * 2),
    "k4": (4, [100, 100, 300, 300, 700, 700, 1400, 1400]),
}

m = est_model()
for name, (seed, kexs) in SCENARIOS.items():
    rng = np.random.default_rng(seed)
    lines = [f"{B0}", f"{T}", f"{V1} 0.0", "#offset(ppm) Intensity error"]
    for i, (dG, dw, kex) in enumerate(zip(DGS, DWS, kexs)):
        kab = PB * kex
        kba = kex - kab
        offs = np.linspace(dG - 8, dG + 8, 200)
        ical = m.Baldwin(kab, kba, dG, dw, R1, R2A, R2B, offs, T, V1, 0.0, B0)
        noisy = ical + rng.normal(0, SIG, len(offs))
        # full residue format so fits start from sane R2a/R2b/dw (as in real usage, e.g. syn10)
        lines.append(f"# B{i + 1} R2a: {R2A} R2b: {R2B} dw: {dw}")
        lines += [f"{o:9.4f} {y:9.4f} {SIG:9.4f}" for o, y in zip(offs, noisy)]
    with open(f"{OUT}/{name}v2.txt", "w") as f:
        f.write("\n".join(lines) + "\n")
    print(name, "written: kex =", kexs)
