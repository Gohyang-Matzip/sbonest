"""Plot the measured trade-off from saved, reoptimized profile fits."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
profiles = json.loads((HERE / "profiles.json").read_text())["R1H"]
best = json.loads((HERE / "fit_both.json").read_text())["chi2"]
# Show the low-chi-square range; all grid results remain in profiles.json.
rows = [r for r in profiles if r["parameters"]["R1H"]["value"] <= 20]
x = [r["parameters"]["R1H"]["value"] for r in rows]
y = [r["parameters"]["R2H"]["value"] for r in rows]
delta = np.array([r["chi2"]-best for r in rows])
assert np.all(delta >= -1e-6)
fig, ax = plt.subplots(figsize=(6.6, 4.5), layout="constrained")
ax.plot(x, y, color="0.7", linewidth=1.2, zorder=1)
points = ax.scatter(x, y, c=delta, cmap="viridis", s=70, zorder=2,
                    vmin=0, vmax=1.6, edgecolor="black", linewidth=.4)
ax.scatter([2], [25], marker="*", s=170, color="#D55E00", edgecolor="black",
           linewidth=.4, label="Generating values (2, 25)", zorder=3)
ax.set(xlabel=r"Fixed $R_{1H}$ (s$^{-1}$)", ylabel=r"Refitted $R_{2H}$ (s$^{-1}$)",
       title="A1: proton relaxation trade-off at 1.2 GHz", xlim=(-1,21), ylim=(14.5,27))
ax.text(.03,.05,"All other fit parameters reoptimized\n30 ppm; 25/50/100 Hz; intensity error = 0.001",
        transform=ax.transAxes,fontsize=9)
ax.legend(loc="upper right",frameon=False,fontsize=9)
ax.spines[["top","right"]].set_visible(False)
fig.colorbar(points,ax=ax,label=r"$\Delta\chi^2$ from best fit",pad=.03)
fig.savefig(HERE/"R1H_R2H_tradeoff.png",dpi=220)
fig.savefig(HERE/"R1H_R2H_tradeoff.pdf")
plt.close(fig)
