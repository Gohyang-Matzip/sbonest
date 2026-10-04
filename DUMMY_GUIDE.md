# SBONEST: step-by-step dummy-data guide

[한국어](DUMMY_GUIDE.ko.md) · [Full manual](SBONEST_MANUAL.md) · [README](README.md)

Start here for your first fit. You will copy the bundled **synthetic** data,
fit three peaks, inspect the results and residuals, and optionally try restarts
and a likelihood scan. No experimental data or neighboring OC checkout is needed.
Commands were checked with Python 3.12 on 2026-10-04.

The example uses one 1.2 GHz proton field, two nitrogen RF amplitudes (25/100 Hz),
and 147 offsets over 105–135 ppm for each of A1, G2 and S3 at each RF amplitude.
It fits all 882 observations and 24 parameters, including peakwise R1H/R2H.
“Two RF” means two amplitudes at the same field. Synthetic recovery checks the
software workflow; it does not validate an experimental sample.

## 1. Open a terminal and install

Use macOS or Linux with Git and Python 3.12 installed (Windows users can use
Linux through WSL). Run these commands in a terminal, not inside a Python prompt.
Copy code blocks without adding a `$` prompt. If you already have this repository,
skip cloning and open its root folder, where `run.py` is located.

```bash
git clone https://github.com/Gohyang-Matzip/sbonest.git
cd sbonest
```

Create the environment below only if `.venv` does not already contain a working
Python 3.12 environment. Do not overwrite another environment in place.
If `.venv` uses a different Python version, follow this guide in a fresh checkout.

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements-sideband.txt -c constraints-sideband.txt
.venv/bin/python --version
```

The last line should show Python 3.12.x. The install includes `optimalcontrol-nmr`.
Keep this terminal open and stay in the repository root for every following step.
No `activate` command is needed because all commands use `.venv/bin/python`.

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg
```

## 2. Make a fresh practice folder

Paste the entire block, including the final `PY` and `)`. It creates a uniquely
named folder under `session_artifacts/`, copies the two input files, and writes
`fit.json`. The originals are preserved. `SBONEST_DEMO_DIR` remembers this folder
in the current terminal; the final command displays its full path.

```bash
SBONEST_DEMO_DIR=$(.venv/bin/python - <<'PY'
import json
import shutil
import tempfile
from pathlib import Path

source = Path("example/sideband_auto_H/two_RF.json").resolve()
config = json.loads(source.read_text())
parent = Path("session_artifacts").resolve()
parent.mkdir(exist_ok=True)
folder = Path(tempfile.mkdtemp(prefix="dummy_", dir=parent))
datasets = []
for name in config["datasets"]:
    original = (source.parent / name).resolve()
    shutil.copy2(original, folder / original.name)
    datasets.append(original.name)
config["datasets"] = datasets
config["Project Name"] = str(folder / "fit")
(folder / "fit.json").write_text(json.dumps(config, indent=2) + "\n")
print(folder)
PY
)
export SBONEST_DEMO_DIR
printf '%s\n' "$SBONEST_DEMO_DIR"
```

You now have `fit.json`, `noisy_25.txt`, and `noisy_100.txt` in that folder.
Dataset paths are relative to the JSON file. `Project Name` is an output prefix,
not a filename ending in `.json`; here it is absolute so results stay together.
The copied configuration is for this bundled example, which has no waveform file.

## 3. Run the fit

```bash
.venv/bin/python run.py "$SBONEST_DEMO_DIR/fit.json"
```

Wait until the command returns to the terminal prompt. The fit prints iterations
and may report weak proton-parameter or boundary warnings. Read those warnings;
they can accompany successful convergence in this noisy example.

For a successful run, the same folder now contains five outputs:

| File | What to do with it |
|---|---|
| `fit_result.json` | Read parameters, uncertainties, diagnostics and provenance |
| `fit_result.txt` | Read the plain-text report; intensity tables are rounded |
| `fit_predictions.csv` | Analyze full-precision predictions and residuals |
| `fit.pdf` | Open the fitted profiles and standardized residual panels |
| `fit_data.pdf` | Open data-only profiles |

To skip PDFs, add `--no-pdf` **before the first run**; only the three numerical
files are then created. To run again, repeat step 2 for a fresh folder and then
step 3. Reusing a completed output prefix is rejected. Keep old folders as records.

## 4. Check that the result makes sense

```bash
.venv/bin/python - <<'PY'
import json
import os
from pathlib import Path

folder = Path(os.environ["SBONEST_DEMO_DIR"])
r = json.loads((folder / "fit_result.json").read_text())
assert r["success"], r.get("message")
assert (r["n_points"], r["n_parameters"], r["dof"]) == (882, 24, 858)
print("kex (s^-1):", r["kex"], "SE:", r["derived_se"]["kex"])
pb_se = r["derived_se"]["pB"]
print("pB (%):", 100 * r["pB"], "SE (percentage points):",
      None if pb_se is None else 100 * pb_se)
print("RF scale:", r["parameters"]["v1n_scale"]["value"])
print("Actual RF (Hz):", r["v1n_hz"])
print("Reduced chi2:", r["chi2"] / r["dof"])
print("Rank:", r["jacobian_rank"], "of", r["n_parameters"])
print("At bounds:", r["at_bounds"])
print("Warnings:", r["warnings"])
print("Provenance fields:", sorted(r["provenance"]))
PY
```

Expect approximately these values; small platform differences are possible.
Reference: [archived two-RF fit](results/auto_H_refit/fits/two_RF_result.json)
and [synthetic generating parameters](results/peakwise_H_fit/summary.json).

| Quantity | Typical fit | Synthetic generating value |
|---|---:|---:|
| kex (s⁻¹) | 299.68 | 300 |
| pB (%) | 4.983 | 5.0 |
| RF scale | 1.08227 | 1.08 |
| Actual RF (Hz), in input-file order | 27.057 / 108.227 | 27 / 108 |
| Reduced χ² | 0.93569 | — |
| Jacobian rank | 24 / 24 | — |

JSON stores `pB` and its SE as **fractions**: 0.05 means 5%. Multiply both by
100 to report percent and percentage-point SE. `derived_se` includes the
covariance between the two exchange rates; do not add their errors independently.
Local SE uses the supplied absolute intensity errors, without reduced-χ² scaling.
`null` means unavailable, and a fixed parameter's zero SE expresses an assumption.
Full rank does not mean every H rate is precise. These are within-model errors;
they exclude incorrect fixed inputs and model mismatch.

## 5. Open the plots and check the CSV

Open the practice folder printed in step 2 using your file manager. Double-click
`fit.pdf` and `fit_data.pdf`; each has one page per peak (three pages here).
If you used `--no-pdf`, continue directly to the CSV check below.
The fit PDF shows full profiles above standardized residuals. Its RF labels are
fitted amplitudes; the data-only PDF labels nominal amplitudes.

Look for persistent residual patterns near both dips, sidebands, and the remote
baseline. A small reduced χ² alone does not establish a correct model. Lines join
predictions at measured offsets; they are not a denser simulated spectrum.

Check that the full-precision CSV reproduces the reported χ²:

```bash
.venv/bin/python - <<'PY'
import csv
import json
import math
import os
from pathlib import Path

folder = Path(os.environ["SBONEST_DEMO_DIR"])
r = json.loads((folder / "fit_result.json").read_text())
with (folder / "fit_predictions.csv").open(newline="") as stream:
    rows = list(csv.DictReader(stream))
chi2 = math.fsum(float(row["residual_sigma"]) ** 2 for row in rows)
assert len(rows) == r["n_points"]
assert math.isclose(chi2, r["chi2"], rel_tol=1e-10, abs_tol=1e-8)
print("CSV rows:", len(rows), "chi2:", chi2)
PY
```

`residual_sigma = (observed - predicted) / sigma`. The CSV also records residue,
dataset index, field, duration, nominal/fitted RF, and offset. Keep the entire
practice folder: the input copies, configuration, and JSON provenance belong with
the outputs. `session_artifacts/` is Git-ignored, so back it up separately if needed.

## 6. Optional: try another start and scan kex

This takes several extra fits. It retains all 24 free parameters and changes no
data. The configured initial fit plus one explicit restart make **two attempts**;
each of three kex grid points then refits the remaining varying parameters.
Use this block once per practice folder; it protects an existing `analysis.json`.

```bash
.venv/bin/python - <<'PY'
import json
import os
from pathlib import Path

folder = Path(os.environ["SBONEST_DEMO_DIR"])
config = json.loads((folder / "fit.json").read_text())
config["Project Name"] = str(folder / "analysis")
config["init"]["multistart"] = {
    "starts": [{"kab": 20.0, "kba": 380.0, "v1n_scale": 1.0}]
}
config["init"]["profile"] = {"kex": [295.0, 300.0, 305.0]}
with (folder / "analysis.json").open("x") as stream:
    json.dump(config, stream, indent=2)
    stream.write("\n")
PY
.venv/bin/python run.py "$SBONEST_DEMO_DIR/analysis.json" --no-pdf
```

Read every attempt and grid-point status:

```bash
.venv/bin/python - <<'PY'
import json
import os
from pathlib import Path

folder = Path(os.environ["SBONEST_DEMO_DIR"])
r = json.loads((folder / "analysis_result.json").read_text())
print("Fit success:", r["success"])
for attempt in r["multistart"]:
    print("Start:", attempt["index"], "success:", attempt["success"],
          "selected:", attempt["selected"], "chi2:", attempt["chi2"],
          "message:", attempt["message"])
if "profiles" in r:
    for point in r["profiles"]["kex"]:
        print("kex:", point["target"], "success:", point["success"],
              "delta chi2:", point["delta_chi2"], "message:", point["message"])
    print("Profile warnings:", r["profiles"]["warnings"])
PY
```

`selected=true` marks the lowest χ² among converged attempts. `success=true` for
the main fit does not guarantee every profile point succeeded. Failed points have
null χ²; negative Δχ² means the scan found a better solution than the selected
baseline. This three-point scan is a demonstration, **not a confidence interval**
or proof of a global minimum. See manual section 8.1 for pB/RF grids and seeded
random starts. If every start fails, the command exits nonzero and saves attempt
records in `analysis_result.json`; inspect them before retrying with a new prefix.

## 7. Optional: measure a single fit

```bash
.venv/bin/python benchmark.py "$SBONEST_DEMO_DIR/fit.json" profile \
  --profile-output "$SBONEST_DEMO_DIR/fit.prof"
.venv/bin/python - <<'PY'
import os
import pstats
from pathlib import Path

path = Path(os.environ["SBONEST_DEMO_DIR"]) / "fit.prof"
pstats.Stats(str(path)).sort_stats("cumulative").print_stats(10)
PY
```

The benchmark runs one fit, saves profiling data, and leaves the fitting outputs
unchanged. It does not run optional multistart/profile-likelihood analyses.
An existing `.prof` file is protected; choose a fresh name for another measurement.
Profiler overhead and the machine affect timing, so do not treat this as an
uninstrumented runtime or compare unlike environments.

## 8. When something goes wrong

| Symptom | Next action |
|---|---|
| `python3.12: command not found` | Install Python 3.12, then repeat step 1 |
| `.venv/bin/python` or `run.py` not found | Return to the repository root and check step 1 |
| `No module named optimalcontrol` | Run step 1's install using the same `.venv/bin/python` |
| `SBONEST_DEMO_DIR` is empty after reopening the terminal | Return to the repo root, repeat the thread export, and set `export SBONEST_DEMO_DIR="/absolute/path/printed/in/step2"` to inspect that run; use step 2 for a new run |
| `Output already exists` or `FileExistsError` | Preserve the old run; repeat step 2 and use the new folder |
| Input file not found after moving JSON | Keep the copied input files beside the JSON, or update its dataset paths |
| Maximum evaluations exceeded | Inspect starts, units and warnings; see manual section 10 before increasing the limit |
| H-rate, rank or boundary warnings | Read uncertainty limits in step 4 and inspect the residuals; convergence alone is insufficient |

## 9. Move on to real data

Use a separate configuration/output folder and follow manual sections 2–6 and 9.
Replace synthetic intensities, positive absolute errors, residue labels, ¹H shifts,
field frequencies, saturation duration and pulse settings with measured values.
The first data header is **¹⁵N MHz**, offsets are absolute ¹⁵N ppm, duration is
seconds, and RF is Hz. Explicitly mark unwanted loaded residues `off`.
Omitting R1H/R2H fits them; supplying both fixes them. Automatic H-rate fitting
currently supports one proton field. Preserve inputs and assess model adequacy
before interpreting experimental exchange parameters.
