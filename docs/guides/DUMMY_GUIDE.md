# SBONEST: step-by-step dummy-data guide

[한국어](DUMMY_GUIDE.ko.md) · [Full manual](../manual/SBONEST_MANUAL.md) · [README](../../README.md)

Start here for your first fit. You will copy the bundled **synthetic** data,
fit three peaks, inspect the results and residuals, and optionally try restarts,
a likelihood scan, and a small bootstrap. No experimental data or neighboring OC
checkout is needed. Updated for the Python 3.12 workflow on 2026-10-04.

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

The helper copies the two bundled inputs into a **new** folder and writes
`fit.json`, preserving the originals. Keep `SBONEST_DEMO_DIR` set in this terminal.
If `dummy_01` already exists, change it to `dummy_02` (or another unused name).

```bash
export SBONEST_DEMO_DIR="$PWD/session_artifacts/dummy_01"
.venv/bin/python sb_workflow.py init-demo --out "$SBONEST_DEMO_DIR"
```

The folder now contains `fit.json`, `data/noisy_25.txt`, and `data/noisy_100.txt`.
Dataset paths resolve from the JSON directory. `Project Name` is an output prefix,
not a filename ending in `.json`; the helper makes it absolute so results stay
together. This helper uses the bundled example, which has no waveform file.

Check the inputs and planned outputs before fitting:

```bash
.venv/bin/python run.py "$SBONEST_DEMO_DIR/fit.json" --check
```

The JSON summary should contain `valid: true`, 882 points and 24 free parameters.
It lists resolved inputs, fields/RF, free/fixed parameters, bounds and output
conflicts. It writes no files and runs no optimizer; an initialization grid can
still evaluate the model. A valid check does not guarantee convergence.
Adding `--identifiability` to the check evaluates one Jacobian at the initial
values and lists weakly determined parameters and strong correlations; for this
example it already names the three `R1H` rates and their `R1H/R2H` pairs.

## 3. Run the fit

```bash
.venv/bin/python run.py "$SBONEST_DEMO_DIR/fit.json"
```

Wait until the command returns to the terminal prompt. The fit prints iterations
and may report weak proton-parameter or boundary warnings. Read those warnings;
they can accompany successful convergence in this noisy example.

A successful run adds these outputs and a checkpoint directory:

| File | What to do with it |
|---|---|
| `fit_result.json` | Read parameters, uncertainties, diagnostics and provenance |
| `fit_result.txt` | Read the plain-text report; intensity tables are rounded |
| `fit_predictions.csv` | Analyze full-precision predictions and residuals |
| `fit.pdf` | Open the fitted profiles and standardized residual panels |
| `fit_data.pdf` | Open data-only profiles |
| `fit_checkpoint/` | Preserve the run identity and completed fitting/analysis records |

To skip PDFs, use `--no-pdf` for both the check and the first fit; the three
numerical files and checkpoint are still created. To start a different run,
repeat step 2 with an unused name such as `dummy_02`. Preserve previous folders.
On a multi-core computer, add `--workers 4` (or the number of cores you can
spare) to the fit command; the results are identical and the fit, restarts,
scans and bootstrap replicates finish sooner. Keep the thread variables from
step 1 at one.

After an interruption, resume the unchanged run with:

```bash
.venv/bin/python run.py "$SBONEST_DEMO_DIR/fit.json" --resume
```

If the first fit used `--no-pdf`, add it to the resume command too. Resume requires
the same configuration, input/waveform files, executed source, Python/packages,
platform, thread settings and PDF mode. It reuses the completed baseline,
restarts, profile points and bootstrap replicates. Unfinished work continues;
existing unrelated output files are protected. Keep `fit_checkpoint/` intact;
its records are checksum-validated and must not be edited. An interrupted fit
can continue, but a fit that already finished with failure is restored and exits
nonzero without retrying. Correct the settings and use a fresh `Project Name`
to retry a failed fit. Changing the environment also requires a fresh prefix.

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
Reference: [archived two-RF fit](../../results/auto_H_refit/fits/two_RF_result.json)
and [synthetic generating parameters](../../results/peakwise_H_fit/summary.json).

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

## 5. Inspect plots and regenerate a saved-result report

Open the practice folder printed in step 2 using your file manager. Double-click
`fit.pdf` and `fit_data.pdf`; each has one page per peak (three pages here).
If you used `--no-pdf`, the report command below can still create a PDF.
The fit PDF shows full profiles above standardized residuals. Its RF labels are
fitted amplitudes; the data-only PDF labels nominal amplitudes.

Look for persistent residual patterns near both dips, sidebands, and the remote
baseline. A small reduced χ² alone does not establish a correct model. Lines join
predictions at measured offsets; they are not a denser simulated spectrum.

Create a report directly from the saved JSON and matching predictions CSV:

```bash
.venv/bin/python sb_workflow.py report "$SBONEST_DEMO_DIR/fit_result.json" \
  --out "$SBONEST_DEMO_DIR/report_01"
```

This checks the CSV residual calculations and their agreement with the JSON
point count, χ² and fitted RF. It creates `report_01_summary.json`, `report_01_summary.txt`
and `report_01.pdf`, with residual statistics by residue and dataset. It does not
refit or need the original fitting inputs. Keep the result JSON and predictions
CSV together; if the CSV was renamed, pass its path with `--predictions`.
Existing report files are protected; use `report_02` to regenerate again.

`residual_sigma = (observed - predicted) / sigma`. The CSV also records residue,
dataset index, field, duration, nominal/fitted RF, and offset. Keep the entire
practice folder: the input copies, configuration, and JSON provenance belong with
the outputs. `session_artifacts/` is Git-ignored, so back it up separately if needed.

## 6. Optional: restarts, a kex scan and bootstrap

This runs extra fits with all 24 free parameters: **two attempts** (the configured
initial fit and one restart), three constrained kex refits, and five parametric
bootstrap replicates. The original data remain unchanged. Use this block once per
practice folder; it protects an existing `analysis.json` and uses a new prefix.

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
config["init"]["bootstrap"] = {"replicates": 5, "seed": 20261004, "confidence": 0.95}
with (folder / "analysis.json").open("x") as stream:
    json.dump(config, stream, indent=2)
    stream.write("\n")
PY
.venv/bin/python run.py "$SBONEST_DEMO_DIR/analysis.json" --check --no-pdf
.venv/bin/python run.py "$SBONEST_DEMO_DIR/analysis.json" --no-pdf
```

Create and read a report containing every attempt, grid point and replicate:

```bash
.venv/bin/python sb_workflow.py report "$SBONEST_DEMO_DIR/analysis_result.json" \
  --out "$SBONEST_DEMO_DIR/analysis_report_01"
```

Open `analysis_report_01_summary.txt` or its PDF. `selected=true` marks the lowest
χ² among converged starts. Main-fit success does not guarantee every profile point
or bootstrap replicate succeeded. Failed profile points retain null χ²; negative
Δχ² means the scan found a better solution than the selected baseline. Three grid
points do not define a confidence interval or prove a global minimum.

Bootstrap uses seeded Gaussian draws around the selected fit with the supplied
absolute sigma, then refits each draw. The result's `bootstrap.samples` retains
successes and failures; `bootstrap.intervals` contains successful-sample percentiles.
Five replicates only demonstrate the workflow: fewer than 100 successful replicates
trigger a tail-estimate warning. These intervals condition on the model and fixed
inputs, exclude model mismatch, and have no guaranteed confidence-interval coverage.
See manual sections 8.1–8.2 for additional settings and synthetic local-SE coverage.

If interrupted, use `--resume --no-pdf` with `analysis.json` and unchanged settings.
If every start fails, the command exits nonzero and preserves the attempt records
in `analysis_result.json` and the checkpoint. `--resume` restores that failure
and exits nonzero without repeating the attempts. Inspect the records, correct
the settings and choose a fresh prefix; a failure result has no prediction CSV
to report.

Three more optional tools use the same files. Adding
`"profile_interval": {"parameters": ["kex"]}` to `init` in a fresh copy of the
configuration locates the 95% likelihood-ratio interval of kex on the profile
(each evaluation is one refit; `--workers` helps). To ask before measuring how
many offsets or RF levels a planned experiment needs, write a design file as in
manual section 8.6 and run `sb_workflow.py design DESIGN_JSON --out NEW_DIR`;
it fits nothing and reports expected standard errors per scenario. To test
whether the three peaks share one exchange process, run
`sb_workflow.py compare "$SBONEST_DEMO_DIR/fit.json" --out "$SBONEST_DEMO_DIR/compare_01"`;
it fits the shared model and each peak alone and compares them with AICc and an
F-test (manual section 8.7). Expect several full fits.

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
unchanged. It does not run optional multistart, profile-likelihood or bootstrap analyses.
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
| Input file not found after moving JSON | Keep the copied `data/` directory with the JSON, or update its dataset paths |
| Checkpoint identity mismatch | Restore the original settings/environment to resume, or preserve the run and use a fresh prefix |
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

Workflow sources: [sb_workflow.py](../../sb_workflow.py), [sbfit.py](../../sbfit.py),
[sb_checkpoint.py](../../sb_checkpoint.py), [sb_report.py](../../sb_report.py), and
[sb_bootstrap.py](../../sb_bootstrap.py).
