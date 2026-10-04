# SBONEST fitting manual

[한국어](SBONEST_MANUAL.ko.md) · [Step-by-step dummy guide](DUMMY_GUIDE.md) ·
[README](README.md) · [Technical notes](SIDEBAND.md)

This manual covers the command-line `Sideband` model in this checkout, updated on
2026-10-04. It fits two-state ¹⁵N CEST profiles, including ¹H decoupling sidebands,
using NH spin operators from OC. The main example uses repeated
90°x–240°y–90°x decoupling and a complete 30 ppm offset window at 1.2 GHz.
All supplied demonstration data are **synthetic**, not experimental measurements.

**R1H and R2H are optional inputs.** Omit both from the decoupling settings to
fit them automatically and independently for each peak. See the
[new full-profile refit](results/auto_H_refit/fits/REFIT_REPORT.md) and its
[portable two-RF configuration](example/sideband_auto_H/two_RF.json), which needs no H-rate inputs.

The inherited [ONEST manual](MANUAL.md) describes other models. Its web interface,
`prepare.py`, and `mcrun.py` are not the supported workflow for `Sideband` fitting.

## 1. Environment and first fit

Use Python 3.12 (the CI target). Install the dependencies in a new checkout:

For a first run with input copies, a new output directory, and result
checks, use the [dummy-data guide](DUMMY_GUIDE.md). The shorter commands below
use the example's fixed output prefix and are intended for its first execution.

```bash
git clone https://github.com/Gohyang-Matzip/sbonest.git
cd sbonest
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements-sideband.txt -c constraints-sideband.txt
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg
.venv/bin/python run.py example/sideband_auto_H/two_RF.json --check
.venv/bin/python run.py example/sideband_auto_H/two_RF.json
```

OC is installed as `optimalcontrol-nmr`; no neighboring OC checkout is needed.
Reuse a working environment when one already exists. Add `--no-pdf` to both the
check and fit commands if only numerical outputs are needed.

`--check` prints a JSON summary with resolved dataset/waveform paths, point counts,
fields/RF, free/fixed parameters, initial values, bounds and output conflicts.
It exits nonzero for invalid settings or conflicts. It writes no files and does
not run an optimizer, although an initialization grid may evaluate the model.
Optional analyses are validated before fitting. This checks preparation, not
convergence or identifiability; do it before creating a new run's outputs.

`--check --identifiability` additionally evaluates one grouped finite-difference
Jacobian at the initial point (about one to two seconds for this example) and
reports, under `identifiability`, the column-scaled singular values, rank,
condition number, the expected standard error and relative error of every free
parameter from `(J^T J)^-1` with the supplied absolute sigma, parameters with
zero sensitivity, weakly determined parameters (relative error above 100% or no
finite error), pairwise correlations at or above 0.95 and the derived kex/pB
errors. These are properties of the initial point and the sampling design, so
they can warn about poorly separated R1H/R2H before any optimizer runs; they do
not replace the fitted covariance. See section 8.4 for parallel execution.

This example fits three peaks at 25 and 100 Hz nitrogen RF, each sampled at
147 uniformly spaced offsets from 105 to 135 ppm. The spacing is approximately
24.985 Hz. All 882 observations, including sidebands, are fitted with 24 free
parameters: two shared exchange rates, one RF scale, and seven parameters per
peak (five nitrogen parameters plus R1H/R2H). **No H-rate inputs are required.**
“Two RF” refers to two amplitudes at the same 1.2 GHz proton field, not two
magnetic fields.

| Quantity | Archived fit | Synthetic generating value |
|---|---:|---:|
| `kex` (s⁻¹) | 299.68025 | 300 |
| `pB` (fraction) | 0.0498302 | 0.05 |
| `v1n_scale` | 1.0822687 | 1.08 |
| Reduced χ² | 0.93569 | — |
| Jacobian rank | 24 of 24 | — |

Small numerical differences between environments are possible. Weak H-rate
identifiability and boundary warnings are expected in this noisy example;
full Jacobian rank does not imply precise H-rate estimates. Reproducing this
example checks execution, not an experimental sample.

The output prefix is `results/auto_H_two_RF`; see section 8 for its files.
To rerun, change `Project Name` in a copy of the JSON to a new prefix. If the copy
is moved to another folder, also rebase dataset/waveform paths or copy those inputs
with it; these paths resolve from the JSON directory. Relative output prefixes
resolve from the working directory. The dummy guide's step 2 handles input copies
and a fresh absolute output prefix for the bundled example. Existing outputs are
protected. For 25/50/100 Hz, run
`example/sideband_auto_H/three_RF.json` instead (1323 observations, 24 parameters).

The older single-peak, fixed-H example remains in
[full.json](manuscript/sideband_30ppm/results/1200/full.json). It has 441
observations and eight fitted parameters. Archived JSON files under `results/`
may contain original-machine absolute paths; use the portable examples above
for a fresh checkout.

## 2. Input data and units

Use one text file per combination of field, saturation time, and nominal ¹⁵N RF.
A file can contain several residue blocks. The beginning of
[full_25.txt](manuscript/sideband_30ppm/results/1200/full_25.txt) is:

```text
121.5949416
0.4
25 0
# offset intensity error
# A1 R2a: 13 R2b: 18 dw: 2.7
105 0.54861876445228 0.001
105.20547945205 0.5469790572768 0.001
105.41095890411 0.54623604725812 0.001
```

This is an excerpt, not a complete fitting dataset.

| Entry | Meaning and unit |
|---|---|
| Line 1 | Positive **¹⁵N Larmor frequency in MHz**, not the ¹H instrument frequency |
| Line 2 | Saturation duration `T`, in seconds |
| Line 3, first value | Nominal ¹⁵N RF amplitude `v1`, in Hz |
| Line 3, second value | RF inhomogeneity width `v1err`, in Hz; zero disables averaging |
| Line 4 | Required column-header line; the parser skips this line |
| Residue header | Exact residue label and initial `R2a`, `R2b` (s⁻¹), `dw` (ppm) |
| Data column 1 | **Absolute ¹⁵N chemical-shift offset in ppm** |
| Data column 2 | Normalized intensity `I/I0` |
| Data column 3 | Positive absolute standard error in the same intensity units |

Thus `0.001` means an uncertainty of 0.001 in `I/I0`, not 0.001 percent.
Field, `T`, `v1`, and intensity errors must be positive; `v1err` must be
nonnegative. All numerical input must be finite.

For another residue, append a header such as `# G2 R2a: 12 R2b: 15 dw: -2`
followed by its rows. Use each residue once per file. Labels such as `A1` or
`G23` must match the JSON exactly. Use the full residue-header format consistently
within a file. A short header such as `# A1` is accepted, but defaults to
`R2a=10`, `R2b=100`, and `dw=0.1`; explicit starting values are preferable.
Avoid inserting standalone comments inside a residue's data block: a line
starting with `#` ends that block. Data rows after such a comment need a new
residue header; otherwise loading fails instead of silently dropping points.
The loader also rejects empty datasets, residue headers without observations,
and a residue header or data row in place of the required fourth column-header
line. Comments between complete residue blocks and blank lines remain allowed.

The fitted model reports `Nz(A,T)/pA`, initialized with `(pA Nz, pB Nz)`.
It has no equilibrium-recovery source term. Far from saturation, the signal
approaches `exp(-R1*T)` under this model. Check that the experimental reference
and normalization match this convention; do not automatically force the remote
baseline to one.

## 3. Complete configuration

Save the following as `config.json` beside `full_25.txt`, `full_50.txt`, and
`full_100.txt`, or replace their paths with your input files. Values below are
the synthetic example settings and must be replaced with measured conditions.

```json
{
  "Project Name": "results/my_sideband_fit",
  "datasets": ["full_25.txt", "full_50.txt", "full_100.txt"],
  "residues": [{"name": "A1", "flag": "on"}],
  "init": {
    "Method": "Sideband",
    "kex": {"min": 200, "max": 400, "nsteps": 3},
    "pB": {"min": 0.03, "max": 0.07, "nsteps": 3},
    "initial": {"A1.peak_ppm": 120.05, "A1.R1": 1.2},
    "max_nfev": 500
  },
  "sideband": {
    "decoupling": {
      "h_larmor_mhz": 1200.0,
      "h_carrier_ppm": 8.5,
      "p90_s": 0.000070,
      "cycle": "RR",
      "b1_scale": 1.0,
      "J_hz": 92.0
    },
    "residues": {
      "A1": {"h_ppm_a": 6.2, "h_ppm_b": 6.5}
    },
    "v1n": {"mode": "scale", "initial": 1.0, "bounds": [0.8, 1.2]}
  }
}
```

Run from the project directory:

```bash
.venv/bin/python run.py /absolute/path/to/config.json
```

`sbfit.py` accepts the same command. Append `--no-pdf` to save only numerical
reports. JSON must have valid syntax, without comments or trailing commas.

**Path rules:** input data and waveform paths are relative to the configuration
file. A relative `Project Name` is resolved from the **working directory where
the command runs**. It is an output prefix, not a filename with an extension.
Existing output files cause the run to stop; use a new prefix or a new working
directory to preserve earlier fits.

**Residue selection:** loaded residues are active by default. A missing entry
in the top-level `residues` list does not turn that residue off. Explicitly set
unwanted residues to `"flag": "off"`, and supply ¹H shifts for every active one.

## 4. Specify the actual decoupling sequence

For this manual, `R = 90°x–240°y–90°x` and `cycle="RR"` repeats the same R block.
`RRbar` is a different phase sequence; it is not an alternative spelling of RR.

| Setting | Meaning |
|---|---|
| `h_larmor_mhz` | ¹H Larmor frequency in MHz; 1200 for the nominal 1.2 GHz example |
| `h_carrier_ppm` | ¹H decoupler carrier in ppm |
| `p90_s` | Nominal 90° duration in seconds; 70 μs is `0.000070` |
| `cycle` | `RR` for the example; also accepts `R`, `RRbar`, and `MLEV4` |
| `b1_scale` | ¹H RF amplitude multiplier, default 1; distinct from `v1n_scale` |
| `J_hz` | NH scalar coupling in Hz, default 92 |
| `R1H`, `R2H` | Optional legacy fixed rates in s⁻¹; omit both for automatic peakwise fitting |
| `h_ppm_a`, `h_ppm_b` | State A/B ¹H chemical shifts in ppm for each active residue |

The three pulse durations are `p90_s`, `(8/3)*p90_s`, and `p90_s`.
The proton RF amplitude is `b1_scale/(4*p90_s)` Hz; varying `b1_scale` leaves
these durations unchanged. For 70 μs and unit scale, the amplitude is about
3571.43 Hz and the R-block duration is 326.67 μs.

The sideband profile is calculated by propagating the coupled two-state NH
system through the actual segments. RF, offsets, J coupling, relaxation, and
exchange act simultaneously. The remaining fractional pulse at the end of `T`
is included. There is no separate fitted sideband amplitude or a manually
assigned sideband position.

Sequence settings and H chemical shifts are fixed inputs. When neither H rate
is specified in `decoupling` or its dataset overrides, `run.py` estimates both
as peakwise nuisance parameters, starting at R1H=2 and R2H=25 s⁻¹. These are
initial guesses, not measured values or randomly fixed assumptions. Rates are
nonnegative and shared across one peak's RF files and A/B states, never between
different peaks. Automatic fitting currently requires a single proton field.

No additional option is required. Advanced users can set
`"proton_relaxation": {"mode": "fit"}` or `{"mode": "fixed"}` inside `sideband`.
An explicit H rate retains legacy fixed behavior unless a mode is supplied;
combining explicit decoupling rates with `mode="fit"` is rejected to avoid
ambiguity. In fit mode, optional starts/bounds use `A1.R1H` and `A1.R2H` in
`init.initial`/`init.bounds`. Omit `init.vary` for the automatic workflow; an
explicit list still controls exactly which parameters are fitted.

The result records H estimates, their standard errors, H-rate correlations,
and boundary/weak-constraint diagnostics. These rates can be poorly separated
even when exchange/RF estimates are stable. Their covariance is included in
the reported uncertainty of the other parameters. Earlier identifiability
checks are in the [proton-relaxation report](PROTON_RELAXATION_FIT_CHECK.md).
If the minor-state ¹H shift is unknown, explicitly state any assumed value and
repeat the fit over plausible values. Errors in ¹H shifts or decoupling amplitude
can otherwise be absorbed into fitted nitrogen RF or exchange parameters.

For different conditions in different files, add `sideband.datasets` as an
ordered list of decoupling overrides. For three files, an example is:

```json
[{"h_carrier_ppm": 8.4}, {}, {"h_carrier_ppm": 8.6}]
```

The list length must equal the top-level `datasets` length. Each entry overrides
the shared `sideband.decoupling` settings only. The order must match the files.

## 5. Fit or fix nitrogen RF, ν1N

Replace the entire `sideband.v1n` object with one of these choices:

| Mode | Fitted RF parameter | Meaning |
|---|---|---|
| `fixed` | None | Use each input file's nominal `v1` |
| `scale` | `v1n_scale` | One common multiplier for all files and active residues |
| `per_dataset` | `v1n[0]`, `v1n[1]`, … | One RF amplitude in Hz per file, shared by its residues |

Fixed RF:

```json
{"mode": "fixed"}
```

Common scale:

```json
{"mode": "scale", "initial": 1.0, "bounds": [0.8, 1.2]}
```

Independent values for the three 25/50/100 Hz files:

```json
{
  "mode": "per_dataset",
  "initial": [25, 50, 100],
  "bounds": [[20, 30], [40, 60], [80, 120]]
}
```

With `scale`, actual RF is `nominal v1 × v1n_scale`. In `per_dataset` mode,
indices start at zero and follow the input-file order. Positive initial values
must lie inside positive bounds. If omitted, scale bounds are 0.5–1.5; per-file
bounds are 0.5–1.5 times the nominal RF.

`v1err` describes Gaussian RF inhomogeneity sampling. During RF fitting,
`v1err/v1` stays fixed. It is **not** the uncertainty of the fitted mean RF.
The latter is reported as parameter `stderr`.

The synthetic example demonstrates that ν1N can be fitted. Experimental
identifiability still depends on sampling, SNR, the available RF levels, and
the accuracy of fixed proton inputs. Compare fixed and common-scale fits
first. Use per-file RF when independent calibration differences justify the
extra parameters and the diagnostic results support them.

## 6. Parameter names, initialization, bounds, and fixing

| Parameter name | Scope | Unit / definition |
|---|---|---|
| `kab` | All active residues and files | A → B exchange rate, s⁻¹ |
| `kba` | All active residues and files | B → A exchange rate, s⁻¹ |
| `v1n_scale` or `v1n[i]` | As selected above | Dimensionless or Hz |
| `A1.peak_ppm` | Residue A1, all its files | State A ¹⁵N chemical shift, ppm |
| `A1.dw_ppm` | Residue A1, all its files | δB − δA, in ppm; either sign is allowed |
| `A1.R1` | Residue A1, both states | ¹⁵N longitudinal rate, s⁻¹ |
| `A1.R2a`, `A1.R2b` | Residue A1, states A/B | ¹⁵N transverse rates, s⁻¹ |
| `A1.R1H`, `A1.R2H` | Residue A1, both states and its RF files; automatic mode | Internal proton relaxation rates, s⁻¹ |

The reported derived parameters are `kex = kab + kba` and `pB = kab/kex`.
`pB=0.05` means 5%. Neither `kex` nor `pB` is a valid parameter name inside
`init.initial`, `init.bounds`, or `init.vary`; use `kab` and `kba` there.

`init.kex` and `init.pB` define a grid used to choose an initial exchange-rate
pair. Their min/max values **do not constrain the final fit**. Residue starting
values come from the first spectrum and its header. `init.initial` overrides
the resulting starting vector, after the grid step. A reasonable `peak_ppm`
override is useful if a sideband could be mistaken for the main minimum.

To narrow bounds, add, for example, this object as `init.bounds`:

```json
{"kab": [0, 100], "kba": [1, 1000], "A1.dw_ppm": [1, 5]}
```

Default physical limits are `kab ≥ 0`, `kba ≥ 1e-8`, and nonnegative relaxation
rates; shifts are unbounded. User bounds must stay within these limits and any
RF bounds already set in `sideband.v1n`. Use strictly ordered lower/upper bounds,
and place every initial value inside them. Fix parameters using `vary`, not
equal lower and upper bounds.

**Omitting `init.vary` fits every parameter. Providing `initial` alone does not
fix anything.** A nonempty `vary` list fits only the named parameters and fixes
all others at their initialized values. For an RF-only recovery test on the
synthetic example, replace the whole `init` object with:

```json
{
  "Method": "Sideband",
  "kex": {"min": 300, "max": 300, "nsteps": 1},
  "pB": {"min": 0.05, "max": 0.05, "nsteps": 1},
  "initial": {
    "kab": 15,
    "kba": 285,
    "v1n_scale": 1.0,
    "A1.peak_ppm": 120.0,
    "A1.dw_ppm": 3.0,
    "A1.R1": 1.5,
    "A1.R2a": 12.0,
    "A1.R2b": 15.0
  },
  "vary": ["v1n_scale"],
  "max_nfev": 500
}
```

This is a calibration demonstration with assumed known exchange and relaxation,
not a substitute for estimating those quantities from experimental data. There
must be more observations than free parameters. Try several physically plausible
starting values in separate runs; the local optimizer does not guarantee the
global minimum.

## 7. Uniform offsets over a complete 30 ppm window

The supplied 1.2 GHz example covers 105–135 ppm, centered at 120 ppm. A 30 ppm
window has width `30 × νN` Hz, where `νN` is the ¹⁵N frequency in MHz.
Thus it spans approximately 3647.85 Hz for `νN = 121.5949416 MHz`.

The study uses a target spacing of 25 Hz while keeping both endpoints and the
center. For a center `δ0`, choose:

```text
Nintervals = 2 × ceil(15 × νN / 25)
δi = δ0 − 15 + 30 × i / Nintervals,  i = 0, …, Nintervals
actual spacing (Hz) = 30 × νN / Nintervals
```

| Nominal ¹H frequency | Example ¹⁵N frequency (MHz) | Points per RF profile | Actual spacing (Hz) |
|---|---:|---:|---:|
| 600 MHz | 60.7974708 | 75 | 24.648 |
| 800 MHz | 81.0632944 | 99 | 24.815 |
| 1200 MHz | 121.5949416 | 147 | 24.985 |

Use the spectrometer's actual ¹⁵N frequency for experimental conversion. The
25 Hz target is the study's working acquisition grid, not a universal guarantee
that narrow features are resolved. Check a denser pilot profile where needed.
A fixed 30 ppm range with roughly fixed Hz spacing contains more points at
higher field, so equal per-point scans do not mean equal total acquisition time.

The fitter uses the supplied offsets directly. It does not interpolate,
resample, enforce equal spacing, or automatically mask sidebands. Fit all valid
measured points for a full-profile analysis. To compare a mask, make separate
input files with the selected rows removed, retain the same normalization and
error estimates, and use a separate configuration/output prefix. Compare
parameter uncertainty and residual structure; raw χ² values from different
numbers of observations are not directly comparable.

`demo_sideband.py --offset-step 25` is a different demonstration: it spans
±2400 Hz, approximately 39.48 ppm at this field. Use section 1's
`example/sideband_auto_H/two_RF.json` or `three_RF.json` for the **30 ppm** test.

### 7.1 Compare SBONEST and ONEST at 600 and 800 MHz

The [600/800 MHz benchmark](results/field_comparison_600_800_20261003_02/REPORT.txt)
fits the two fields **separately**. At each field, both models use the same noisy
synthetic observations, absolute σ = 0.001, parameter bounds, and optimizer.
Each has eight free parameters: `kab`, `kba`, `v1n_scale`, and the peak's
`peak_ppm`, `dw_ppm`, `R1`, `R2a`, and `R2b`.

The single peak is at 120 ppm with ΔδN = 3 ppm, T = 0.4 s, and nominal nitrogen
RF = 25/50/100 Hz. The 105–135 ppm window contains 75 points per RF at 600 MHz
(225 total) and 99 at 800 MHz (297 total). SBONEST fixes R1H = 2 and R2H = 25 s⁻¹
to the manuscript generating values. This is a **fixed-H benchmark**; section 1
describes the current automatic peakwise H-rate workflow.

From the repository root, after installing dependencies and setting the thread
variables in section 1:

```bash
.venv/bin/python compare_field_models.py \
  --source results/field_comparison_600_800_20261003_02/inputs \
  --out results/field_comparison_repeat_01
```

Always choose a new `--out` directory. Explicit `--source` is required for a fresh
checkout because the script's default refers to historical local manuscript
inputs. The bundled input copies are portable and preserve the original data.

| ¹H field | Model | kex ± SE (s⁻¹) | pB ± SE (%) | Reduced χ² |
|---|---|---:|---:|---:|
| 600 MHz | SBONEST | 297.89 ± 1.63 | 5.0287 ± 0.0299 | 1.058 |
| 600 MHz | ONEST Matrix | 297.95 ± 1.64 | 5.0142 ± 0.0357 | 1.055 |
| 800 MHz | SBONEST | 299.85 ± 1.23 | 4.9994 ± 0.0313 | 1.035 |
| 800 MHz | ONEST Matrix | 301.41 ± 1.23 | 4.9702 ± 0.0330 | 1.058 |

Truth: kex = 300 s⁻¹, pB = 5%, RF scale = 1.08. Errors are local 1 SE calculated
with absolute intensity errors and parameter cross-covariances, without
reduced-χ² rescaling. Population errors in this table are percentage points;
`parameters.csv` stores `pB` and its SE as fractions.

ONEST uses the actual Matrix worker with its zero clamp, while RF-scale handling
and optimization match SBONEST. This is not a run of the unmodified ONEST CLI.
The saved `N_signed_control` removes the clamp to isolate its contribution.
Three starts are retained for each model and noise condition (36 fits total).
The minimum χ² is reported; noisy ONEST fits reach slightly different minima,
so multistart agreement and global optimality must not be assumed.

At 600 MHz the noisy fits are nearly identical. At 800 MHz, noiseless ONEST gives
kex = 301.535 s⁻¹ (+0.51%); the signed control gives 301.551 s⁻¹ (+0.52%).
SBONEST recovers 300 s⁻¹ at both fields. The prescribed |offset from 120 ppm| =
1250–1850 Hz mask removes no acquired points, so duplicate masked fits are omitted.
These findings concern one parameter set and one noise realization per field.
ONEST's local SE excludes model mismatch; different peak positions, windows,
pulse conditions, and experimental data require their own validation.

Inspect `summary.json` for all starts, covariance matrices, diagnostics, and
source hashes; `parameters.csv` for every fitted parameter; `predictions_*.npz`
for plotted data; and `verification.json` for saved-result checks. The
[600 MHz](results/field_comparison_600_800_20261003_02/comparison_600.png) and
[800 MHz](results/field_comparison_600_800_20261003_02/comparison_800.png) figures
show profiles and standardized residuals. New runs generate fits, tables, figures,
and input copies; `REPORT.txt` and `verification.json` are curated records of the
published run and are not generated by the comparison script.

**Joint two-field fits.** `compare_joint_fields.py --source INPUTS --out NEW_DIR
[--workers N]` fits the six "full" datasets of both fields together with the
field-group models of section 12.1 and compares them with the separate fits
above; the archived run is
[`results/field_comparison_joint_20261004/summary.json`](results/field_comparison_joint_20261004/summary.json).
With the generating proton rates fixed (the benchmark assumption) the joint fit
gives kex = 298.9 ± 0.9 s⁻¹ (truth 300) against 297.9 ± 1.6 at 600 MHz alone and
299.9 ± 1.2 at 800 MHz alone. Fitting field-group proton rates leaves kex and
pB unchanged (298.9 ± 1.0, 0.0501) but the proton rates themselves are not
determined at these fields: without upper bounds they drift to the lower bound
(R1H) or to tens of thousands of s⁻¹ (R2H, together with biased nitrogen rates
when nitrogen relaxation is also field-specific), so the bounded variants use
`init.bounds` of 50 and 500 s⁻¹ (`--proton-bounds`). AICc prefers the fixed-rate
description (Δ AICc +7.8 for field-group proton rates, +12.5 with field-group
nitrogen rates), as expected when the extra parameters carry no information.
This is one noise draw per field and a synthetic truth; it shows the cost of the
extra freedom, not its benefit on experimental data.

## 8. Read and inspect the outputs

For `Project Name = results/my_sideband_fit`, the files are:

| File | Contents |
|---|---|
| `results/my_sideband_fit_predictions.csv` | Full-precision observations, predictions, sigma, standardized residuals and dataset/residue identifiers |
| `results/my_sideband_fit_result.json` | Parameters, full covariance, derived standard errors, diagnostics, configuration and provenance |
| `results/my_sideband_fit_result.txt` | Human-readable fit report and observed/calculated intensity tables |
| `results/my_sideband_fit.pdf` | Profile and standardized-residual panels, one page per active residue |
| `results/my_sideband_fit_data.pdf` | Data-only profile plots |
| `results/my_sideband_fit_checkpoint/` | Run identity and completed fitting/analysis records (section 8.3) |

The PDF covers the full supplied offset range. Calculated lines join values
at the measured offsets; they are not an independently dense simulation and
can miss narrow structure between points. Fit-PDF legend RF values are **fitted**;
the data-only PDF uses nominal RF. JSON `v1n_hz` also stores fitted values. Text
intensity tables remain rounded to three decimals; use the full-precision CSV
and original input for quantitative reanalysis. CSV `residual_sigma` is
`(observed-predicted)/sigma`; its sum of squares equals chi2.
`plot_full_profile.py` expects the separate demo's files and `validation.json`;
it is not a general plotter for arbitrary fit output.

After running section 1, inspect the result from the project directory:

```bash
.venv/bin/python - <<'PY'
import json
from pathlib import Path
r = json.loads(Path("results/auto_H_two_RF_result.json").read_text())
assert r["success"], r.get("message")
print("kex (s^-1):", r["kex"], "SE:", r["derived_se"]["kex"])
print("pB (fraction):", r["pB"], "SE:", r["derived_se"]["pB"])
print("actual RF (Hz):", r["v1n_hz"])
print("reduced chi2:", r["chi2"] / r["dof"])
print("rank:", r["jacobian_rank"], "of", r["n_parameters"])
print("warnings:", r["warnings"])
PY
```

Check these before interpreting fitted exchange parameters:

1. Inspect the entire profile, including sidebands. Look for systematic
   deviations around the main dip, minor dip, and remote baseline.
2. Verify `n_points`, the active residues, the RF-file order, and `v1n_hz`.
3. Read `warnings`, `at_bounds`, `jacobian_rank`, and `scaled_condition`.
   A rank smaller than `n_parameters` means separate parameters are not locally
   identifiable. RF correlations above 0.95 in absolute value and scaled
   condition numbers above 10⁶ trigger warnings.
4. Inspect `v1n_correlations` and repeat plausible initial conditions. A solver
   convergence message alone does not establish a unique physical result.

`stderr` is a local linear covariance estimate using the supplied intensity
errors as absolute σ. It is **not rescaled by reduced χ²**. Rank-deficient fits
report `null` errors for free parameters. A fixed parameter has `vary=false`
and `stderr=0`; this expresses an assumption, not experimental precision.
These errors exclude uncertainty in fixed proton inputs and model mismatch.
The JSON exports `covariance` with rows and columns ordered by `parameter_order`.
`derived_se.kex` and `derived_se.pB` include the covariance between `kab` and `kba`. Unavailable
entries/errors are `null`. `schema_version=2` retains existing scalar fields.
The order includes fixed parameters as well as free ones; `n_parameters` counts
only free parameters. `pB` and `derived_se.pB` are fractions: multiply both by 100
to report population in percent and its SE in percentage points. Older archived
results may lack the new fields; preserve them and run a fresh fit if needed.
`provenance` snapshots input, waveform and executed-source SHA-256 hashes before
fitting, a canonical config hash, Python/package versions, platform, thread
settings, Git state and elapsed seconds. Hashes do not replace preserved inputs.
`constraints-sideband.txt` records the tested Python 3.12 dependency combination.
Optional profile likelihood and parametric bootstrap are described in sections
8.1–8.2; neither removes model or fixed-input uncertainty.

## 8.1. Restarts, profile likelihood and performance checks

These optional settings belong inside `init`; omitting them retains a single fit.
The following is a fragment to merge into the existing `init`, not a standalone
configuration. [Dummy-guide step 6](DUMMY_GUIDE.md#6-optional-restarts-a-kex-scan-and-bootstrap)
creates and runs a complete configuration without hand-editing JSON.

```json
"multistart": {
  "starts": [{"kab": 20.0, "kba": 380.0, "v1n_scale": 1.0}],
  "random_starts": 2,
  "seed": 20261004
},
"profile": {
  "kex": [295.0, 300.0, 305.0],
  "pB": [0.045, 0.05, 0.055],
  "v1n_scale": [1.06, 1.08, 1.10]
}
```

`multistart` always includes the configured initial fit. `starts` contains
explicit parameter overrides in the names from section 6. An explicit seed is
required for random starts. Only varying parameters may be restarted; fixed
parameters remain fixed. Random starts sample uniformly inside two finite bounds;
one-sided/unbounded parameters use local perturbations described in `sb_analysis.py`.
This explores starting-point sensitivity and does not prove a global optimum.
Every attempt's full start, fitted vector, status, message and chi2 is retained in
`multistart`; `selected=true` marks the lowest chi2 among converged attempts.
The example fragment requests four attempts: one initial, one explicit and two
random. They can have different convergence statuses; inspect every record.
If all attempts fail, the CLI exits nonzero and preserves their records and
provenance in `_result.json`, with `success=false`; choose a new output prefix
before retrying. Invalid settings are rejected rather than treated as fit attempts.

`profile` requests explicit grid points for kex, pB or v1n_scale (scale mode only).
At each point all remaining varying nuisance parameters are refitted, respecting
physical bounds, user bounds and fixed parameters. kex/pB constraints transform
kab/kba exactly; they do not add a penalty residual. Infeasible or failed points
retain their target and failure message with null chi2. `profiles` stores the
baseline chi2, fitted vectors, nuisance names and raw delta chi2. Negative delta
chi2 is retained and flagged: the scan found a better solution than the baseline.
The scan does not silently replace the reported fit. Each point starts from the
selected baseline, so successful local convergence does not prove the constrained
global minimum. No confidence interval is automatically inferred from a grid:
boundaries, weak identifiability, grid extent and model adequacy require assessment.
These are within-model absolute-sigma likelihood diagnostics, not experimental
validation or uncertainty in fixed inputs. See section 8.2 for parametric bootstrap.

The top-level `success` describes the selected fit, not every requested scan
point. Read each `profiles.<name>` entry's `success` and `message`, and the
profile warnings, before interpreting the scan. `parameter_order` indexes
multistart vectors; `profiles.parameter_names` indexes profile vectors.

Grouped finite differences apply to the baseline and constrained profile refits,
including arbitrary `init.vary` subsets/orders and feasible inward steps at bounds.
Only residue-independent coordinates are grouped; exchange constraints retain
their coupled effect. The benchmark now accepts both Sideband and
ONEST configs and times a single fit (optional analysis settings are not run):

```bash
.venv/bin/python benchmark.py example/sideband_auto_H/two_RF.json
mkdir -p session_artifacts
.venv/bin/python benchmark.py example/sideband_auto_H/two_RF.json profile \
  --profile-output session_artifacts/sideband_01.prof
```

Create the parent directory first. Existing profile files are protected. The
historical `config.json profile` form still writes `benchmark_profile.prof` when
that path is unused. Profiles can be inspected with Python's `pstats` module.
Benchmarking does not write the usual fit JSON/CSV/PDF files or overwrite them.
The reported profiled timing includes profiler/reporting overhead and is not
directly comparable with an uninstrumented fit.

## 8.2. Parametric bootstrap and synthetic uncertainty checks

Add this fragment inside `init` in a new configuration with a fresh output prefix:

```json
"bootstrap": {"replicates": 5, "seed": 20261004, "confidence": 0.95}
```

`replicates` must be a positive integer, `seed` an explicit nonnegative integer,
and `confidence` between 0 and 1 (default 0.95). For each replicate, SBONEST draws
independent Gaussian noise around the selected baseline predictions using the
supplied absolute sigma. Offsets, fields, RF and fixed inputs stay unchanged.
Each synthetic dataset is fitted from the selected baseline; optional restarts
and profile scans are not repeated inside each replicate. Original observations
and the selected fit are restored afterward.

`bootstrap.samples` preserves each index, status, message, fitted vector, local
SE, χ² and boundary flags, including failed fits. `parameter_names` gives the
vector order. `bootstrap.intervals` reports lower/median/upper percentiles over
successful fits, with `n_success` and `fixed` flags for each model parameter and
derived kex/pB. Fixed intervals express assumptions. `successful`, `replicates`
and `warnings` expose failures; failures can bias the retained distribution.
Fewer than 100 successful replicates trigger a warning because tail estimates are
unstable; five are only a workflow demonstration. Even a larger count does not
guarantee nominal confidence-interval coverage. Intervals are conditional on the
chosen model, supplied absolute sigma and fixed inputs, and exclude model mismatch.

For a separate repeated-data study, use `.venv/bin/python validate_uncertainty.py --config PATH
--truth PATH --replicates N --seed N --confidence 0.95 --out NEW_DIRECTORY` with
your synthetic configuration and known generating truth. The truth JSON must be
a map from **every** model parameter name to a finite value, or `{"truth": {...}}`,
with no extra names. Values must satisfy configured bounds, and fixed parameters
must match their configured initial values. Use `--check` to inspect parameter
names; a fitted estimate is not independent known truth.

The study draws fresh Gaussian observations at that truth and fits from the truth
for each replicate. `coverage.json` measures **local normal-SE interval coverage**,
not bootstrap-percentile interval coverage. For each varying quantity it reports
the eligible count, coverage fraction, binomial SE, Wilson 95% interval, unavailable
or zero SE counts, and boundary frequency; failures are reported separately.
Coverage denominators exclude failed fits and nonpositive/unavailable SEs, so
inspect those exclusions. Fixed quantities have no coverage claim. The output
also preserves `study.json`, config/truth copies, `samples.json` and individual
`samples/` records. Small studies have large sampling uncertainty and establish
neither experimental validity nor universal coverage.

`--profile-interval kex pB` and `--inner-bootstrap B` extend the study: every
replicate then also receives likelihood-ratio intervals (section 8.5) and a
B-replicate parametric bootstrap on its own data (seed `seed·1000003 + index`),
and `coverage.json` gains an `intervals` block with the coverage fraction,
binomial SE, Wilson interval and mean width of the profile and bootstrap
intervals next to the local-SE coverage. Open or one-sided profile intervals and
failed inner bootstraps are counted separately and excluded from the fraction.
Each replicate costs one fit plus about six refits per interval quantity plus B
fits, so use `--workers N` (replicates run in parallel and the serial result is
identical). Fewer than 100 replicates remain a demonstration.

## 8.3. Checkpoints, resume and saved-result reports

A normal fit creates `PROJECT_checkpoint/` by default. It preserves an immutable
`manifest.json` identity, `provenance.json`, and a completed `baseline.json`
snapshot. Optional analyses add `attempt-0.json`, `profile-kex-0.json`,
`bootstrap-0.json` and subsequent numbered records. Failed fits retain diagnostics
in `failure.json` and the failure result JSON when available. Export staging
folders (`export-*`), `exports.json` and `complete.json` track output publication.
Files appear as their stages complete; keep the entire directory with the results.
Checkpoint records are checksum-validated. Do not edit them; changed records are
rejected on resume.

After an interruption of the section 1 example:

```bash
.venv/bin/python run.py example/sideband_auto_H/two_RF.json --resume
```

Resume skips completed baseline, restart, profile and bootstrap work. Work left
unfinished by an interruption (including Ctrl-C) runs again. A fit that already
finished with failure, including all-failed multistart, restores its failure JSON
and exits nonzero without retrying. Inspect the failure, correct the settings and
use a fresh `Project Name` for another attempt. Resume requires matching
configuration, resolved input/waveform
paths and hashes, executed-source hashes, Python/package versions, platform,
thread settings and `--no-pdf` mode. Add `--no-pdf` if it was used initially.
`--check` and `--resume` cannot be combined; a check of an existing run reports
its protected paths as conflicts. Changed settings, source or environment require
a fresh prefix. Resume does not overwrite unrelated or changed output files.
Ordinary reruns also reject an existing checkpoint, including one from a failed
run. Preserve failure records before changing settings and starting a new run.

Regenerate a report from an existing successful result and its predictions CSV:

```bash
.venv/bin/python sb_workflow.py report results/auto_H_two_RF_result.json \
  --out results/auto_H_two_RF_report_01
```

This writes `_summary.json`, `_summary.txt` and `.pdf` using only saved JSON/CSV,
without optimization or the original data files. It checks their consistency,
including point counts and χ², and adds residual statistics by residue/dataset,
fit diagnostics and any saved restarts, profiles or bootstrap records. Failed
profile points and negative Δχ² remain visible. With renamed files, use
`--predictions PATH` to select the CSV explicitly. A failure-only JSON has no
predictions to report. The command preserves existing files; choose a new report
prefix for each export. This also creates a PDF after an original `--no-pdf` fit.

## 8.4. Parallel execution with `--workers`

`run.py CONFIG --workers N` (also accepted by `sbfit.py` and by `--check
--identifiability`) starts N worker processes, each holding one copy of the model
with the same configuration and data. The main process keeps the optimizer and
the checkpoint; workers evaluate the grouped Jacobian columns of every fit, the
explicit and random restarts, the profile points and the bootstrap replicates.
Rows are collected in index order, so checkpoint records remain an ordered
prefix and `--resume` behaves exactly as in a serial run. The worker count is
not part of the checkpoint identity because results do not depend on it: the
parallel run reproduces every number of the serial run, which the regression
`test_sb_parallel.py` verifies on result JSON and checkpoint records.

Keep `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS` and `VECLIB_MAXIMUM_THREADS` at
one; the pool supplies the parallelism. On a 10-core laptop the bundled
882-point fit took about 31 s with one worker and 10 s with eight. Every worker
is started when the pool opens, the Jacobian columns are evaluated one vector
per worker, and the optimizer's own residual evaluations are split into
residue/dataset blocks across the workers; restarts, profile points and
bootstrap replicates scale almost linearly with the worker count up to the
number of items. Workers are started with the spawn method, so Python scripts
that call `run_config(..., workers=N)` must guard their entry point with
`if __name__ == "__main__":`. Progress lines with an elapsed time and a
remaining-time estimate are printed for restarts, profile points and bootstrap
replicates. After Ctrl-C, workers finish their current item before exiting;
completed items are already in the checkpoint.

## 8.5. Likelihood-ratio intervals from the profile

`init.profile` (section 8.1) evaluates a fixed grid. `init.profile_interval`
instead searches for the two points where the exact nuisance-refit profile
crosses a chi-square threshold, which gives a within-model confidence interval
without assuming a quadratic likelihood:

```json
"profile_interval": {
  "parameters": ["kex", "pB"],
  "confidence": 0.95,
  "max_evaluations": 40,
  "relative_tolerance": 0.001,
  "max_doublings": 8
}
```

`parameters` accepts `kex`, `pB` and, in scale RF mode, `v1n_scale`. The threshold
is the chi-square quantile with one degree of freedom (3.84 for 95%). Starting
from the fitted value, each side is bracketed outward in doubling steps from
`z × local SE` (10% of the estimate when no finite local error exists), then
Brent's method locates the crossing to `relative_tolerance × |estimate|`. Every
evaluated profile point is retained in evaluation order, both in the result
under `profile_intervals.<name>.points` and in the checkpoint as
`profile_interval-<name>-N` records, so `--resume` replays the search without
refitting. If the profile stays below the threshold up to a parameter bound or
within `max_doublings`, the interval is reported as open on that side with a
message instead of a number; a profile point below the base chi2 marks the base
fit as not optimal and the interval as unreliable. A failed nuisance refit is
recorded in `message` and the remaining parameters still run. Each evaluation is
one constrained refit (about 10–17 s for the bundled example serially; use
`--workers` to parallelize its Jacobian columns). The report lists the interval,
all evaluated points and a plot of the evaluated profile with the threshold.
These intervals are conditional on the model, the fixed inputs and the supplied
absolute sigma; they do not guarantee a global optimum.

## 8.6. Experimental design: expected errors before measuring

`python sb_workflow.py design DESIGN_JSON --out NEW_DIRECTORY [--workers N]`
evaluates how well a planned acquisition would determine the parameters, using
only the model and the local Fisher information at a known truth. No optimizer
runs. A design file names a base configuration (its decoupling, residues, RF
and proton settings, bounds and `vary` are reused), a truth and scenarios:

```json
{
  "config": "fit.json",
  "truth_result": "results/auto_H_two_RF_result.json",
  "scenarios": [
    {"name": "two_rf_147", "datasets": [
      {"v1n_hz": 25, "T": 0.4, "sigma": 0.01, "offsets_ppm": {"min": 105, "max": 135, "n": 147}},
      {"v1n_hz": 100, "T": 0.4, "sigma": 0.01, "offsets_ppm": {"min": 105, "max": 135, "n": 147}}]},
    {"name": "three_rf_60", "datasets": [
      {"v1n_hz": 25, "T": 0.4, "sigma": 0.01, "offsets_rel_ppm": {"min": -15, "max": 15, "n": 60}},
      {"v1n_hz": 50, "T": 0.4, "sigma": 0.01, "offsets_rel_ppm": {"min": -15, "max": 15, "n": 60}},
      {"v1n_hz": 100, "T": 0.4, "sigma": 0.01, "offsets_rel_ppm": {"min": -15, "max": 15, "n": 60}}]}
  ]
}
```

`truth` maps every model parameter name to a value, or `truth_result` points to
a saved result JSON whose fitted values are used as the truth (a fitted estimate
is a planning assumption, not independent knowledge). Each scenario dataset
gives the nitrogen RF amplitude, saturation time, absolute sigma and an offset
grid, either absolute (`offsets_ppm`) or relative to each residue's `peak_ppm`
(`offsets_rel_ppm`); `v1err_hz`, `field_mhz` and a `decoupling` override are
optional. For every scenario the command writes complete noise-free synthetic
inputs and `design_config.json` under `scenarios/<name>/`, evaluates the grouped
Jacobian at the truth and reports, in `design.json`, `design.txt` and
`design.pdf`, the point count, the acquisition proxy Σ(points × T), rank and
condition, the expected standard error and relative error of every free
parameter, the derived kex/pB errors, weakly determined parameters and strong
correlations. Expected errors scale exactly with sigma and are local linear
values at the truth; they rank designs under the model and do not validate a
sample. For the bundled example at sigma 0.01, two RF levels with 147 offsets
give an expected kex error of 6.7 s⁻¹, 60 offsets per level 10.8 s⁻¹ at 41% of
the saturation time, and a single 100 Hz level 116 s⁻¹.

**Optimizing where to measure.** An optional `optimize` section selects a
measurement budget from one scenario's dense candidate grids:

```json
"optimize": {"scenario": "two_rf_147", "budget": 60, "criterion": "kex", "min_per_dataset": 4}
```

The candidate scenario must use absolute `offsets_ppm` grids, because one
saturation offset is one spectrum row that yields every residue at once; the
selection unit is therefore the whole row. Starting from all candidate rows, the
row whose removal increases the criterion least is dropped repeatedly, using
exact Woodbury downdates of the Fisher information, until `budget` rows remain
and each dataset keeps at least `min_per_dataset`. `criterion` is the expected
variance of `kex` (default), `pB`, a free parameter name, or `D` for the
determinant criterion. The result adds two scenarios, `<name>_optimized_<budget>`
and `<name>_uniform_<budget>` (the same budget spread evenly), both evaluated
exactly like ordinary scenarios so the gain can be read off directly, plus the
selected offsets per dataset, the criterion path and a plot. For the bundled
example at sigma 0.01, 60 of 294 candidate rows (two RF levels) give an expected
kex error of 7.6 s⁻¹ against 17.5 s⁻¹ for 60 uniformly spaced rows and 6.7 s⁻¹ for
all 294. The selection is local to the truth and the model; it is a planning
aid, and a design optimized for one quantity can be worse for another.

## 8.7. Shared versus per-residue exchange

`python sb_workflow.py compare CONFIG --out NEW_DIRECTORY [--workers N] [--pdf]`
fits the configured model twice: once as configured with kab/kba (and the RF
scale) shared by all active residues (global), and once per residue with every
other residue switched off (individual). Per-residue `initial`, `bounds`,
`vary` and multistart starts are kept for the fitted residue only; profile,
bootstrap and interval analyses are not repeated. Each sub-fit runs through the
normal checkpointed path under `global/` and `individual/<residue>/`.
`comparison.json`, `comparison.txt` and `comparison.pdf` list chi2, parameter
counts, kex and pB with local errors for every model, the summed individual
chi2, AICc and BIC for both descriptions (Gaussian with the supplied absolute
sigma, up to a common constant), the preferred description by AICc and the
nested F-test of the shared model against the per-residue model. A failed
sub-fit is reported and the comparison is left empty. The statistics compare
descriptions under the supplied sigma and fixed inputs; a preference for shared
rates is consistent with one exchange process but does not prove it, and a
preference for individual rates can also reflect model mismatch or
miscalibrated errors.

**Two states or three.** `sb_workflow.py compare CONFIG --out DIR --models Sideband
Sideband_3st_Linear [--h-ppm-c A1=7.1 ...]` fits the same data with a two-state
configuration and the three-state models derived from it (section 12.2). The
derived configuration reuses the state-B proton shift for state C unless
`--h-ppm-c` gives one, starts `kbc`/`kcb` and `dwC_ppm` from five explicit
multistart combinations around the two-state exchange rate, and drops the
two-state analyses. `comparison.json`/`.txt` list chi², AICc, BIC, populations,
rates and boundary flags for every model and the AICc/BIC differences against
the two-state fit. No F-test is reported because the three-state model reduces
to two states only at the boundary of its parameter space. A return rate at its
lower bound (`kcb`) makes state C absorbing; the warning marks such a fit as
degenerate and its populations as meaningless. On synthetic three-state data the
three-state model is preferred by Δ AICc ≈ −26,000 with recovered rates; on
two-state data it degenerates and AICc prefers two states (+6.9).

## 8.8. Residual and sigma diagnostics

Every uncertainty statement in this manual assumes the supplied absolute sigma
and a model that describes the data. After each fit the result JSON contains
`residual_diagnostics`, computed from the standardized residuals alone:

- `reduced_chi2` with its expected spread `sqrt(2/dof)` and the implied
  `sigma_scale_estimate = sqrt(chi2/dof)`; a value more than three spreads from 1
  warns that sigma may be too small (or the model misses structure) or too large.
- Per residue, per dataset and per residue/dataset block: point counts and reduced
  chi-square; per block, a Wald–Wolfowitz runs test on the residual signs ordered by
  offset, the lag-1 autocorrelation with its `2/sqrt(n)` flag level and the largest
  standardized residual. Sign structure or autocorrelation calls for inspection
  of the model, noise and acquisition; neither establishes a cause.
  `runs_direction` is `alternating` for positive runs z, `clustered` for negative z,
  `balanced` for zero z, and `unavailable` for null z. Independently,
  `correlation_direction` is `positive`, `negative`, `zero`, or `unavailable`
  according to lag-1 correlation. Direction alone is not a significance test:
  the two-sided runs p-value and the absolute-correlation threshold are unchanged.
  Alternation and clustering can both trigger warnings. `overall_runs_test`
  is descriptive because concatenated block boundaries are artificial; inspect
  the offset-sorted per-block statistics.
- Outlier counts beyond 3 sigma against the Gaussian expectation.
- `rescaled_stderr` and, per varied parameter, `stderr_rescaled`: the local
  standard errors multiplied by `sqrt(chi2/dof)`. This is the conventional
  alternative when sigma is believed to be uniformly mis-estimated; it assumes a
  correct model and is not the default.

The text report and `sb_workflow.py report` print the same diagnostics (the
report recomputes them from the predictions CSV), and their warnings join the
result warnings. The tests are indicators: a reduced chi-square far from 1 does
not say whether sigma or the model is wrong, and the runs test loses power for
short blocks.

## 9. Experimental fitting workflow and limits

1. Record field frequencies, carrier, pulse phases and timing, calibrated ¹H
   amplitude, saturation time, nominal ¹⁵N RF levels, reference normalization,
   residue shifts, and intensity errors.
2. Start with an isolated residue and several RF levels. Inspect the complete
   profile before choosing a fixed-RF or common-scale fit.
3. Compare fixed and fitted RF using distinct output prefixes. Examine
   residuals, parameter shifts, rank, correlations, and calibration agreement.
4. Repeat fits over plausible fixed ¹H shifts, amplitude, and coupling. In
   automatic mode, try different peakwise H-rate starts and inspect proton
   diagnostics; if you deliberately select fixed-H mode, also vary the fixed
   H rates. Report these sensitivities separately from local `stderr`.
5. Combine residues only when a common `kab` and `kba` is scientifically
   justified. Preserve input data, configurations, and all compared results.

For multiple datasets, this implementation shares each residue's nitrogen
shifts across every file and, by default, `R1`, `R2a`, `R2b` as well. Section 12
describes field groups, field-specific proton and nitrogen relaxation and the
three-state models; compare shared and per-field fits before combining fields.

The model includes one N and one H spin in each exchanging state (two by
default; three with the models of section 12.2),
with phenomenological product-operator relaxation. State exchange can change
both N and H chemical shifts. Additional protons, a separate proton/water
exchange process, CSA–DD cross-correlation, and time-dependent RF drift are
not included. An exact segment propagator is not an exact model of every
experimental effect.

The current 1.2 GHz results support a synthetic proof of principle for using
sideband information. They do not establish that experimental fitting is
always better at 1.2 GHz. Sampling, SNR, acquisition time, and the adequacy of
the decoupling model must be tested with the planned experiment.

## 10. Troubleshooting

| Symptom | Check or action |
|---|---|
| `No module named optimalcontrol` | Use `.venv/bin/python`; install `requirements-sideband.txt` or the local OC checkout |
| `Output already exists` | Use a fresh `Project Name`, or `--resume` for the unchanged interrupted run |
| Checkpoint identity mismatch | Restore the original configuration/environment to resume, or preserve it and use a new prefix |
| File not found | Resolve dataset paths from the JSON directory; output paths use the working directory |
| Comparison cannot find `results/600/full.json` or `results/800/full.json` | Use the explicit bundled `--source` path in section 7.1 |
| Missing residue or missing ¹H shift | Match labels exactly and explicitly turn unwanted residues off |
| `Missing column-header line` | Keep the fourth column-header line; place the first residue header after it |
| `outside a residue block` | Check residue headers and remove standalone comments interrupting data rows |
| `No residue data found` or `No data points for residue` | Supply observations for each residue block and at least one block per file |
| Invalid `v1n` initial/bounds | Match mode and array length; clear scale settings when switching to `fixed` |
| Parameter name error | Use names in section 6; use `kab`/`kba`, not `kex`/`pB`, in `vary` or bounds |
| Initial values outside bounds | Check header `dw`/R2 values and explicit overrides against all limits |
| Maximum evaluations exceeded | Inspect units, starts, and identifiability first; then adjust `init.max_nfev` if justified |
| Rank warning or `stderr: null` | Examine correlations and sampling; reduce unjustified free parameters or add informative data |
| `No multistart attempt converged` | Inspect statuses/messages in the failure result JSON, then retry with a fresh prefix and justified starts/bounds |
| Negative profile Δχ² | The scan improved on the baseline; examine its vector and refit from that solution with a fresh output prefix |
| Missing `derived_se`, `covariance` or `provenance` | Check whether this is an older archived result; preserve it and produce a new fit with current code |
| `Profile already exists` | Choose a fresh `--profile-output` filename; do not remove earlier measurements |
| Long runtime | Keep numerical-library threads at one; use `--no-pdf` and `--workers N` (section 8.4); nonzero `v1err` adds RF averaging work |

The numerical regression checks can be rerun from the project directory:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  MPLBACKEND=Agg .venv/bin/python test_sideband.py
```

See [VALIDATION.md](VALIDATION.md) for the scope of those checks.

## 11. Optional OC waveform input

For a measured or designed waveform, the existing interface accepts OC
`optimalcontrol.io.export_json` output. Replace the entire decoupling object,
for example:

```json
{
  "h_larmor_mhz": 1200.0,
  "h_carrier_ppm": 8.5,
  "waveform_json": "my_oc_period.json",
  "rf_hz": 5000.0,
  "b1_scale": 1.0,
  "J_hz": 92.0
}
```

The file must contain the **complete repeated period**, including its intended
supercycle, exactly two channels `x`/`y`, uniform times starting at zero, and
positive `metadata.pulse_dt`. For `units="a.u."`, `rf_hz` supplies the amplitude
conversion; for `units="Hz"`, omit `rf_hz`. Do not combine `waveform_json` with
`p90_s` or `cycle`, including inherited per-file settings. More distinct waveform
segments increase computation time. `cestdec.Scheme` JSON is a different format
and cannot be passed directly as an OC waveform.

## 12. Several proton fields and three-state exchange

### 12.1. Field groups and field-specific relaxation

Datasets whose decoupling entries share one `h_larmor_mhz` form a **field
group**; groups are numbered in ascending field order and listed in the result
(`field_groups`) and the text report. With one group nothing changes. With
several groups:

- Automatic proton relaxation (`proton_relaxation.mode = "fit"`) fits one
  `R1H`/`R2H` pair per residue **and per group**, named `A1.R1H[0]`,
  `A1.R2H[0]`, `A1.R1H[1]`, … Fitting several fields jointly with automatic
  proton rates is therefore allowed; the earlier single-field restriction is gone.
- `"nitrogen_relaxation": {"mode": "per_field"}` in the `sideband` section gives
  each group its own `R1`, `R2a`, `R2b` (and `R2c`) per residue, named
  `A1.R1[0]`, `A1.R2a[1]`, … The default `"shared"` keeps one set for all fields,
  as before. Peak positions, chemical-shift differences, exchange rates and the
  RF parameters are always shared.
- In `init.initial`, `init.bounds`, `init.vary` and multistart starts an
  ungrouped name such as `A1.R1` addresses every group of that key; grouped
  names may still be given explicitly and take precedence.

```json
"sideband": {
  "decoupling": {"h_carrier_ppm": 8.5, "p90_s": 7e-05},
  "datasets": [{"h_larmor_mhz": 600}, {"h_larmor_mhz": 600},
               {"h_larmor_mhz": 800}, {"h_larmor_mhz": 800}],
  "proton_relaxation": {"mode": "fit"},
  "nitrogen_relaxation": {"mode": "per_field"},
  "residues": {"A1": {"h_ppm_a": 6.2, "h_ppm_b": 6.5}},
  "v1n": {"mode": "scale"}
}
```

Field-specific relaxation adds parameters; check their identifiability with
`--check --identifiability` and compare shared and per-field fits with distinct
output prefixes before drawing conclusions. Section 7.1 reports the joint
600/800 MHz benchmark with these options: kex and pB are unchanged, the proton
rates are not determined at those fields and need bounds, and AICc prefers the
fixed-rate description.

### 12.2. Three exchanging states

`init.Method = "Sideband_3st_Linear"` (A ⇄ B ⇄ C) or `"Sideband_3st_Triangle"`
(additionally A ⇄ C) propagates a 48-dimensional NH Liouvillian with the same
decoupling, relaxation and RF treatment as the two-state model. Rates are
`kab, kba, kbc, kcb` (plus `kca, kac` for the triangle); each residue gains
`dwC_ppm` (shift of state C relative to A) and `R2c`, and its
`sideband.residues` entry needs `h_ppm_c`. Populations follow from the
stationary distribution of the rate network and are reported under
`exchange.populations`; `kex` and `pB` are null because they summarize the
two-state model only, and `exchange.kex_AB`/`kex_BC` give the pairwise sums.
There is no exchange-rate grid search for three states: supply starting rates in
`init.initial`. Profiles and profile intervals of `kex`/`pB` are rejected for
three-state models; bootstrap percentiles are reported for every rate instead.
A three-site minor-state fit needs more RF levels or fields than the two-state
case and remains sensitive to starting values; use `multistart`, inspect
`--check --identifiability` and treat near-zero populations as unsupported by
the data. Setting `kbc` and `kcb` so that state C is unpopulated reproduces the
two-state result; `test_sb_models.py` verifies this and the recovery of a
synthetic linear three-state truth.

```json
"init": {
  "Method": "Sideband_3st_Linear",
  "initial": {"kab": 12, "kba": 300, "kbc": 80, "kcb": 60,
              "A1.dwC_ppm": -3.5, "A1.R2c": 18}
},
"sideband": {"residues": {"A1": {"h_ppm_a": 6.2, "h_ppm_b": 6.5, "h_ppm_c": 7.1}}, ...}
```

## 13. Installing the `sbonest` command

The repository is also an installable package. Install noneditably from the
checkout in a Python 3.12 environment (`-e .` remains available for development):

```bash
python -m pip install . -c constraints-sideband.txt
sbonest --help
sbonest version
```

Wheel and source-distribution installs include the synthetic two-RF demo and
runtime diagnostics. No checkout or `example/` directory is needed afterward.
In a writable working directory, with the numerical thread settings from section 1
and a new `demo` directory name, run:

```bash
sbonest init-demo --out demo
sbonest check demo/fit.json --no-pdf
sbonest fit demo/fit.json --no-pdf --workers 1
sbonest resume demo/fit.json --no-pdf --workers 1
sbonest report demo/fit_result.json --out demo/report
```

`init-demo` copies the bundled inputs into `demo/data/`, writes `demo/fit.json`
with relative dataset paths and an absolute output prefix, and rejects an existing
destination. A second ordinary fit also rejects existing outputs; use resume only
with matching configuration, inputs, source/dependency identity, thread settings
and PDF mode (section 8.3). `version` reads source hashes without loading a dataset.

The console command groups every tool: `sbonest check CONFIG [--identifiability]`,
`sbonest fit CONFIG [--no-pdf] [--workers N]`, `sbonest resume CONFIG`,
`sbonest report RESULT_JSON --out PREFIX`, `sbonest init-demo --out DIR`,
`sbonest design DESIGN_JSON --out DIR`, `sbonest compare CONFIG --out DIR`,
`sbonest import-bruker ...` (section 14), `sbonest serve` (section 15),
`sbonest benchmark CONFIG [profile]` and `sbonest version`, which prints the
package version and the executed-source hashes recorded in provenance. Every
command calls the same functions as the scripts (`run.py`, `sb_workflow.py`, …),
so outputs, checkpoints and provenance are identical; the scripts remain usable
from a plain checkout without installation. The modules stay flat at the
repository root, which keeps the source hashes in `provenance` meaningful.

Installed `compare` accepts the script's `--models` and `--h-ppm-c` options.
With `--models Sideband Sideband_3st_Linear`, supply a two-state `Sideband`
configuration; three-state configurations are derived from it. `--h-ppm-c A1=7.1`
supplies state C's proton shift for A1 (otherwise its state-B shift is used).
Without `--models`, comparison remains shared versus per-residue exchange.
The same `--workers N` and `--pdf` controls apply. Model comparisons remain
within-model statements under the supplied absolute sigma, not experimental validation.

## 14. Importing Bruker pseudo-2D data

`python sb_import.py PDATA OFFSETS --out NEW_FILE --peak LABEL=ppm[:dw_ppm] ...`
(or `sbonest import-bruker ...`) converts a processed Bruker pseudo-2D CEST
experiment into one SBONEST dataset file. `PDATA` is the processed directory
(`.../pdata/1`) with `procs`, `proc2s` and `2rr`; the reader handles the
submatrix layout, both byte orders, int32 and float64 storage and the
`NC_proc` scaling with NumPy alone. The detected dimension is the proton axis
(ppm from `OFFSET`, `SW_p`, `SF`); each row is one saturation offset.

Required inputs are explicit: `OFFSETS` lists one saturation offset per row
(ppm on the saturated nucleus, or Hz with `--offset-unit hz --carrier-ppm X`);
`--peak A1=8.30:2.5` gives each residue's proton position and optional
starting `dw` (repeatable); `--reference-row K` names the reference spectrum
that normalizes intensities and is excluded from the output; `--noise-region
LO HI` is a signal-free proton range whose standard deviation in the reference
row becomes the absolute error; `--saturation-s` and `--v1-hz` fill the header.
`--half-width` (ppm) and `--mode max|sum` select the window statistic, `--r2a`,
`--r2b` set the header starting values, `--exclude-row` removes further rows,
and the saturated-nucleus field defaults to `SF × γ(15N)/γ(1H)` (`--nucleus`,
`--field-mhz`). The output is the text format of section 2 and a JSON summary is
printed. Phase and baseline quality, peak overlap and the choice of reference
row are the user's responsibility; inspect the converted profiles before fitting.

Three helpers make that inspection easier. Bruker frequency lists (`fq1list`,
…) are read directly: a first line such as `bf ppm`, `sfo hz` or `P` sets the
unit and `O1`/`O2` lines are skipped, so `--offset-unit` is needed only for a
plain numeric list. `--peaks-from-reference [--peak-snr 10]` replaces `--peak`
entries by the local maxima of the reference row above the threshold times the
noise (labels `P1`, `P2`, … in order of height, at least 0.03 ppm apart); edit
the labels and starting `dw` values in the written file or rerun with explicit
`--peak` entries once the assignments are known. `--qa-pdf NEW_FILE` writes a
figure with the reference row, the peak windows, the noise region and every
extracted profile with its error bars, which is the quickest way to see a wrong
reference row, an overlapping window or a baseline problem.

## 15. Web runner for Sideband fits

`python sb_server.py [--host 127.0.0.1 --port 5050]` (or `sbonest serve`) starts
a Flask page for people who prefer a browser to the terminal. A job is created
by uploading a Sideband configuration and its data files (and an optional OC
waveform); the server rewrites `datasets` to the uploaded names and the output
prefix to `fit`, stores everything under `SB_JOBS/<job id>/` (or
`SBONEST_JOBS_DIR`), runs the preflight check with identifiability and shows
its JSON. Buttons start the fit as a background `run.py` process with the
chosen worker count, resume an interrupted job through its checkpoint, and
regenerate reports; the page polls the job status (process state, checkpoint
records, result summary, log tail) and links to every output for download. Only
files inside a job directory are served.

The upload form also accepts optional analyses that are written into `init`:
random restarts with a seed, a kex profile grid, profile intervals for kex/pB
and a seeded bootstrap. After a fit the page shows a PNG preview of every
residue (`/jobs/<id>/preview.png`, rendered from the predictions CSV and
cached), and a job can be moved to `SB_JOBS/archive/` with the archive button;
`--max-age-days D` archives idle jobs older than D days when the server starts
or when `POST /jobs/archive-expired` is called. Nothing is ever deleted. Set
`SBONEST_TOKEN` (or `--token`) to require an access token on every request
except the page itself; the page asks for it and sends it as the
`X-SBONEST-Token` header (or `?token=` for downloads). Without a token the
runner is meant for a trusted local network only.

The installed command accepts the same token and archive options, for example:

```bash
sbonest serve --host 127.0.0.1 --port 5057 --token qa-token --max-age-days 30
```

Use your own access token instead of this demonstration value. Archive age must
be positive. Fit, resume, report and archive errors are displayed separately from
the job log and remain visible across status refreshes. Failed archive preserves
the selected job and its controls; stale status responses cannot restore an
archived selection. Authentication, network and invalid-response errors are
shown rather than treated as success. PNG previews work with and without a token.

## 16. Module layout and API reference

`sbfit.py` holds the model (`SidebandModel`, configuration validation);
`sb_run.py` holds the run workflow (`check_config`, `run_config`, the command
line shared by `run.py` and `sbfit.py`), and both names remain importable from
`sbfit`. `docs/API_REFERENCE.md` lists every public function and class of the
Sideband modules with its signature and first docstring line; regenerate it with
`python generate_api_reference.py` after changing a public signature or docstring (CI runs
`--check`).

Implementation references: [sbfit.py](sbfit.py), [sideband.py](sideband.py),
[run.py](run.py), [est_data.py](est_data.py), [sb_analysis.py](sb_analysis.py),
[sb_checkpoint.py](sb_checkpoint.py), [sb_workflow.py](sb_workflow.py),
[sb_report.py](sb_report.py), [sb_bootstrap.py](sb_bootstrap.py),
[sb_design.py](sb_design.py), [sb_compare.py](sb_compare.py), [sb_import.py](sb_import.py),
[sb_server.py](sb_server.py), [sb_cli.py](sb_cli.py) and [validate_uncertainty.py](validate_uncertainty.py).
