# SBONEST fitting manual

[한국어](SBONEST_MANUAL.ko.md) · [README](README.md) · [Technical notes](SIDEBAND.md)

This manual covers the command-line `Sideband` model in this checkout, checked on
2026-10-02. It fits two-state ¹⁵N CEST profiles, including ¹H decoupling sidebands,
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

```bash
git clone https://github.com/Gohyang-Matzip/sbonest.git
cd sbonest
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-sideband.txt
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg
.venv/bin/python run.py example/sideband_auto_H/two_RF.json
```

OC is installed as `optimalcontrol-nmr`; no neighboring OC checkout is needed.
Reuse a working environment when one already exists. Add `--no-pdf` to the last
command if only numerical outputs are needed.

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
To rerun, change `Project Name` in a copy of the JSON to a new prefix. Existing
outputs are protected. For 25/50/100 Hz, run
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
starting with `#` ends that block.

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

## 8. Read and inspect the outputs

For `Project Name = results/my_sideband_fit`, the files are:

| File | Contents |
|---|---|
| `results/my_sideband_fit_result.json` | Parameters, RF in Hz, diagnostics, and input configuration |
| `results/my_sideband_fit_result.txt` | Human-readable fit report and observed/calculated intensity tables |
| `results/my_sideband_fit.pdf` | Experimental data and calculated profiles, one page per active residue |
| `results/my_sideband_fit_data.pdf` | Data-only profile plots |

The PDF covers the full supplied offset range. Calculated lines join values
at the measured offsets; they are not an independently dense simulation and
can miss narrow structure between points. Legend RF values are **nominal**;
read actual fitted RF from JSON `v1n_hz`. Text intensity tables are rounded to
three decimals; retain the original input for quantitative reanalysis.
`plot_full_profile.py` expects the separate demo's files and `validation.json`;
it is not a general plotter for arbitrary fit output.

After running section 1, inspect the result from the project directory:

```bash
.venv/bin/python - <<'PY'
import json
from pathlib import Path
r = json.loads(Path("results/auto_H_two_RF_result.json").read_text())
print("kex:", r["kex"], "pB:", r["pB"])
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
The JSON does not export the full covariance matrix; do not derive errors for
`kex` or `pB` by assuming `kab` and `kba` are independent. Profile-likelihood and
bootstrap procedures are not built-in CLI commands in this version.

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
shifts **and `R1`, `R2a`, `R2b` across every file**. It has no independent
per-field nitrogen relaxation parameters. Combining fields therefore imposes
that constraint; separate field fits may be needed to assess its effect.
The `sideband.datasets` override cannot change this parameter-sharing rule.

The model includes one N and one H spin in each of two exchanging states,
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
| `Output already exists` | Choose a new `Project Name` or run in a fresh directory |
| File not found | Resolve dataset paths from the JSON directory; output paths use the working directory |
| Missing residue or missing ¹H shift | Match labels exactly and explicitly turn unwanted residues off |
| First residue disappears | Keep the fourth header line; do not place the first residue header there |
| Unexpectedly few points | Check residue headers and comments inside data blocks; inspect JSON `n_points` |
| Invalid `v1n` initial/bounds | Match mode and array length; clear scale settings when switching to `fixed` |
| Parameter name error | Use names in section 6; use `kab`/`kba`, not `kex`/`pB`, in `vary` or bounds |
| Initial values outside bounds | Check header `dw`/R2 values and explicit overrides against all limits |
| Maximum evaluations exceeded | Inspect units, starts, and identifiability first; then adjust `init.max_nfev` if justified |
| Rank warning or `stderr: null` | Examine correlations and sampling; reduce unjustified free parameters or add informative data |
| Long runtime | Keep numerical-library threads at one; use `--no-pdf`; nonzero `v1err` adds RF averaging work |

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

Implementation references: [sbfit.py](sbfit.py), [sideband.py](sideband.py),
[run.py](run.py), and [est_data.py](est_data.py).
