# SBONEST

Fit two-state ¹⁵N CEST profiles **including 90°x–240°y–90°x ¹H decoupling
sidebands**, using NH spin operators from OC. Fit nitrogen RF (`v1n`) as a
common scale or independently for each dataset, or keep it fixed.

**R1H and R2H need no user input.** If both are omitted, SBONEST estimates them
as independent nuisance parameters for each peak, alongside exchange, nitrogen
relaxation, and RF parameters. It starts at (2, 25) s⁻¹; these are starting
values, not arbitrary fixed rates. H rates are shared across RF datasets and
states A/B for the same peak. Several proton fields can be fitted jointly with
field-specific proton (and optionally nitrogen) relaxation, and linear or
triangular three-state exchange is available (manual section 12).

[English manual](SBONEST_MANUAL.md) · [한국어 매뉴얼](SBONEST_MANUAL.ko.md) ·
[Model and configuration](SIDEBAND.md) · [Contributor guide](AGENTS.md) ·
[Development handoff](HANDOFF.md)

**First time here?** Follow the [step-by-step dummy-data guide](DUMMY_GUIDE.md)
or [한국어 따라하기](DUMMY_GUIDE.ko.md): install, prepare synthetic inputs with
`init-demo`, check and fit, inspect saved reports, and try restarts, scans or bootstrap.

## Install and run

Python 3.12 is the CI target. From a new checkout:

```bash
git clone https://github.com/Gohyang-Matzip/sbonest.git
cd sbonest
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements-sideband.txt -c constraints-sideband.txt
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg
.venv/bin/python run.py example/sideband_auto_H/two_RF.json --check
.venv/bin/python run.py example/sideband_auto_H/two_RF.json
```

OC is installed as `optimalcontrol-nmr`; a neighboring OC checkout is not
required. Add `--no-pdf` to both check and fit to skip figures while keeping numerical
outputs. `--check` prints a JSON summary of inputs, parameters and output conflicts
without optimization or file writes; an initialization grid may evaluate the model.
`--check --identifiability` adds one grouped Jacobian at the initial point: scaled
singular values, rank, condition, expected standard errors, weakly determined
parameters and strong correlations, so weak proton rates are visible before a fit.
`--workers N` runs Jacobian columns, restarts, profile points and bootstrap
replicates in N worker processes with results identical to `--workers 1`
(on this 10-core Mac the 882-point fit takes about 33 s serially and 16 s with 8
workers; bootstrap replicates and profile points scale almost linearly). Keep
the numerical-library thread variables at one in both cases.

The first example fits three peaks at **1.2 GHz, 25/100 Hz nitrogen RF,
105–135 ppm**, with 147 equally spaced offsets per peak and RF
(about 24.985 Hz between offsets). All 882 observations, including sidebands,
are used; all 24 model parameters are fitted. No R1H/R2H inputs are present.

Outputs are `results/auto_H_two_RF_result.json`,
`results/auto_H_two_RF_result.txt`, `results/auto_H_two_RF_predictions.csv`,
`results/auto_H_two_RF.pdf`, and
`results/auto_H_two_RF_data.pdf`, plus the run-owned
`results/auto_H_two_RF_checkpoint/` directory. Existing outputs are protected: change
`Project Name` to a new output prefix before rerunning. Dataset paths are
relative to the configuration file; output prefixes are relative to the
working directory.

For a fresh practice folder, saved-result report, or interrupted run:

```bash
.venv/bin/python sb_workflow.py init-demo --out session_artifacts/dummy_01
.venv/bin/python sb_workflow.py report results/auto_H_two_RF_result.json \
  --out results/auto_H_two_RF_report_01
.venv/bin/python run.py example/sideband_auto_H/two_RF.json --resume
```

`init-demo` requires a new directory and creates `fit.json` plus `data/` input
copies. `report` uses only the saved JSON and its matching predictions CSV to
write `_summary.json`, `_summary.txt` and `.pdf` without refitting; use a new report
prefix each time. Resume reuses completed baseline/restart/profile/bootstrap
records. It requires unchanged config, inputs/waveforms, executed source,
Python/packages, platform, thread settings and `--no-pdf` mode. Preserve the
checkpoint without editing its checksum-validated records. Interrupted work can
continue, but a completed failed fit is restored and exits nonzero without retrying.
Corrected settings require a fresh `Project Name`. See manual section 8.3 for
checkpoint contents and recovery limits.

For three RF amplitudes (25/50/100 Hz), run:

```bash
.venv/bin/python run.py example/sideband_auto_H/three_RF.json
```

Here “two RF” means two saturation amplitudes at the **same magnetic field**.
It does not mean two spectrometers. All included demonstrations are synthetic;
experimental validation remains to be performed.

## Results and diagnostics

The archived automatic fits gave:

| Nominal RF (Hz) | Observations | Parameters | kex (s⁻¹) | pB | RF scale | Reduced χ² |
|---|---:|---:|---:|---:|---:|---:|
| 25 / 100 | 882 | 24 | 299.68025 | 0.0498302 | 1.0822687 | 0.93569 |
| 25 / 50 / 100 | 1323 | 24 | 299.49797 | 0.0499607 | 1.0806102 | 0.94625 |

Generating values: kex = 300 s⁻¹, pB = 0.05, RF scale = 1.08.
Three random peakwise H-rate starts per design converged to the same χ² within
2 × 10⁻⁷. Individual H rates remained weakly constrained; their uncertainty
is included in the local covariance of the other fitted parameters.

![Full 30 ppm synthetic profiles with automatic peakwise proton relaxation](results/auto_H_refit/fits/full_profile_fit.png)

Read `warnings`, `at_bounds`, `jacobian_rank`, `scaled_condition`,
`v1n_correlations`, and proton diagnostics in the result JSON. A converged fit
does not establish that every parameter is identifiable. Local standard errors
use the supplied absolute intensity errors and are not scaled by reduced χ².
The JSON includes the full `covariance` matrix ordered by `parameter_order`,
`derived_se` for kex/pB accounting for rate covariance, and calculation-time input/source hashes and environment
versions in `provenance`. The CSV preserves full precision; the fit PDF includes
standardized residuals. Optional `init.multistart` and `init.profile` provide
reproducible restarts and constrained nuisance refits. Optional `init.bootstrap`
requires `replicates` and an explicit `seed` (with `confidence`, default 0.95);
results retain every replicate and successful-sample percentile intervals. Fewer
than 100 successful replicates trigger a warning. These intervals are conditional
on the selected model, fixed inputs and supplied absolute sigma, with no guarantee
of nominal coverage. `validate_uncertainty.py` separately measures synthetic
**local normal-SE interval coverage**, including failures and Wilson sampling
intervals; it does not measure bootstrap-interval coverage. See manual sections
8.1–8.3 and [workflow sources](sb_workflow.py), [checkpoint sources](sb_checkpoint.py),
[report sources](sb_report.py), [bootstrap sources](sb_bootstrap.py) and
[coverage-study sources](validate_uncertainty.py).

Three further tools (manual sections 8.5–8.7): `init.profile_interval` locates
the likelihood-ratio interval of kex, pB or the RF scale on the exact nuisance-refit
profile (every evaluated point retained and resumable);
`sb_workflow.py design DESIGN_JSON --out DIR` reports expected standard errors,
rank and correlations of planned acquisitions at a known truth without fitting
([design sources](sb_design.py)); `sb_workflow.py compare CONFIG --out DIR` fits the
shared-exchange model and per-residue models and compares chi2, AICc, BIC and a
nested F-test ([comparison sources](sb_compare.py)). All three are within-model
statements under the supplied absolute sigma.

[Preview of the new residual report](docs/diagnostics-preview.png)
(A1 from the same synthetic two-RF example).

- [Full-profile refit report, parameters, and random-start checks](results/auto_H_refit/fits/REFIT_REPORT.md)
- [Numerical validation](VALIDATION.md)
- [Peakwise proton relaxation checks](PROTON_RELAXATION_FIT_CHECK.md)
- [Random fixed H-rate sensitivity](RANDOM_PROTON_RELAXATION_CHECK.md)
- [Two versus three RF amplitudes](TWO_RF_FIT_CHECK.md)

The last two reports describe earlier sensitivity studies. The current default
fits peakwise H rates. Explicit legacy `R1H`/`R2H` settings retain fixed-rate
behavior; see the manuals before reusing an older configuration.

## Fit experimental data

Use one text file per field, saturation time, and nominal nitrogen RF, with
residue blocks containing absolute nitrogen offsets (ppm), normalized
intensities, and positive absolute standard errors. The first file header is
the **¹⁵N Larmor frequency in MHz**, not the proton instrument frequency.
Copy a portable example configuration and replace the data, residue labels,
proton shifts, and measured pulse conditions. The manuals cover the complete
format, bounds, parameter sharing, normalization, and RF calibration choices.

The supported Sideband workflow is `run.py` (or the installed `sbonest`
command) with a Sideband `init.Method`. Bruker processed pseudo-2D data can be
converted with `sb_import.py` (explicit offsets, peaks, reference row and noise
region; manual section 14), and `sb_server.py` provides a browser runner for
Sideband jobs with background fits, resume and reports (manual section 15).
The inherited ONEST web interface (`server_run.py`), `prepare.py` and `mcrun.py`
are for ONEST models.

To install the command line into an environment:

```bash
python -m pip install -e . -c constraints-sideband.txt
sbonest version
sbonest check example/sideband_auto_H/two_RF.json --identifiability
```

`sbonest check|fit|resume|report|init-demo|design|compare|import-bruker|serve|benchmark|version`
call the same functions as the scripts, so outputs and provenance are identical
(manual section 13). The scripts keep working from a plain checkout.

Exact pulse-segment propagation is used, but the model has one N/H pair per
state and phenomenological relaxation. Automatic proton relaxation does not
remove sensitivity to incorrect proton shifts, RF calibration, or pulse timing.
Synthetic results do not establish that 1.2 GHz always improves experimental
fits; compare sampling, SNR, acquisition time, and model adequacy experimentally.

## Reproduce the checks

After setting the thread variables above:

```bash
.venv/bin/python test_performance.py
.venv/bin/python test_debugging.py
.venv/bin/python verify_3state.py
.venv/bin/python test_sideband.py
.venv/bin/python test_grouped_jacobian.py
.venv/bin/python test_benchmark.py
.venv/bin/python test_sb_analysis.py
.venv/bin/python test_sb_output.py
.venv/bin/python test_sb_check.py
.venv/bin/python test_sb_checkpoint.py
.venv/bin/python test_sb_workflow.py
.venv/bin/python test_profile_jacobian.py
.venv/bin/python test_sb_bootstrap.py
.venv/bin/python test_sb_parallel.py
.venv/bin/python test_profile_interval.py
.venv/bin/python test_sb_design.py
.venv/bin/python test_sb_compare.py
.venv/bin/python test_sb_models.py
.venv/bin/python test_sb_import.py
.venv/bin/python test_sb_server.py
.venv/bin/python test_sb_cli.py
.venv/bin/python demo_sideband.py --out session_artifacts/sideband_demo
```

Use a new output directory for each demo run. The demo covers ±2400 Hz,
about 39.48 ppm at 1.2 GHz; the `example/sideband_auto_H` configurations are
the full **30 ppm** examples. To repeat the archived ten automatic H-rate fits
(two/three RF, exact controls, and three random starts), use a new directory:

```bash
.venv/bin/python refit_automatic_proton.py --out results/auto_H_repeat_01
```

This is a longer study, not necessary for a first fit. Archived `results/`
files preserve calculation-time absolute paths and source hashes. Use the
portable configurations for new fits; do not rewrite archived provenance.
Raw experimental data, historical manuscript drafts, environments, and local
backups are excluded. Selected manuscript DOCX files and their numerical evidence
are published in [the manuscript bundle](manuscript/sideband_30ppm/README.md).
Small synthetic inputs under `manuscript/sideband_30ppm/results/1200/` support
the older manual examples; the 600/800 MHz comparison bundles its own inputs below.

## Compare SBONEST and ONEST at 600/800 MHz

The [comparison report](results/field_comparison_600_800_20261003_02/REPORT.txt)
compares both models on identical synthetic 105–135 ppm data at 25/50/100 Hz
nitrogen RF. Each model fits eight parameters, including a shared RF scale;
SBONEST fixes R1H = 2 and R2H = 25 s⁻¹, matching the original manuscript benchmark.
The two magnetic fields are fitted separately, with 225 observations at 600 MHz
and 297 at 800 MHz. These are separate benchmarks from the automatic H-rate
examples above.

| ¹H field | Model | kex ± SE (s⁻¹) | pB ± SE (%) | Reduced χ² |
|---|---|---:|---:|---:|
| 600 MHz | SBONEST | 297.89 ± 1.63 | 5.0287 ± 0.0299 | 1.058 |
| 600 MHz | ONEST Matrix | 297.95 ± 1.64 | 5.0142 ± 0.0357 | 1.055 |
| 800 MHz | SBONEST | 299.85 ± 1.23 | 4.9994 ± 0.0313 | 1.035 |
| 800 MHz | ONEST Matrix | 301.41 ± 1.23 | 4.9702 ± 0.0330 | 1.058 |

Generating values are kex = 300 s⁻¹ and pB = 5%. Errors are local 1 SE with
absolute σ = 0.001, without reduced-χ² rescaling; pB errors are percentage points.
The 600 MHz fits are nearly identical. At 800 MHz, the noiseless ONEST fit retains
a +0.51% kex displacement, while SBONEST recovers the truth at both fields.
The prescribed 1250–1850 Hz sideband mask removes no points in these windows.

The bundle includes [all parameters and local standard errors](results/field_comparison_600_800_20261003_02/parameters.csv),
[600 MHz](results/field_comparison_600_800_20261003_02/comparison_600.png) and
[800 MHz](results/field_comparison_600_800_20261003_02/comparison_800.png) figures,
three-start fits, noiseless controls, input copies, and source hashes.

For a fresh checkout, explicitly select the bundled inputs; the script's default
source points to historical local manuscript data that are not distributed:

```bash
.venv/bin/python compare_field_models.py \
  --source results/field_comparison_600_800_20261003_02/inputs \
  --out results/field_comparison_repeat_01
```

Choose a new output directory each time. ONEST uses its actual Matrix forward
model, including its zero clamp, with RF-scale handling and optimization matched
to SBONEST. This is not a comparison with default ONEST CLI settings. The lowest
χ² among three starts is reported; noisy ONEST fits reached slightly different
minima. A signed nitrogen-only control is also retained. These results describe
one synthetic noise realization per field; local errors exclude model mismatch,
and performance at other peak positions or pulse conditions remains untested.
See section 7.1 of the [English](SBONEST_MANUAL.md) or
[Korean](SBONEST_MANUAL.ko.md) manual for the full procedure.

## ONEST and license

SBONEST is based on [ONEST](https://github.com/dleess/ONEST), originally
[jhyeokchoi/ONEST](https://github.com/jhyeokchoi/ONEST). Its Git history and
[GNU GPL v3 license](LICENSE) are preserved. The inherited ONEST documentation
is available in [English](MANUAL.md) and [한국어](MANUAL.ko.md).

The inherited ONEST citation is: Choi, J., Lee, SY., Han, K., Carneiro, M. G.,
Ryu, KS., and Lee, D., “ONEST: A Web-Based Platform for the Rapid and Robust
Analysis of Protein Excited States through CEST Spectroscopy” (listed as
in preparation in the inherited documentation).
