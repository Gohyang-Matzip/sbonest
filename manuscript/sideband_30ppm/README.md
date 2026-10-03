# Sideband CEST manuscript and fitting uncertainties

- [Current manuscript](Sideband_CEST_1p2GHz_parameters_errors.docx)
- [Tracked uncertainty revision, author Donghan Lee](Sideband_CEST_1p2GHz_parameters_errors_tracked.docx)
- [Baseline before adding uncertainties](Sideband_CEST_1p2GHz_parameters_13C.docx)

The manuscript consistently uses **1.2 GHz** for the proton instrument frequency.
Tables 5 and 6 compare conventional nitrogen-only fitting with SBONEST on the
same sideband-containing synthetic observations, both with and without the
specified sideband mask. Table 6 includes all eight independent parameters and
the derived minor-state shift. Table 7 reports the isolated ¹³C–¹H configuration
check, including every fitted parameter and derived exchange quantities.

Tables 3, 5, 6 and 7 contain 67 estimates with **± 1 local standard error**.
The covariance uses the supplied absolute intensity error, σ = 0.001, without
rescaling by reduced χ². The errors for kex, pB and δB include parameter
cross-covariances. Population errors beside percentages are percentage points.
Generating values and fixed inputs have no fitted uncertainty.

The conventional-model errors describe within-model local curvature and exclude
model discrepancy. The noiseless carbon check uses assigned noise weights;
its standard errors describe hypothetical measurement noise, not experimental
replicate scatter. The manuscript benchmark fixes proton relaxation explicitly;
it is separate from the repository's automatic proton-relaxation examples.
All results here are synthetic, and experimental validation remains prospective.

## Evidence

| Artifact | Contents |
| --- | --- |
| `results/1200/` | Original 1.2 GHz NH data, configurations and fits |
| `results/summary.json` | Original field-study summary and stored standard errors |
| `results/conventional_comparison_20261003/` | Matched nitrogen-only fits, controls, source hashes and predictions |
| `results/carbon_configuration_20261003/` | Independent CH data, configuration and recovered parameters |
| `results/fit_errors_20261003/summary.json` | Nine fitted-point covariance matrices, propagated errors and source hashes |
| `qa/fit_errors_20261003/revision_audit.json` | Document revision ledger and source/output hashes |
| `qa/fit_errors_20261003/validation.json` | Verification of 67 table entries, tracked changes and 19-page document renders |

The four Python scripts are preserved as executed research sources, including
their formatting, so hashes in the numerical records continue to identify the
exact sources. The DOCX files contain the figures and equations. Historical
drafts, Word lock files, page-render images and local backups remain ignored.
Absolute paths in the audit record identify the original local run; the files
are available at the repository-relative paths listed above.

## Reproduction

Install the repository's `requirements-sideband.txt` in a Python 3.12 environment.
The scripts resolve inputs relative to this directory. Their output paths are
fixed and intentionally refuse to overwrite existing artifacts. In a separate
checkout, move the output being regenerated to an archival path before running:

```bash
python manuscript/sideband_30ppm/compare_conventional.py
python manuscript/sideband_30ppm/check_carbon_configuration.py
python manuscript/sideband_30ppm/collect_fit_errors.py
```

The comparison also creates a new `figures/conventional_comparison_20261003/`
directory. The covariance script reconstructs errors at saved optima, checks
the original standard errors and objectives, and does not refit the data.

For document reproduction, install `lxml`, archive both current output DOCX
files, and run `python manuscript/sideband_30ppm/revise_fit_errors.py`. It uses
the included baseline DOCX and covariance summary and checks that rejecting
the uncertainty revision restores the baseline text. Document render checks
must be repeated after regeneration; the archived QA hashes describe the
delivered files from 2026-10-03.

## Separate 600/800 MHz comparison

The subsequent [field comparison](../../results/field_comparison_600_800_20261003_02/REPORT.txt)
uses the original synthetic 30 ppm inputs, fixed proton relaxation, and eight
free parameters per model. Each field is fitted independently. Its figures and
results are a separate repository study; they have not been added to the
1.2 GHz DOCX tables listed above.

At 600 MHz, SBONEST and ONEST Matrix give nearly identical noisy fits. At
800 MHz, ONEST's noiseless kex is 301.535 s⁻¹ (+0.51% from truth); SBONEST
recovers 300 s⁻¹ at both fields. The actual ONEST worker includes its zero clamp;
a signed nitrogen-only control is retained. RF-scale handling, bounds, and
optimization match SBONEST. All three starts are preserved, including slightly
different noisy ONEST minima. These are synthetic comparisons with within-model
local standard errors, not experimental accuracy estimates.

From the repository root, use the bundled input copies and a new output path:

```bash
.venv/bin/python compare_field_models.py \
  --source results/field_comparison_600_800_20261003_02/inputs \
  --out results/field_comparison_repeat_01
```

The original local `results/600/` and `results/800/` directories are not
distributed. [The bundle](../../results/field_comparison_600_800_20261003_02)
contains their input copies, all fitted parameters, covariances, predictions,
figures, and source hashes. See the [SBONEST manual](../../SBONEST_MANUAL.md),
section 7.1, for units, uncertainty conventions, and interpretation limits.
