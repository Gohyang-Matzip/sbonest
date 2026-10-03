# HANDOFF: SBONEST public release and fitting workflow

**Updated:** 2026-10-03 · **Working directory:** `/Users/donghanlee/work/projects/sbonest`

**Repository:** https://github.com/Gohyang-Matzip/sbonest · **Published branch:** `main`

## Purpose

Maintain the SBONEST program, portable synthetic examples, numerical evidence,
and English/Korean documentation. This handoff records the published scientific
baseline through PR #3; inspect GitHub for subsequent documentation or code work.
The original release workflow is complete, not a pending task to repeat.

## Implementation status

The public repository includes the SBONEST release ([PR #1](https://github.com/Gohyang-Matzip/sbonest/pull/1)),
the manuscript with fitting uncertainties ([PR #2](https://github.com/Gohyang-Matzip/sbonest/pull/2)),
and the 600/800 MHz model comparison ([PR #3](https://github.com/Gohyang-Matzip/sbonest/pull/3)).
PR #3 merged as `8d484de94a198b9f943ec75cb053b3f1007aaf34`; its
[post-merge CI](https://github.com/Gohyang-Matzip/sbonest/actions/runs/37114856838)
passed, and local `main` was synchronized on 2026-10-03.
No experimental dataset has been fitted; supplied demonstrations are synthetic.

## What is applied

- **[still applied]** `sideband.py`: OC-based two-state NH propagation through
  the actual repeated 90°x–240°y–90°x pulse segments, including partial periods.
  `sbfit.py`: shared exchange, peakwise nitrogen shifts/relaxation, fixed/common
  scale/per-dataset nitrogen RF, and diagnostics. `run.py` routes `Sideband`.
- **[still applied]** Omitting both R1H/R2H selects automatic peakwise nuisance
  fitting, initialized at (2, 25) s⁻¹. Each peak shares its H rates across RF
  files and states A/B. Different peaks do not share H rates. The current
  automatic mode accepts one proton magnetic field. Explicit legacy H inputs
  preserve fixed-rate behavior. Users need not enter or measure H rates.
- **[still applied]** `fit.py` accepts a grouped finite-difference step used by
  the sideband model. The dense Jacobian retains covariance estimation.
  `estmodel.py` saves PDF figures with the external legend included.
- **[still applied]** `example/sideband_auto_H/{two_RF,three_RF}.json` use
  checkout-relative datasets, contain no H-rate inputs, and fit all three peaks
  over 105–135 ppm at 1.2 GHz. Two RF means 25/100 Hz at one magnetic field;
  three RF means 25/50/100 Hz. Each peak/profile has 147 points (~24.985 Hz).
- **[still applied]** `refit_automatic_proton.py` relocates archived input
  basenames to this checkout, so reproducing the ten-fit study needs no original
  machine paths. It writes a new output directory and records source/input hashes.
- **[still applied]** README, both `SBONEST_MANUAL` files, and `SIDEBAND.md`
  describe the portable automatic workflow. Inherited `MANUAL` files point to
  the SBONEST manuals. The original ONEST handoff is preserved locally as
  `.archive/onest-handoff-before-sbonest-20261002.md`.
- **[still applied]** `.gitignore` excludes environments, backups, raw Bruker
  files, diagnostics, and historical manuscript drafts/renderings. Selected
  manuscript DOCX files, evidence, and QA records are published; consult
  `manuscript/sideband_30ppm/README.md`. Small synthetic fixtures under
  `manuscript/sideband_30ppm/results/1200/` remain for older examples.
  `results/auto_H_refit/fits/executed_source/` preserves scientific provenance.
- **[still applied]** `compare_field_models.py` and
  `results/field_comparison_600_800_20261003_02/` preserve the separate 600/800 MHz
  study, including portable input copies, three-start fits, noiseless controls,
  full covariance, source hashes, and PNG/PDF figures. Root and bilingual manuals
  explain its fixed-H conditions separately from automatic H-rate fitting.
- **[still applied]** GitHub CI installs Python 3.12 and PyPI `optimalcontrol-nmr`,
  checks root Python files, runs inherited and sideband regression checks, the
  synthetic demo, and the portable two-RF CLI fit.

## Evidence and commands

From a new checkout, create a Python 3.12 environment and install
`requirements-sideband.txt`. A sibling OC checkout is unnecessary.

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg
.venv/bin/python run.py example/sideband_auto_H/two_RF.json
.venv/bin/python test_performance.py
.venv/bin/python test_debugging.py
.venv/bin/python verify_3state.py
.venv/bin/python test_sideband.py
.venv/bin/python demo_sideband.py --out session_artifacts/demo_new
```

Choose a new `Project Name`/output directory before rerunning; existing fit
outputs are protected. Dataset paths resolve from the JSON directory; output
prefixes resolve from the working directory.

`results/auto_H_refit/fits/REFIT_REPORT.md` records the completed ten-fit study:

| RF (Hz) | Points | Free parameters | kex (s⁻¹) | RF scale | Reduced χ² |
|---|---:|---:|---:|---:|---:|
| 25/100 | 882 | 24 | 299.68025 | 1.0822687 | 0.93569 |
| 25/50/100 | 1323 | 24 | 299.49797 | 1.0806102 | 0.94625 |

Three random H-rate starts per RF design agree in χ² within 2 × 10⁻⁷.
Two noiseless controls recovered all 24 generating parameters. The report,
`summary.json`, `predictions.csv`, `metadata.json`, and `full_profile_fit.png`
are under `results/auto_H_refit/fits/`. H rates are weakly determined and some
reach bounds; do not present them as measured precise proton relaxation rates.

A clean Python 3.12 environment with the PyPI OC distribution was installed in
`session_artifacts/github_publish_venv`. The index-only checkout is in
`session_artifacts/github_publish_checkout`; its validation logs are in
`session_artifacts/github_publish_checks/`. These local verification artifacts
are intentionally not distributed. The clean-checkout lint/compile checks,
`test_performance.py`, `test_debugging.py`, `verify_3state.py`, `test_sideband.py`,
synthetic demo, and portable two-RF CLI (including PDF output) all passed.
The portable fit recovered kex = 299.6802480565636 s⁻¹ and χ² = 802.8203212916861,
with 882 observations and rank 24/24. README/manual local links were checked.
Those paths describe the initial release's local checks. CI is the portable
verification record. Archived scientific inputs and executed sources retain
their calculation-time bytes and hashes.

## 600/800 MHz comparison

The fields are fitted separately using the same 105–135 ppm observations within
each model comparison: 225 points at 600 MHz and 297 at 800 MHz, with nominal
nitrogen RF = 25/50/100 Hz. R1H = 2 and R2H = 25 s⁻¹ are fixed. All models share
eight free parameters, error weights, bounds, RF-scale handling, and optimizer.
This is distinct from the current 24-parameter automatic-H examples.

```bash
.venv/bin/python compare_field_models.py \
  --source results/field_comparison_600_800_20261003_02/inputs \
  --out results/field_comparison_repeat_01
```

The explicit `--source` uses tracked copies; default historical manuscript
600/800 MHz inputs are local only. Always use a new output directory. The script
generates fit JSON, CSV, figures, predictions, and input copies. `REPORT.txt` and
`verification.json` in the published bundle were written separately after review.

The noisy 600 MHz fits are nearly identical. At 800 MHz, noisy ONEST kex is
301.414 s⁻¹ versus SBONEST 299.850 s⁻¹; the noiseless ONEST displacement is
+0.51% from 300 s⁻¹. SBONEST recovers the noiseless truth at both fields.
ONEST uses its actual zero-clamped Matrix worker; `N_signed_control` isolates
that clamp. Three starts per model/noise/field produce 36 fits. Noisy ONEST starts
reach slightly different minima; the lowest χ² is retained without asserting
global optimality. See `REPORT.txt` and manual section 7.1 for uncertainty and limits.

All 36 fits were independently reproduced from bundled inputs. Input/generator
and J=0/clamp controls, full-rank covariance, saved objectives/predictions, and
source hashes were checked. A tracked-files-only snapshot passed input/control
checks. The four regression scripts and PR/main CI passed. The full field-study
driver is not part of CI. Published CSV line endings were normalized to LF with
cell-value equality verified; the original CRLF file remains in local `.archive/`.

## Pitfalls already identified

- Archived configs/results may contain original-machine absolute paths and
  calculation-time hashes. **[still applied: archives preserved]** Use the
  portable example configs; do not rewrite scientific provenance to current hashes.
- Older `check_two_rf.py` checks historical source hashes and is not the release
  smoke test. Use `test_sideband.py`, the portable CLI, or the new refit driver.
- The first field-comparison attempt stopped because a temporary assertion
  required all ONEST starts to reach identical minima. **[reverted]** That
  assumption was removed; the completed driver records objective spreads and
  `multistart_agreement`. The partial run is preserved in
  `.archive/field_comparison_600_800_initial_partial_20261003/`.
- `git push` as `dleess` was rejected for this repository on 2026-10-03.
  `Gohyang-Matzip` had write access and completed the publication using scoped
  authentication. Check account permissions before retrying; avoid changing
  global credentials for unrelated repositories. No credentials are stored here.
- Local `.venv/bin/python` is Python 3.14 and lacks the `pip` module. The clean
  publication environment was installed with `uv pip install --python ...`.
  **[still applied: original environment preserved]** Do not assume every local
  interpreter has pip; the documented fresh `python3 -m venv` workflow supplies it.
- The inherited web interface, `prepare.py`, and `mcrun.py` are not supported
  Sideband workflows. Profile likelihood/bootstrap are not built-in Sideband CLI commands.
- The ±2400 Hz demo spans ~39.48 ppm at 1.2 GHz; use `sideband_auto_H` for 30 ppm.

## Future maintenance

1. Run `git status --short --branch` and `git remote -v`. Preserve untracked
   work, local manuscripts, and `.archive/`; do not destructively delete artifacts.
2. Inspect `gh pr list --repo Gohyang-Matzip/sbonest --state all` and
   `gh run list --repo Gohyang-Matzip/sbonest`. The release and comparison PRs
   above are complete; act on subsequent PRs only within the current user task.
   Require passing CI before merging a reviewed head.
3. Fetch `origin`, switch to `main`, and fast-forward to `origin/main` after merge.
   Verify the merged SHA's CI, repository visibility, and documentation links.
4. For documentation changes, verify local links, bilingual agreement, numbers
   against result JSON, and archived source hashes. Follow `AGENTS.md`.
5. No further work is authorized by this handoff alone.
   Experimental validation is a future user task, not a completed result.

## Scientific limits

Nitrogen relaxation is shared across files for a peak; there is no independent
per-field nitrogen relaxation. The NH model omits extra proton spins, separate
water/proton exchange, and CSA–DD cross-correlation. Automatic H-rate fitting
cannot repair wrong proton shifts, RF calibration, or pulse timing. Whether
sidebands improve real 1.2 GHz fits remains unverified experimentally. The
local covariance uses absolute intensity errors, with no reduced-χ² rescaling;
boundary errors are not formal confidence intervals.
