# HANDOFF: SBONEST public release and fitting workflow

**Written:** 2026-10-02 · **Working directory:** `/Users/donghanlee/work/projects/sbonest`

**Repository:** https://github.com/Gohyang-Matzip/sbonest · **Release branch:** `sbonest` · **Target:** `main`

## Goal

Publish the SBONEST program, synthetic examples, numerical evidence, English/Korean
Markdown manuals, and this handoff in a public repository. Complete commit → push
→ PR → merge, verify CI on the merged commit, and synchronize local `main`.
The user explicitly authorized public visibility and the complete merge workflow.

## Implementation status

The implementation and documentation are applied. The public repository exists;
the ONEST baseline `6d178f3e5d6dc82b4cfc7a546d9004c17ff492d9` was pushed to `main`.
This handoff accompanies the SBONEST release PR. Its merge/CI state is external:
read GitHub before deciding whether publication work remains. Once the PR is
merged, successful CI and matching local/remote `main` complete this task.
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
  files, process diagnostics, and manuscript drafts/renderings. Small synthetic
  fixtures under `manuscript/sideband_30ppm/results/1200/` remain for older examples.
  `results/auto_H_refit/fits/executed_source/` preserves scientific provenance.
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
CI is the portable verification record. Source/docs pass `git diff --check`;
archived CSV CRLF endings and solver-log trailing spaces remain untouched.

## Pitfalls already identified

- Archived configs/results may contain original-machine absolute paths and
  calculation-time hashes. **[still applied: archives preserved]** Use the
  portable example configs; do not rewrite scientific provenance to current hashes.
- Older `check_two_rf.py` checks historical source hashes and is not the release
  smoke test. Use `test_sideband.py`, the portable CLI, or the new refit driver.
- Local `.venv/bin/python` is Python 3.14 and lacks the `pip` module. The clean
  publication environment was installed with `uv pip install --python ...`.
  **[still applied: original environment preserved]** Do not assume every local
  interpreter has pip; the documented fresh `python3 -m venv` workflow supplies it.
- The inherited web interface, `prepare.py`, and `mcrun.py` are not supported
  Sideband workflows. Profile likelihood/bootstrap are not built-in Sideband CLI commands.
- The ±2400 Hz demo spans ~39.48 ppm at 1.2 GHz; use `sideband_auto_H` for 30 ppm.

## Resume or confirm publication

1. Run `git status --short --branch` and `git remote -v`. Preserve untracked
   work, local manuscripts, and `.archive/`; do not destructively delete artifacts.
2. Inspect `gh pr list --repo Gohyang-Matzip/sbonest --state all` and
   `gh run list --repo Gohyang-Matzip/sbonest`. If the release PR is still open,
   inspect its diff, require passing CI, and merge the reviewed head.
3. Fetch `origin`, switch to `main`, and fast-forward to `origin/main` after merge.
   Verify the merged SHA's CI, repository visibility, and documentation links.
4. If publication is complete, no further work is authorized by this handoff.
   Experimental validation is a future user task, not a completed result.

## Scientific limits

Nitrogen relaxation is shared across files for a peak; there is no independent
per-field nitrogen relaxation. The NH model omits extra proton spins, separate
water/proton exchange, and CSA–DD cross-correlation. Automatic H-rate fitting
cannot repair wrong proton shifts, RF calibration, or pulse timing. Whether
sidebands improve real 1.2 GHz fits remains unverified experimentally. The
local covariance uses absolute intensity errors, with no reduced-χ² rescaling;
boundary errors are not formal confidence intervals.
