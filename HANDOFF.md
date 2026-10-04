# HANDOFF: SBONEST reliable analysis workflow

**Written:** 2026-10-04 (Asia/Seoul)
**Repository:** https://github.com/Gohyang-Matzip/sbonest
**Main checkout:** `/Users/donghanlee/work/projects/sbonest`
**Implementation worktree:** `/Users/donghanlee/work/projects/sbonest/session_artifacts/workflow_20261004/tree`
**Branch:** `codex/reliable-analysis-workflow`
**Base:** `e19fbb14b0ded0a988d9191db597aef0051307a2` (main and origin/main at handoff)

## Goal and current status

The user approved six proposed improvements (preflight checks, durable and
resumable analysis, one-command dummy setup, grouped profile derivatives,
saved-result reporting, and a reproducible parametric bootstrap with synthetic
uncertainty checks), followed by debugging, refactoring, this handoff, then
commit, push, PR and merge. This authorizes the complete remote workflow,
including passing PR CI, merge, merged-commit CI and local main sync.

Implementation, debugging/refactoring, an independent branch review and local
validation are complete. At the time this handoff was written the changes are
uncommitted in the worktree above; remote delivery is the remaining step. This
document is part of the change to publish, so it is not evidence that delivery
has happened. Inspect live Git/GitHub state before repeating any delivery action.

## Applied changes

Every item below is **still applied**.

- **[still applied] `sbfit.py`, `run.py`, `test_sb_check.py`:**
  `SidebandModel.prepare_fit` initializes/validates vectors, bounds and free
  indices without optimizing; `check_config` and `run.py CONFIG --check` print a
  JSON summary (validity, errors, warnings, resolved datasets/waveforms, planned
  outputs and conflicts, point/residue/parameter counts, free/fixed names and
  bounds) and write nothing. Invalid optional `init.multistart`, `init.profile`
  and `init.bootstrap` settings fail before any optimizer runs. `--check` and
  `--resume` exit with an error for non-Sideband methods and cannot be combined.
- **[still applied] `sb_checkpoint.py`, `sbfit.py`, `test_sb_checkpoint.py`:**
  every Sideband run owns `<Project Name>_checkpoint/` with an flock'd `.lock`,
  a manifest holding the execution identity (config/input/waveform/source
  SHA-256, packages, Python, platform, threads, PDF mode), provenance, baseline
  snapshot, per-attempt/per-point/per-replicate records, staged exports and a
  completion marker. Records are immutable JSON with checksums; staged outputs
  are published exclusively, and resume only accepts identical already-published
  files. `--resume` requires the same identity, skips the completed baseline,
  restarts, profile points and replicates, and a run that already finished with
  failure is restored and exits nonzero without refitting. Failed fits retain
  optimizer diagnostics and provenance in `<prefix>_result.json`.
- **[still applied] `sb_analysis.py`, `test_sb_analysis.py`,
  `test_profile_jacobian.py`:** multistart and profile scans accept `completed`
  records and `on_complete` callbacks (callback exceptions and interruptions
  propagate); `snapshot_fit`/`restore_fit` restore a selected fit exactly
  (covariance, rank, diagnostics, predictions). Profile refits group verified
  residue-independent coordinates and keep the dense reference for unknown
  layouts; constrained exchange coordinates, reordered/subset vary, fixed
  parameters and inward bound stencils are covered against an independent dense
  reference.
- **[still applied] `sb_bootstrap.py`, `test_sb_bootstrap.py`:**
  `init.bootstrap = {replicates, seed, confidence}` draws seeded Gaussian
  replicates around the selected predictions with the supplied absolute sigma,
  refits each from the selected baseline, retains every replicate and status,
  summarizes percentiles only across successful varied estimates, labels fixed
  inputs as assumptions, and warns on failures, boundary hits and fewer than 100
  successful replicates. Model/data state is restored even after interruption.
- **[still applied] `validate_uncertainty.py`, `test_uncertainty_validation.py`:**
  optional seeded synthetic coverage study of local standard-error intervals
  with explicit truth, independent noise draws, retained failures and Wilson
  sampling intervals. Separate from CI's small executable checks.
- **[still applied] `sb_workflow.py`, `sb_report.py`, `test_sb_workflow.py`,
  `test_sb_output.py`:** `python sb_workflow.py init-demo --out NEW_DIRECTORY`
  copies the bundled synthetic two-RF inputs into a new folder and writes
  `fit.json` with an absolute `Project Name`; existing targets are rejected.
  `python sb_workflow.py report RESULT_JSON --out NEW_PREFIX` reads the result
  and its predictions CSV (hashes both, validates consistency), and writes
  `<prefix>_summary.json`, `<prefix>_summary.txt` and `<prefix>.pdf` with
  profile/residual panels, restart/profile/replicate summaries and residue/dataset
  residual statistics. No optimization and no new dependency.
- **[still applied] `.github/workflows/ci.yml`, `AGENTS.md`:** CI runs the six
  new executable regressions and the actual dummy-guide commands (init-demo,
  `--check`, fit, `--resume`, report) once with `--no-pdf`. Contributor guide
  lists the new tests and the checkpoint/resume rules.
- **[still applied] `README.md`, `SBONEST_MANUAL.md`, `SBONEST_MANUAL.ko.md`,
  `DUMMY_GUIDE.md`, `DUMMY_GUIDE.ko.md`, `SIDEBAND.md`:** synchronized
  documentation of the helper CLI, preflight, checkpoints/resume, reports and
  bootstrap. Both guides contain the same 12 bash blocks; blocks 3–12 were
  executed verbatim (see evidence).
- **[still applied] `docs/superpowers/plans/2026-10-04-reliable-analysis.md`:**
  the implementation plan with its execution record; a previous plan is already
  tracked under the same directory.
- **[still applied] `HANDOFF.md`:** this briefing. The previous handoff
  (PR #6 diagnostics work) is archived locally, outside Git, at
  `.archive/HANDOFF.before-workflow-20261004.md` in the worktree.

No files in historical `results/`, manuscript bundles or raw inputs were changed.
The root `.venv` of the main checkout was preserved. Worktree environments, run
folders and evidence under `session_artifacts/` are ignored and must not be added.

## Debugging and refactoring findings

- **Relative output prefix regression (fixed, still applied).** The final
  verification run of `demo_sideband.py` failed with
  `'/.../noisefree_scale_checkpoint/export-5ixb_21o/noisefree_scale_result.txt' is not in the subpath of 'noisefree_scale_checkpoint'`.
  Cause: Python 3.12 `tempfile.mkdtemp` returns absolute paths while a
  cwd-relative `Project Name` left `Checkpoint.path` relative, so
  `output_plan`'s `relative_to` failed. The regression came in with the
  checkpoint-integrity hardening after the earlier verification pass.
  `Checkpoint.__init__` now stores `Path(path).absolute()`;
  `check_relative_prefix` in `test_sb_checkpoint.py` reproduces the failure on
  the old behavior (verified red) and passes now. The test compares resolved
  paths because macOS `/var` is a symlink to `/private/var`.
- Independent branch review (worker `workflow_review`) found checkpoint
  integrity, failed-run recovery and fixed-population labeling issues; all were
  fixed before handoff and the final verdict was approve. Logs:
  `session_artifacts/workflow_20261004/checkpoint_red.log`,
  `checkpoint_green.log`, `checkpoint_fix.log`, `checkpoint_integrity.log`.
- Grouped profile derivatives were investigated against the real 882-point
  layout (`investigate_profile_layout.py`, `verify_full_profiles.py` in
  `session_artifacts/workflow_20261004/`); grouping is only applied to verified
  residue-independent coordinates, and the dense path remains the reference.

## Verified evidence and commands

Use the worktree's Python 3.12.14 environment (NumPy 2.5.3, SciPy 1.18.1):

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg
cd /Users/donghanlee/work/projects/sbonest/session_artifacts/workflow_20261004/tree
ruff check *.py --select F
.venv/bin/python -m py_compile *.py
for t in test_performance test_debugging verify_3state test_sideband test_grouped_jacobian \
         test_benchmark test_sb_analysis test_sb_output test_sb_check test_sb_checkpoint \
         test_profile_jacobian test_sb_bootstrap test_sb_workflow test_uncertainty_validation; do
  .venv/bin/python $t.py || break
done
.venv/bin/python demo_sideband.py --out session_artifacts/final_demo_03
```

Evidence directory: `/Users/donghanlee/work/projects/sbonest/session_artifacts/workflow_20261004/`.

- `verification_final/summary.json` + `NN.log`: Ruff, compile, 14 regression
  scripts and `git diff --check` exited 0. Entry 16 (`demo_sideband.py`) recorded
  the relative-prefix failure above; `16_rerun.log` (fresh
  `session_artifacts/final_demo_02`) and `11_rerun.log` (`test_sb_checkpoint.py`)
  passed after the fix. `post_fix_summary2.txt` lists the reruns of every script
  that exercises `run_config`/`Checkpoint` after the fix (all exit 0;
  `post_fix_summary.txt` is an aborted earlier attempt).
- `guide_final.log`: DUMMY_GUIDE blocks 3–12 executed verbatim with `set -eu`
  in a fresh `tree/session_artifacts/dummy_01` (check → fit → resume → result
  check → report → multistart/profile/bootstrap analysis → analysis report →
  profiling). The fit gives chi2/dof 0.93569, kex 299.68 ± 0.67 s⁻¹, rank 24 of
  24, `S3.R1H` at its bound with the documented proton-rate warnings.
  Rendered pages of `analysis_report_01.pdf` and `fit.pdf` are in
  `verification_final/*.png` and were inspected.
- `tree/session_artifacts/full_workflow_01/`: the default 882-point/24-parameter
  fit reproduces the archived baseline exactly (chi2 = 802.8203212916861, rank 24;
  `results/auto_H_refit/fits/two_RF_result.json`, SHA-256
  `c60a488b…b060b`), plus report outputs and `full_report.log`.
- `full_profile_comparison_fixed.json`: dense vs grouped kex profiles at
  295/300/305 s⁻¹ from the archived baseline give identical parameters and chi2
  (differences 0.0, constraint error 0.0) with 2.27× lower wall time and 2.26×
  fewer residual evaluations for the grouped path.
  `full_profile_numerics.json` confirms identical refit rows and nfev
  (24/16/16) and gradient agreement ≤ 3.7e-12 across memory layouts.
- Earlier passes: `verification_01/summary.json` (all scripts before the review
  fixes), `full_check.json` (`--check` output for the full example), `full_fit.log`.

Fresh checkouts: create a Python 3.12 environment and install
`python -m pip install -r requirements-sideband.txt -c constraints-sideband.txt ruff`.
The main checkout's root `.venv/bin/python` is Python 3.14; do not mistake it
for CI. New output prefixes are mandatory; dataset/waveform paths resolve from
the config directory and output prefixes from the working directory.

## Delivery steps — perform only those still missing

1. Check the worktree and remote state (`git status --short --branch`,
   `git fetch origin`). Main and origin/main were the base SHA above at handoff;
   the main checkout was clean. Do not reset or discard user edits that appear.
2. Stage only: `.github/workflows/ci.yml`, `AGENTS.md`, `DUMMY_GUIDE.ko.md`,
   `DUMMY_GUIDE.md`, `HANDOFF.md`, `README.md`, `SBONEST_MANUAL.ko.md`,
   `SBONEST_MANUAL.md`, `SIDEBAND.md`, `run.py`, `sb_analysis.py`,
   `sb_bootstrap.py`, `sb_checkpoint.py`, `sb_report.py`, `sb_workflow.py`,
   `sbfit.py`, `validate_uncertainty.py`, `docs/superpowers/plans/2026-10-04-reliable-analysis.md`
   and the seven `test_*.py` files (`test_profile_jacobian`, `test_sb_analysis`,
   `test_sb_bootstrap`, `test_sb_check`, `test_sb_checkpoint`, `test_sb_output`,
   `test_sb_workflow`, `test_uncertainty_validation`). Never add `.venv` or
   `session_artifacts`. Commit with a `feat:` prefix and push the branch.
3. GitHub account `Gohyang-Matzip` has ADMIN access to the repository; the global
   gh account is `dleess` and must stay active. The helper
   `/Users/donghanlee/work/projects/sbonest/session_artifacts/improve_20261004/github.py`
   obtains that account's token into a subprocess environment without printing
   it: `python3 github.py gh pr ...` and
   `python3 github.py git -c credential.helper= -c 'credential.helper=!gh auth git-credential' push ...`.
4. Look for an existing PR for `codex/reliable-analysis-workflow` before creating
   one. Use a body file; include the validation above. Require successful CI for
   the exact head SHA, then merge into main.
5. Verify successful push CI for the merge SHA, fetch, and fast-forward the clean
   main checkout to origin/main. Report the PR URL, merge SHA and CI outcome.
6. Keep the worktree and evidence; archive instead of deleting. Stop after
   delivery; no further research work is authorized.

Useful read-only status commands:

```bash
git status --short --branch
gh pr list --repo Gohyang-Matzip/sbonest --state all \
  --head codex/reliable-analysis-workflow \
  --json number,state,url,headRefOid,mergeCommit
gh run list --repo Gohyang-Matzip/sbonest --branch main --limit 5
```

## Scientific boundaries

All validation uses the bundled synthetic example or small independently
generated synthetic inputs; it checks the software workflow, not an experimental
sample. Local errors, profile differences and bootstrap percentile intervals are
conditional on the selected model, fixed inputs and supplied absolute sigma; few
replicates are demonstrations, not coverage evidence. Weak H-rate
identifiability, boundary parameters and model mismatch still require
scientific interpretation. The historical fixed-H 600/800 MHz benchmark was not
rerun or altered.
