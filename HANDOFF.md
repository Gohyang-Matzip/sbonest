# HANDOFF: SBONEST diagnostics, uncertainty and fitting performance

**Written:** 2026-10-04 (Asia/Seoul)
**Repository:** https://github.com/Gohyang-Matzip/sbonest
**Main checkout:** `/Users/donghanlee/work/projects/sbonest`
**Implementation worktree:** `/Users/donghanlee/work/projects/sbonest/session_artifacts/improve_20261004/tree`
**Branch:** `codex/fit-diagnostics-and-uncertainty`
**Base:** `96d6dd68ae3991e840e8f90696edbb6498c5f0dc`

## Goal and current status

The user approved six proposed improvements, followed by debugging, refactoring,
handoff, then commit, push, PR and merge. This authorizes the complete remote
workflow, including passing PR CI, merge, merged-commit CI and local main sync.

Implementation, debugging/refactoring, independent review and local validation
are complete. At the time this handoff was written the changes are uncommitted
in the worktree above; remote delivery is the remaining step. This document is
part of the change to publish, so it is not evidence that delivery has happened.
Inspect live Git/GitHub state before repeating any delivery action.

## Applied changes

Every item below is **still applied**.

- **[still applied] `sbfit.py`, `sb_report.py`:** result schema version 2 keeps
  existing scalar fields and adds ordered full covariance, covariance-aware
  `derived_se.kex`/`pB`, calculation-time source/input/waveform SHA-256 hashes,
  config hash, package/Python/platform/thread/Git metadata and elapsed time.
  Nonfinite covariance entries and unavailable errors are JSON null. Fixed
  parameters remain zero-uncertainty assumptions. CSV preserves all observation,
  prediction, sigma and standardized-residual digits; its residual square sum
  matches chi2. CLI fit PDFs include residual panels and fitted-RF labels.
- **[still applied] `fit.py`, `sbfit.py`:** grouped finite differences support
  ordered free subsets and feasible bound stencils. Natural full-order fits keep
  the prior 2-point scheme; partial/reordered fits keep prior 3-point steps and
  accuracy while grouping independent residue parameters. The measured 23-free
  case uses 21 residual evaluations per Jacobian instead of 47; full 24-parameter
  fits retain 11. No physical forward model or numerical tolerance was loosened.
- **[still applied] `sb_analysis.py`:** optional `init.multistart` accepts explicit
  starts and seeded random starts, preserves all attempts and restores the lowest
  converged solution. `init.profile` scans kex, pB or v1n_scale with exact
  constrained nuisance refits, fixed parameters and user/physical bounds.
  Failed points and raw negative delta-chi2 are retained and warned about. No
  automatic confidence-interval or global-optimum claim is made. All-failed
  multistart runs exit nonzero but preserve attempts/provenance in result JSON.
- **[still applied] `benchmark.py`:** dispatches Sideband configs correctly;
  retains legacy ONEST and `<config> [profile]`; adds `--profile-output PATH`.
  Existing profiles are protected with early checks and exclusive creation.
- **[still applied] `constraints-sideband.txt`, CI:** tested Python 3.12 dependency
  snapshot and four new executable regression scripts. No new runtime library
  family was introduced.
- **[still applied] `test_grouped_jacobian.py`, `test_benchmark.py`,
  `test_sb_analysis.py`, `test_sb_output.py`:** independent derivative/reference
  tests, actual CLI checks, analytical constrained fits, failure/state/seed/bound
  tests, provenance bytes, correlated errors, CSV reconciliation and overwrite
  protection. Real-data partial fitting is compared to the prior SciPy 3-point
  path at rtol=1e-8/atol=1e-8 for parameters and chi2.
- **[still applied] README, both SBONEST manuals, SIDEBAND, AGENTS:** synchronized
  setup/CLI/output/uncertainty documentation. `docs/diagnostics-preview.png` is
  the A1 page from the existing synthetic two-RF example; it is not new
  experimental evidence. All three generated fit-PDF pages were rendered and
  inspected. Local documentation links were checked.
- **[still applied] `HANDOFF.md`:** this briefing. The prior loader handoff is
  archived in the worktree at `.archive/HANDOFF.before-improve-20261004.md`.

No files in historical `results/`, manuscript bundles or raw inputs were changed.
The root `.venv` was preserved. Temporary environments, runs and worktree are
ignored local artifacts and must not be added to Git.

## Debugging and refactoring findings

- An initial grouped-forward derivative for partial fits was **replaced**, not
  shipped: it changed chi2 by 3.06115e-5 and rank in a real-data comparison. Its
  evidence remains in `session_artifacts/improve_20261004/partial_after*` and
  `performance.json`. Do not use those as the final performance/parity result.
  The final grouped 3-point comparison changes chi2 by only 1.53079e-9 and the
  largest parameter by 1.37803e-10; both paths report rank 22 of 23 for that
  particular boundary-start example. Null covariance for this example is an
  identifiability diagnostic, not an error to suppress.
- Feasible profile rates can reconstruct one ULP outside a bound. Only actual
  IEEE-scale rate/interval roundoff is repaired; infeasible small-rate targets
  are rejected. A `max(1,rate)` tolerance floor was found too broad and removed.
- A rejected default start originally blocked valid explicit restarts, and a
  reused model could retain stale fixed values. The model now publishes fresh
  initialization/bounds before the start-feasibility check; stale preparation
  cannot be consumed. New/reused real-model regressions verify this behavior.
- Dangling output symlinks were not caught by `Path.exists()`. They now fail
  before writing; numerical reports and CSV use exclusive creation.
- Report/provenance logic lives in a small separate module; the original NH
  propagation and inherited ONEST reporting remain intact. No bulk formatting of
  executed research sources was performed.

## Verified evidence and commands

Use the existing Python 3.12.14 environment:

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg
PY=/Users/donghanlee/work/projects/sbonest/session_artifacts/repo_review_20261004_01/.venv/bin/python
cd /Users/donghanlee/work/projects/sbonest/session_artifacts/improve_20261004/tree
"$PY" test_performance.py
"$PY" test_debugging.py
"$PY" verify_3state.py
"$PY" test_sideband.py
"$PY" test_grouped_jacobian.py
"$PY" test_benchmark.py
"$PY" test_sb_analysis.py
"$PY" test_sb_output.py
```

The 12 entries in the main checkout's
`session_artifacts/improve_20261004/verification_03/summary.json` all exited zero:
Ruff F, compilation, the eight scripts above, demo, and real automatic-H CLI
including PDF output. Logs and fresh outputs are in that directory. The full
fit has 882 points, 24 parameters, rank 24 and chi2=802.8203212916861, exactly the
prior local baseline. CSV contents equal the preceding inspected PDF run.

`session_artifacts/improve_20261004/partial_final_parity.json` stores the final
partial-fit comparison; `measure_final.py` and `performance_final.json` record
Jacobian timings and evaluation counts. Timing applies to one Jacobian on this
machine, not every full fit. Independent reviews found no unresolved actionable
issue, including a final review of grouped 3-point wiring and small-rate bounds.

A fresh checkout should create a Python 3.12 environment and install:
`python -m pip install -r requirements-sideband.txt -c constraints-sideband.txt ruff`.
The original root `.venv/bin/python` is Python 3.14; do not mistake it for CI.
New output prefixes are mandatory. Dataset/waveform paths resolve from config;
output prefixes resolve from the current working directory.

## Delivery steps — perform only those still missing

1. Check both worktrees and remote state. At handoff, main and origin/main were
   the base SHA above; no user-owned dirty files existed. Do not reset or discard
   subsequently appearing user edits.
2. Stage only implementation, test, constraints, CI, docs/preview and this handoff;
   commit using the repository's conventional prefix, then push the branch.
3. GitHub account `Gohyang-Matzip` has verified ADMIN access. Global gh account
   is `dleess`; preserve it. The local helper
   `/Users/donghanlee/work/projects/sbonest/session_artifacts/improve_20261004/github.py`
   obtains the correct account token into a subprocess environment without
   printing or saving it. Use it to run scoped gh commands and Git with
   `-c credential.helper= -c 'credential.helper=!gh auth git-credential'`.
4. Find an existing PR for `codex/fit-diagnostics-and-uncertainty` before creating
   one. Use a body file, include numerical validation and the committed preview.
   Require successful CI for its exact head SHA, then merge into main.
5. Verify successful push CI for the merge SHA, fetch, and fast-forward the clean
   main checkout to origin/main. Report actual PR URL, merge SHA and CI outcome.
6. Retain or archive the ignored worktree and evidence; do not destructively
   remove artifacts. Stop after delivery; no further research work is authorized.

Useful read-only status commands:

```bash
git status --short --branch
gh pr list --repo Gohyang-Matzip/sbonest --state all \
  --head codex/fit-diagnostics-and-uncertainty \
  --json number,state,url,headRefOid,mergeCommit
gh run list --repo Gohyang-Matzip/sbonest --branch main --limit 5
```

## Scientific boundaries

All new validation uses existing or small independently generated synthetic
inputs. Local errors/profile differences are conditional on the model and fixed
inputs. Weak H-rate identifiability, boundaries and model mismatch still require
scientific interpretation. Bootstrap and automatic profile-based confidence
intervals are not implemented. The historical fixed-H 600/800 MHz benchmark was
not rerun or altered.
