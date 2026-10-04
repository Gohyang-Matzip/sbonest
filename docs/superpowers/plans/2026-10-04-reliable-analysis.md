# Reliable Sideband analysis implementation plan

> **For agentic workers:** Use focused workers under dispatching-parallel-agents for independent file ownership. The coordinator integrates dependent changes and runs whole-branch review. Do not commit until debugging, refactoring and HANDOFF are complete, as requested by the user.

**Goal:** Implement the six approved improvements: preflight checks, durable/resumable analysis, one-command dummy setup, grouped profile derivatives, saved-result reporting, and reproducible parametric bootstrap with synthetic uncertainty checks.

**Architecture:** Keep the physical propagator and default optimizer unchanged. Extract reusable preparation from fitting, persist completed work in a run-owned checkpoint directory, and extend existing analysis functions through optional callbacks/completed records. Small CLI helpers consume existing JSON/CSV; no new runtime dependency.

**Tech stack:** Python 3.12, NumPy/SciPy/Matplotlib, standard library, executable regression scripts.

**Spec:** User-approved proposal in this conversation (2026-10-04): order 1 → 2 → 4 → 3 → 5 → 6, then debugging/refactoring, handoff, commit/push/PR/merge. Independent implementation may overlap; integration and acceptance follow dependencies.

## Constraints and review focus

- Preserve public parameters, units, ordering, physical model and numerical tolerances. Do not alter historical research inputs/results.
- Keep Python 3.12 constraints and both manual/guide languages synchronized.
- Preserve existing output protection. Resume is explicit and requires matching config/input/waveform/source/version identity; reject stale or malformed checkpoints.
- Invalid optional settings must fail before an optimizer runs; `--check` creates no results.
- Interruptions retain each completed attempt/point; a restored fit must reproduce covariance, rank, diagnostics and predictions without re-optimizing completed work.
- Profile grouping must handle constrained exchange coordinates, reordered/subset vary, fixed parameters, and inward bound stencils, with an independent dense reference.
- Bootstrap uses supplied absolute sigma, an explicit seed, fixed data design, fresh synthetic observations per replicate, retained failures, and restores original model/data state even after interruption.
- Bootstrap percentile intervals are conditional on the selected model/fixed inputs; few replicates are demos, not coverage evidence. A separate seeded synthetic coverage experiment reports failures and binomial sampling error.
- Reports validate CSV/JSON consistency and preserve failed profile points/negative deltas. Regeneration must not refit or overwrite prior artifacts.
- Numerical-library threads stay at one. No destructive artifact cleanup. Main workspace and original virtualenv remain untouched.

## Tasks and ownership

- [x] **1. Preflight** — `sbfit.py`, `run.py`, `test_sb_check.py`: extract `SidebandModel.prepare_fit(p0=None, fitting_config=None)` to initialize/validate vectors, bounds and free indices without optimizing; add `check_config` and `--check`. Validate all optional settings early. Report data counts, fields/RF, free/fixed names, resolved paths and output conflicts. Tests fail on current late validation and check side effects/zero fitting.
- [x] **2. Persistence** — coordinator owns `sb_checkpoint.py`, `test_sb_checkpoint.py`, later `run_config` integration. A checkpoint directory owns its manifest and completed stage records; atomically persist records, validate provenance identity, preserve a baseline snapshot before scans. `--resume` skips finished work; output export is recoverable and existing unrelated files remain protected. General fitting failures also retain diagnostics.
- [x] **3. Dummy setup and saved-result report** — `sb_workflow.py`, `sb_report.py`, `test_sb_workflow.py`: `python sb_workflow.py init-demo --out NEW_DIRECTORY` copies bundled inputs and writes fit.json, rejecting existing targets. `python sb_workflow.py report RESULT_JSON --out NEW_PREFIX` reads its predictions CSV, writes summary JSON/TXT and PDF including profile/restart/replicate summaries and residue/dataset residual statistics. No optimization or new dependency.
- [x] **4. Analysis hooks and grouped profile Jacobian** — `sb_analysis.py`, `test_sb_analysis.py`, `test_profile_jacobian.py`: optional `completed` records and `on_complete` callbacks for multistart and profiles; exported `snapshot_fit`/`restore_fit` helpers serialize enough state to restore a selected fit exactly. Callback exceptions/interruptions must propagate. Group only verified residue-independent coordinates; retain dense reference for unknown model layouts. Numerical parity and evaluation-count evidence required.
- [x] **5. Bootstrap** — coordinator owns `sb_bootstrap.py`, `test_sb_bootstrap.py`, `validate_uncertainty.py`: `init.bootstrap={replicates,seed,confidence}`; all settings validated before fitting. Seeded parametric Gaussian draws around selected predictions with supplied sigma; fit from selected baseline, retain every replicate/status, percentile summaries only across successful varied estimates, fixed assumptions identified, failure/low-count warnings. Reuse completion hooks/checkpoints. Synthetic coverage runner uses independent noise draws, explicit truth, and reports empirical local-SE coverage with sampling uncertainty.
- [x] **6. Integration/docs/CI** — update README, both manuals and dummy guides to use helper CLI and explain preflight, checkpoints/resume, reports/bootstrap. Add targeted executable tests and actual guide commands to CI; avoid duplicate expensive full fits.
- [x] **7. Debug/refactor/review** — run full existing/new regressions, full 882-point default fit parity, dense/grouped profile comparisons/timings, interruption/resume byte/numerical checks, seeded bootstrap, rendered report QA; simplify changed code; independent branch review, fix verified findings.
- [x] **8. Handoff/delivery** — archive previous HANDOFF locally, write verified self-contained HANDOFF before committing. Commit scoped changes, push, PR; require exact-head CI before merge, merged-SHA CI afterward, synchronize clean main.

## Execution record

- 2026-10-04: isolated worktree created at `session_artifacts/workflow_20261004/tree`, branch `codex/reliable-analysis-workflow`, base `e19fbb14b0ded0a988d9191db597aef0051307a2`.
- Existing validated Python 3.12 environment: `/Users/donghanlee/work/projects/sbonest/session_artifacts/manual_20261004/tree/.venv/bin/python`; reference full-fit/profiling evidence remains in that worktree's `session_artifacts/dummy_k1fnomib`.
- Authorization: all six improvements and remote merge explicitly approved. Routine design and reversible implementation decisions proceed without another approval gate.
- 2026-10-04: all eight tasks completed in the worktree. Final verification found
  and fixed a relative output-prefix regression in `sb_checkpoint.py` (see
  `HANDOFF.md`); independent branch review approved. Delivery (commit, push, PR,
  CI, merge, main sync) follows the handoff.
