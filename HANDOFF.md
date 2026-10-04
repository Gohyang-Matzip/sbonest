# HANDOFF: SBONEST improvement programme — PR 1 (parallel workers and diagnostics)

**Written:** 2026-10-04 (Asia/Seoul)
**Repository:** https://github.com/Gohyang-Matzip/sbonest
**Main checkout:** `/Users/donghanlee/work/projects/sbonest`
**Implementation worktree:** `/Users/donghanlee/work/projects/sbonest/session_artifacts/improve2_20261004/tree`
**Branch:** `codex/parallel-and-diagnostics`
**Base:** `791fb1dee0f1d8b2e1a0d34fd43b264044c18a73` (main after PR #8)

## Goal and current status

The user approved the complete 2026-10-04 improvement list ("제안대로 전부 해").
The programme is delivered as four themed PRs described in
`docs/superpowers/plans/2026-10-04-full-improvements.md`; this handoff covers
PR 1. Each PR follows the same delivery rule: executable regressions, bilingual
documentation, local verification, PR CI, merge, merged-commit CI, main sync.

PR 1 implementation and local validation are complete; at the time of writing
the changes are uncommitted in the worktree above. This file is part of the
change, so it is not evidence that delivery happened. Inspect live Git/GitHub
state before repeating any delivery action. The previous handoff (PR #8,
reliable analysis workflow) is archived locally at
`.archive/HANDOFF.before-parallel-20261004.md`.

## Applied changes (PR 1)

Every item below is **still applied**.

- **[still applied] `sb_parallel.py` (new):** spawn-based `WorkerPool` holding one
  `SidebandModel` per worker with cached prepared data; module-level tasks
  `residual_task`, `multistart_task`, `profile_task`, `bootstrap_task` reuse the
  serial code paths; `map` yields in submission order; `Progress` prints
  elapsed time and an ETA. `validate_workers` rejects non-positive/bool/float.
- **[still applied] `fit._block_jacobian(evaluate_many=...)`:** the base vector is
  evaluated first (its dtype selects the 3-point epsilon exactly as before), then
  all perturbed vectors are evaluated in one batch. Stencils, steps and arithmetic
  are unchanged; `test_grouped_jacobian.py` and `test_profile_jacobian.py` pass.
- **[still applied] `SidebandModel.fit`:** when `model.pool` is set, Jacobian
  columns are evaluated in workers. The residual closure memoizes the last
  evaluated vector: SciPy evaluates the residual at an accepted step and then
  requests the Jacobian at the same point, so the base evaluation is reused.
  Results are identical (882-point example: chi2 802.8203212916861, same
  parameters and covariance); serial time 36 s → 33 s, 8 workers 16 s.
- **[still applied] `sb_analysis.py`:** `fit_attempt`, `profile_point`,
  `base_fit_config`, `local_jacobian`, `identifiability`; `fit_multistart` and
  `profile_likelihood` accept `pool=` and run pending items in index order
  (the configured initial fit always runs in the main process). The profile
  residual closure also memoizes its last vector. Completed-prefix and
  `on_complete` semantics are unchanged.
- **[still applied] `sb_bootstrap.py`:** `bootstrap_draws` (SeedSequence([seed,
  index]) draws in the main process) and `bootstrap_replicate` (fit one replicate,
  restore data); `bootstrap_fit(pool=)` sends drawn observations to workers.
- **[still applied] `sbfit.run_config(workers=)`, `check_config(identifiability=,
  workers=)`, `run.py`/`sbfit.py` flags `--workers N` and `--identifiability`
  (requires `--check`).** Provenance records `workers`; the checkpoint identity
  excludes it, so a parallel run resumes serially and vice versa.
- **[still applied] Identifiability report:** column-scaled singular values, rank,
  condition, expected and relative SE from `(J^T J)^-1`, zero-sensitivity and weak
  parameters (relative SE > 100% or no finite SE), correlations ≥ 0.95 and derived
  kex/pB SE at the initial point; `--check` adds warnings for rank deficiency,
  weak parameters and strong correlations. On the bundled example it names
  `A1.R1H, G2.R1H, S3.R1H` and the three R1H/R2H pairs before any fit.
- **[still applied] `sb_report._paginate`:** section-aware, balanced pagination of
  the report summary text (no near-empty orphan page).
- **[still applied] CI:** `checks` job (ruff, compile, 15 fast scripts) on Python
  3.12 and 3.13; `workflow` job (test_sb_output, demo, dummy-guide commands with
  `--check --identifiability`, `--workers 2` fit, serial `--resume`, report) on 3.12.
  The pinned constraints install and the fast scripts pass locally on 3.13.15.
- **[still applied] `test_sb_parallel.py` (new):** worker validation, Progress
  output, exact Jacobian/identifiability parity through a 2-worker pool, task
  error propagation and pool reuse, serial-vs-parallel `run_config` with
  multistart/profile/bootstrap (result JSON and every checkpoint record identical),
  resume across worker counts without refitting, CLI flag errors and
  `--check --identifiability --workers 2` output.
- **[still applied] Docs:** README, both SBONEST manuals (section 8.4, --check
  identifiability, troubleshooting), both dummy guides, SIDEBAND.md, AGENTS.md.

## Findings

- The first full-fit timing of the serial path showed 394 s wall time; the
  result's own `elapsed_s` was 38 s. A concurrent Python 3.13 environment build
  was saturating the machine. A clean rerun took 36 s (33 s with memoization).
- Spawned workers re-import the caller's `__main__`; a script without an
  `if __name__ == "__main__":` guard recursively starts pools. All repository
  entry points have guards; the manual documents the requirement.
- Single-fit speedup is bounded by the optimizer's serial residual evaluations and
  the base Jacobian evaluation; bootstrap/profile/multistart items scale almost
  linearly up to the item count.

## Verified evidence and commands

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg
cd /Users/donghanlee/work/projects/sbonest/session_artifacts/improve2_20261004/tree
ruff check *.py --select F
.venv/bin/python -m py_compile *.py
for t in test_performance test_debugging verify_3state test_sideband \
         test_grouped_jacobian test_benchmark test_sb_analysis test_sb_output \
         test_sb_check test_sb_checkpoint test_profile_jacobian test_sb_bootstrap \
         test_sb_workflow test_uncertainty_validation test_sb_parallel; do .venv/bin/python $t.py; done
.venv/bin/python demo_sideband.py --out session_artifacts/final_demo_NEW
```

Evidence under `session_artifacts/improve2_20261004/` (main checkout, ignored):
`verification_01/summary.json` (all commands exit 0), `parity_small.py`
(serial vs 3 workers identical result JSON), `tree/session_artifacts/par_full/`
(882-point fits with 1/4/8 workers, all chi2 802.8203212916861 and identical
parameters), `py313_*.log` (fast scripts on Python 3.13.15), `guide_workers/`
(dummy-guide workflow with `--workers`, rendered report pages).

## Delivery steps — perform only those still missing

1. Check worktree and remote state; do not discard user edits in the main checkout.
2. Stage implementation, tests, CI, docs, plan and this handoff; commit with
   `perf:`/`feat:` prefix; push `codex/parallel-and-diagnostics`.
3. Use `session_artifacts/improve_20261004/github.py` for `Gohyang-Matzip`-scoped
   gh/Git commands; preserve the global `dleess` account.
4. Find or create the PR; require CI success for the exact head SHA; merge.
5. Verify merge-commit CI, fast-forward the clean main checkout, report.
6. Continue with PR 2 (`codex/design-intervals-comparison`) from the merged main.

## Scientific boundaries

Parallel execution changes only where residuals are evaluated. Identifiability
diagnostics are local linear properties of the initial point and the design with
the supplied absolute sigma; they do not guarantee convergence or replace the
fitted covariance. All evidence is synthetic.
