# HANDOFF: SBONEST round 2 — PR 5 (parallel 2, residual diagnostics, optimal design)

**Written:** 2026-10-04 (Asia/Seoul)
**Repository:** https://github.com/Gohyang-Matzip/sbonest
**Main checkout:** `/Users/donghanlee/work/projects/sbonest`
**Implementation worktree:** `/Users/donghanlee/work/projects/sbonest/session_artifacts/improve6_20261004/tree`
**Branch:** `codex/perf2-diagnostics-design` (base: main `4f88171`)
**Programme plan:** `docs/superpowers/plans/2026-10-04-full-improvements.md` (Round 2 section)

## Goal and current status

The user approved the second improvement list ("제안 한대로 하고 debugging,
refactoring 하고 commit push PR merge 해"). Round 2 is delivered as three PRs;
this handoff covers PR 5. Implementation and local validation are complete; the
branch is to be pushed, CI-checked and merged. The previous handoff is archived
at `.archive/HANDOFF.before-round2-20261004.md`.

## Applied changes (PR 5)

Every item below is **still applied**.

- **[still applied] `sb_parallel.py`:** `WorkerPool.__enter__` submits one
  `warm_task` per worker and waits, so all processes exist before the first
  Jacobian (previously the executor spawned them lazily: first batch of ten
  evaluations 0.73 s against 0.18 s in steady state). `predict_task` and
  `predict_blocks` evaluate residue/dataset blocks with balanced point counts.
- **[still applied] `estmodel.py`:** matplotlib is imported inside
  `_plot_residues`; `import sbfit` no longer imports matplotlib (0.45 → 0.37 s,
  and the same saving per worker).
- **[still applied] `SidebandModel.errFunc` override:** with a pool and at least
  two data blocks, block predictions come from the workers and the residual
  arithmetic stays in the main process; the `evaluate_many` closure routes a
  single vector through the memoized/block-parallel residual. The 882-point
  example: 31 s with one worker, 10 s with eight (was 16 s), chi2
  802.8203212916861 and parameters identical to the archived result.
- **[still applied] `sb_diagnostics.py` (new):** `runs_test`,
  `lag1_autocorrelation`, `residual_diagnostics(rows, n_parameters, covariance)`,
  `diagnostics_lines`. `run_config` stores `residual_diagnostics` and per-parameter
  `stderr_rescaled` and merges warnings; the text report appends the diagnostics
  (`getLogBuffer` keeps the covariance for rescaled errors); `sb_report`
  recomputes them from the CSV for regenerated reports.
- **[still applied] `sb_design.py`:** `optimize` section (scenario, budget,
  criterion kex/pB/parameter/D, min_per_dataset), `_criterion_gradient`,
  `_measurement_groups`, `select_measurements` (backward elimination with exact
  Woodbury downdates over whole spectrum rows), `optimize_design` producing
  `<name>_optimized_<budget>` and `<name>_uniform_<budget>` scenarios, text and
  plot. Bundled example at sigma 0.01: 60 of 294 rows → kex SE 7.6 s⁻¹ versus
  17.5 s⁻¹ uniform and 6.7 s⁻¹ for all rows; the incremental criterion equals the
  re-evaluated variance.
- **[still applied] CI/metadata:** three check shards (core, analysis, tools) on
  Python 3.12, 3.13 and 3.14 plus the 3.12 workflow job; `CHANGELOG.md`,
  `CITATION.cff`, `pyproject.toml` version 1.2.0. The pinned constraints install
  on Python 3.14.7 locally and the pool/server/model/checkpoint/MC scripts pass
  there.
- **[still applied] Docs:** manual 8.4 (numbers and mechanism), 8.6 (optimize),
  new 8.8 (diagnostics) in both languages; README, SIDEBAND.md, AGENTS.md.
- **[still applied] Tests:** `test_sb_diagnostics.py` (new),
  `test_sb_design.py::check_optimization`.

## Evidence

`session_artifacts/improve6_20261004/verification_01/summary.json` (ruff,
compile, 23 scripts, demo, `git diff --check`); `tree/session_artifacts/perf2`
(1 and 8 worker fits, identical to the archive); `tree/session_artifacts/design_opt`
(optimized design example); `py314_*.log`.

## Delivery steps — perform only those still missing

1. Push, open the PR with `session_artifacts/improve6_20261004/pr5_body.md`,
   require CI for the exact head, merge, verify merge CI, fast-forward main.
2. Continue with PR 6 (`codex/coverage-fields-models`).

## Scientific boundaries

Diagnostics are indicators: a reduced chi-square away from 1 does not identify
whether sigma or the model is wrong, and rescaled errors assume a correct model.
Optimized designs are local to the assumed truth and criterion. All evidence is
synthetic.
