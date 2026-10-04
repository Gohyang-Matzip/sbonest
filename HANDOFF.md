# HANDOFF: SBONEST round 2 — PR 6 (interval coverage, joint fields, model selection)

**Written:** 2026-10-04 (Asia/Seoul)
**Repository:** https://github.com/Gohyang-Matzip/sbonest
**Main checkout:** `/Users/donghanlee/work/projects/sbonest`
**Implementation worktree:** `/Users/donghanlee/work/projects/sbonest/session_artifacts/improve7_20261004/tree`
**Branch:** `codex/coverage-fields-models` (stacked on `codex/perf2-diagnostics-design`, PR #13)
**Programme plan:** `docs/superpowers/plans/2026-10-04-full-improvements.md` (Round 2)

## Goal and current status

Round 2 of the approved improvements. PR 5 is PR #13. This handoff covers PR 6:
coverage of profile and bootstrap intervals, the joint 600/800 MHz study and
two-versus-three-state model selection. Implementation and local validation are
complete; after PR #13 merges this branch is rebased, pushed, CI-checked and
merged. Previous handoff: `.archive/HANDOFF.before-round2-pr6-20261004.md`.

## Applied changes (PR 6)

Every item below is **still applied**.

- **[still applied] `sb_bootstrap.sample_row`** factored out; `bootstrap_replicate(analyses=)`
  runs extra analyses on a replicate's data while it is in place.
- **[still applied] `validate_uncertainty.py`:** `--profile-interval NAME...`,
  `--inner-bootstrap B`, `--workers N`; `study_replicate`, `_replicate_analyses`,
  `_interval_coverage`; `coverage.json` schema 2 with an `intervals` block (profile
  and bootstrap coverage, Wilson intervals, mean widths, open/failed counts);
  replicates run in the pool through `sb_parallel.study_task`; parallel and serial
  results identical.
- **[still applied] `compare_joint_fields.py` (new) and
  `results/field_comparison_joint_20261004/`:** joint fits of the archived 600/800
  benchmark with four variants (fixed H; automatic H unbounded; automatic H bounded
  50/500 s⁻¹; plus field-group N). Joint fixed-H fit kex 298.9 ± 0.9 (separate
  297.9 ± 1.6 / 299.9 ± 1.2); automatic H leaves kex/pB unchanged but the proton
  rates are undetermined at these fields (unbounded per-field-N run drifted to R2H
  ≈ 2.5e4 s⁻¹ with biased N rates, kept in
  `session_artifacts/improve7_20261004/joint_try_unbounded`); AICc prefers fixed H
  (Δ +7.8 / +12.5). `.gitignore` excludes checkpoints and input copies of the
  archived run.
- **[still applied] `sb_compare.three_state_config`, `compare_models`, `models_lines`;
  `sb_workflow.py compare --models ... [--h-ppm-c LABEL=ppm]`:** derived three-state
  configurations with five explicit multistart starts, AICc/BIC comparison, no
  F-test (boundary nesting), degeneracy warning when a return rate hits its bound.
  Synthetic three-state data: Δ AICc −26,307 and recovered rates; two-state data:
  degenerate three-state fit flagged, AICc +6.9 for two states.
- **[still applied] Tests:** `test_uncertainty_validation.py::check_profile_and_bootstrap_coverage`,
  `test_sb_compare.py::check_model_comparison`.
- **[still applied] Docs:** manual 7.1 (joint study), 8.2 (coverage extension), 8.7
  (model selection), 12.1 (both languages); README, SIDEBAND.md, AGENTS.md, CHANGELOG.

## Findings

- A first look at the unbounded per-field-N variant suggested nondeterminism; a
  direct check (`session_artifacts/improve7_20261004/determinism.py`) showed serial
  and 8-worker runs identical — the runaway values belonged to a different variant.
- Three-state fits on two-state data drive `kcb` to its lower bound, which makes
  state C absorbing and the stationary populations meaningless; the comparison
  warns about this explicitly.

## Evidence

`session_artifacts/improve7_20261004/verification_01/summary.json`;
`results/field_comparison_joint_20261004/summary.json`;
`session_artifacts/improve7_20261004/compare_test.log`.

## Delivery steps — perform only those still missing

1. After PR #13 merges, rebase onto main, rerun the verification, push, open the
   PR with `session_artifacts/improve7_20261004/pr6_body.md`, require CI, merge,
   verify merge CI, fast-forward main.
2. Continue with PR 7 (`codex/import-qa-web-refactor`).

## Scientific boundaries

Coverage studies, the joint benchmark and model selection are synthetic and
conditional on the supplied sigma and the generating/fitted models. The joint
study is one noise draw per field. Model selection by AICc does not prove a
mechanism.
