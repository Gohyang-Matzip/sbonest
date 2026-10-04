# HANDOFF: SBONEST debugging and refactoring round (after the folder reorganization)

**Written:** 2026-10-05 (Asia/Seoul) · **Working dir:** `/Users/donghanlee/work/projects/sbonest` · **Branch:** `codex/debug-refactor` (from `main` at `59ab89a`, the merge of PR #17)
**Repository:** https://github.com/Gohyang-Matzip/sbonest · **Last release:** v1.2.0 · **This round's PR:** #18 (https://github.com/Gohyang-Matzip/sbonest/pull/18)
**Previous handoffs:** `.archive/HANDOFF.before-debug-refactor-20261005.md` (folder reorganization, PR #17), `.archive/HANDOFF.before-restructure-20261005.md` (round 3 proposals). `.archive/` is gitignored.

## Goal

The user asked for "debugging, refactoring, handoff" after the folder
reorganization. Two audit agents (bug hunt with reproductions; behaviour-preserving
refactoring candidates) reported; every confirmed bug was fixed with a regression
check, the safe refactors were applied, and the round is delivered as one PR
following the established pattern (CI green, merge, fast-forward `main`).

## Bugs fixed (all still applied, each with a regression check)

| # | Defect | Fix | Check |
|---|---|---|---|
| 1 | `sb_compare.run_comparison`/`compare_models`/`three_state_config` took residue labels from `config.residues`; residues present in the data but unlisted are active (manual 8.7), so "individual" fits silently kept them (wrong AICc/BIC/F) and the three-state derivation failed on a missing `h_ppm_c` | `_dataset_labels(config, config_dir)` reads the data with `EstDataSet` + `run.set_residue_flags`; `three_state_config(..., labels=None)` defaults to the shift table minus explicit `off` | `tests/test_sb_compare.py::check_unlisted_residues` (patches `_fit`, no fitting) |
| 2 | `sb_design` with `nitrogen_relaxation.mode = "per_field"` and 2+ field groups always raised "truth must include A1.R2a" (placeholder/dataset writers demanded ungrouped names while the model uses `A1.R2a[g]`) | `_truth_value(mapping, label, key, group)` accepts `label.key` or `label.key[group]`; `_write_dataset` uses `model.dataset_group[index]` | `tests/test_sb_design.py::check_per_field_design` |
| 3 | `sb_server.create_job` returned 500 for `init`/`sideband` that are not objects; the `sideband` case left an orphan job folder listed in `/jobs` | type check before `folder.mkdir` → HTTP 400 | `tests/test_sb_server.py::check_rejections` (no new folder) |
| 4 | Docs claimed every `sbonest` command behaves like the scripts, but `fit` rejects ONEST configurations that `run.py` runs | README, manual 13 (both languages) and the `sb_cli` docstring say `check|fit|resume` are Sideband-only | – |

Rejected candidates (bug-hunt agent, do not re-investigate without new evidence):
checkpoint key regex vs grouped names (profiles limited to kex/pB/v1n_scale),
`derived_errors` division (kba lower bound 1e-8), glob case sensitivity, server
log-handle leak and fit/resume/archive race (lock + 409), `_log_params` swallowed
ValueError, writer/reader JSON keys, numerical formulas, provenance source list
(sb_parallel/sb_compare/sb_design omitted by design).

## Refactors applied (no behaviour change; golden/worker-identity/output tests unchanged)

- `sb_cli.py`: every command except `serve`/`version` is delegated unchanged to a
  module command line (`sb_run.main` for check/fit/resume with `--check`/`--resume`
  appended, `sb_workflow.main` for report/init-demo/design/compare, `sb_import`,
  `benchmark`); the sbonest parser keeps `add_help=False` stubs so `sbonest --help`
  lists everything and `sbonest fit --help` shows `run.py`'s help. Visible change:
  `sbonest init-demo` prints the workflow's three lines (Config/Data/Fit output prefix)
  instead of one `config:` line.
- `sb_run.py`: `build_parser(description=None)`, `parse_arguments(parser, argv)`,
  `dispatch(parser, args, config, config_dir)` (returns the exit status), `main(argv)`
  returns it; `run.py` uses them for the Sideband branch and `sys.exit(main())`.
  `run.validate_config_shape` is the single shape check (`sb_run.check_config` used
  a copy; its error text is now `Invalid config value types or empty dataset list.`).
  `sbfit.__getattr__` re-exports only `check_config`/`run_config`.
- `sideband.py`: `_check_segments` and `_propagate(a, hrf, seg, n_full, remaining, c0)`
  shared by `profile` and `profile_states` (same operations in the same order;
  golden test and `results/auto_H_refit` comparison pass).
- Single definitions: `sb_diagnostics.number_text`, `sb_report._finite` (with
  `positive`), `sb_report.reject_constant`, `sb_report.SOURCES` (used by
  `sbonest version`), `sb_checkpoint._fsync_directory`; removed
  `sbfit.is_sideband_method` (no callers), unread `WorkerPool.started`, dead
  `getattr(model, "rate_names", ...)` default.
- `tests/_env.py` holds the repository-root `sys.path` insert and the
  single-thread BLAS/MPLBACKEND defaults; every test starts with
  `from _env import ROOT` or `import _env`. `scripts/demo_sideband.py` already puts
  `tests/` on `sys.path` before importing `test_sideband`.
- Not applied on purpose: merging the three `evaluate_many` pool closures (sbfit's
  has a single-vector shortcut; numerical path, little gain) and `sb_compare`'s
  `number()` keeps a local wrapper because callers pass `digits` positionally.

## Verification (local, Python 3.12 venv at `session_artifacts/manual_20261004/tree/.venv`)

`ruff check *.py tests scripts --select F`, `py_compile`,
`scripts/generate_api_reference.py` regenerated (`--check` passes), the full CI test
list including `tests/test_sb_golden.py --out` and `scripts/demo_sideband.py --out`;
results are in the PR description. `tests/test_sb_distribution.py` runs in CI only
(needs `build`/`twine`).

## Still-pending round 3 proposals (unchanged, not started)

See `.archive/HANDOFF.before-restructure-20261005.md` or the PR #17 handoff: 2
runs-test direction + multi-seed field benchmark, 3 `sbonest simulate`, 4 global
optimizer for 3-state, 5 MCMC option, 6 Bruker fixture, 7 release automation +
Docker, 8 type hints/extra Ruff rules, 9 docs site. PR #16
(`execute/sbonest-user-trust-c`, another agent's branch) still needs a rebase.

## Key files & commands

- Environment for every run: `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg`.
- Interpreters: `session_artifacts/manual_20261004/tree/.venv/bin/python` (3.12, no ruff/build/twine); `ruff` is on PATH.
- Local check sequence: `ruff check *.py tests scripts --select F`, `python -m py_compile *.py tests/*.py scripts/*.py`, `python scripts/generate_api_reference.py --check`, the `tests/` list from `.github/workflows/ci.yml`, `python scripts/demo_sideband.py --out NEW_DIR`.
- GitHub: push rights are on the `Gohyang-Matzip` account; `python session_artifacts/improve_20261004/github.py gh <args>` wraps `gh` with that token, push with `... github.py git -c credential.helper= -c 'credential.helper=!gh auth git-credential' push`. `gh pr merge N --merge --match-head-commit <full sha>` after all 10 CI jobs pass, then `git fetch && git merge --ff-only origin/main` on `main`.
- Docs to update with any behaviour change: `docs/manual/SBONEST_MANUAL.md` and `.ko.md`, `README.md`, `docs/guides/SIDEBAND.md`, `AGENTS.md`, `CHANGELOG.md`, `docs/API_REFERENCE.md` via `python scripts/generate_api_reference.py`.
- Archiving rule: never delete artifacts; copy to `.archive/` or a new `results/...` folder.

## Next steps

1. If this round's PR is not merged yet: wait for CI, merge, verify the merge commit's CI, fast-forward `main`.
2. Then the round 3 proposals, when the user picks items.

## Open questions / risks

- The golden tolerances set in PR #17 (stderr 1e-4, proton rates 1e-3 / values 1e-2) were chosen from one Linux run; another runner generation may need a small adjustment.
- `sbonest init-demo` output format changed (three lines); no test or document quoted the old line.
