# HANDOFF: SBONEST deep debugging round, input-boundary fixes and `simulate`

**Written:** 2026-11-23 (Asia/Seoul) · **Working dir:** `/Users/donghanlee/project/sbonest` · **Branch:** `codex/deep-round` (from `main` at `14277df`, the merge of PR #18)
**Repository:** https://github.com/Gohyang-Matzip/sbonest · **Last release:** v1.2.0 · **This round's PR:** not opened yet (see "Blocked")
**Previous handoffs:** `.archive/HANDOFF.before-deep-round-20261123.md` (debug/refactor round, PR #18). `.archive/` is gitignored.

## State

`codex/deep-round` is ahead of `main` and has never been pushed. It holds three
pieces of work, all recorded under "Unreleased" in `CHANGELOG.md`:

1. The deep debugging round of 2026-10-22 (ten commits; plan, reproduction scripts
   and the prepared PR text are in `.archive/deep-round-20261022/`, `pr_body.md`).
2. Input-boundary fixes (2026-11-23), each with a regression check:
   - `sb_import.convert` refuses an existing `--qa-pdf` target before writing the
     dataset (the dataset used to be left behind and blocked the rerun);
     repeated `--peak` labels are an error; non-finite offsets are rejected
     (`tests/test_sb_import.py`).
   - `sb_server.RESERVED_NAMES` also covers the job's fit outputs
     (`fit_result.txt`, `fit_predictions.csv`, `fit.pdf`, `fit_data.pdf`,
     `fit_checkpoint`, `preview.png`); such an upload gets HTTP 400 before a job
     folder exists (`tests/test_sb_server.py::check_upload_validation`).
3. Round 3 proposal 3, `sbonest simulate` (2026-11-23):
   `sb_workflow.simulate(config, config_dir, out, *, seed, truth=None)` and the
   `simulate` subcommand of `sb_workflow.py`/`sbonest`. It evaluates the model at
   the configuration's starting point (`prepare_fit`), applies `truth` overrides
   (parameter names, ungrouped aliases, or `kex` with `pB` for two-state models),
   draws noise with `sb_bootstrap.bootstrap_draws` (seed, index 0) and writes
   `data/dataset_N.txt`, `fit.json` and `truth.json` into a new directory.
   Documented in manual 8.2 and 13 (both languages), README, `AGENTS.md`;
   checked by `tests/test_sb_workflow.py::check_simulate` (no fitting).

4. Round 3 proposal 7, release automation (2026-11-23, without Docker):
   `.github/workflows/release.yml` runs on a `v*` tag, refuses a tag that
   disagrees with `pyproject.toml`/`CITATION.cff` or has no changelog section,
   builds, runs `twine check --strict` and creates the GitHub release with the
   distributions and that section as notes. Its version step was run locally
   for a matching and a mismatching tag; the workflow itself has never run on
   GitHub. `tests/test_sb_cli.py::check_packaging_metadata` checks that the
   three version records agree.

## Blocked

Pushing and merging need write access to `Gohyang-Matzip/sbonest`. The `gh`
accounts on this machine (`dleess`, `deepnmr`) have read access only, and the
helper named in the previous handoff
(`session_artifacts/improve_20261004/github.py`) does not exist in this checkout.
Once an account with write access is logged in:
`git push -u origin codex/deep-round`, open the PR (body: the archived
`pr_body.md` plus items 2 and 3 above), wait for the ten CI jobs,
`gh pr merge N --merge --match-head-commit <full sha>`, verify the merge commit's
CI, then `git fetch && git merge --ff-only origin/main` on `main`.

## Verification (local, Python 3.12, `.venv` in the repository root)

`ruff check *.py tests scripts --select F`, `py_compile`,
`scripts/generate_api_reference.py --check`, and the test list of
`.github/workflows/ci.yml` (23 tests) plus `tests/test_sb_golden.py`,
`tests/test_sb_distribution.py` and `scripts/demo_sideband.py` all pass on the
branch head.

## Remaining round 3 proposals (not started)

2 multi-seed field benchmark, 4 global optimizer for three-state fits, 5 MCMC
option, 6 Bruker fixture from real data, 7 Docker image (release automation is
done), 9 documentation site. Proposal 8 (further Ruff rules) was surveyed and
dropped: the bug-prone families (B, ARG, RET, BLE, PLE, PLW) report 37 findings
in the root modules, none a defect (loop closures called inside their iteration,
deliberate broad excepts around worker and optimizer failures), and enforcing
them would churn hashed sources. PR #16
(`execute/sbonest-user-trust-c`, another agent's branch) still needs a rebase.
The original proposal text was in a handoff of the earlier checkout and is not in
this one.

## Key files & commands

- Environment: `uv venv -p 3.12 .venv && uv pip install -p .venv/bin/python -r requirements-sideband.txt -c constraints-sideband.txt ruff build twine` (`.venv/` is gitignored).
- For every run: `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg`.
- Local check sequence: `ruff check *.py tests scripts --select F`, `python -m py_compile *.py tests/*.py scripts/*.py`, `python scripts/generate_api_reference.py --check`, the `tests/` list from `.github/workflows/ci.yml`, `python scripts/demo_sideband.py --out NEW_DIR`.
- Docs to update with any behaviour change: `docs/manual/SBONEST_MANUAL.md` and `.ko.md`, `README.md`, `docs/guides/SIDEBAND.md`, `AGENTS.md`, `CHANGELOG.md`, `docs/API_REFERENCE.md` via `python scripts/generate_api_reference.py`.
- Archiving rule: never delete artifacts; copy to `.archive/` or a new `results/...` folder.

## Open questions / risks

- The golden tolerances set in PR #17 (stderr 1e-4, proton rates 1e-3 / values 1e-2) were chosen from one Linux run; another runner generation may need a small adjustment.
- `simulate` uses the configuration's starting point as the default truth. With a
  `kex`/`pB` grid that point depends on the template data, so pass `--truth` when
  the generating values matter.
- Bug hunting by reading has reached low yield after three rounds; the importer
  and web runner boundaries were the last places with findings.
