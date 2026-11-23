# HANDOFF: SBONEST 1.3.0 and the round 3 proposals

**Written:** 2026-11-23 (Asia/Seoul) · **Working dir:** `/Users/donghanlee/project/sbonest` · **Branch:** `main` after PRs #19 to #22
**Repository:** https://github.com/Gohyang-Matzip/sbonest · **Last release:** v1.3.0 · **PRs this round:** #19, #20, #21, #22
**Previous handoffs:** `.archive/HANDOFF.before-deep-round-20261123.md` (debug/refactor round, PR #18). `.archive/` is gitignored.

## State

PR #19 (`codex/deep-round`) delivered four pieces of work, all in the 1.3.0
section of `CHANGELOG.md`:

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

## Delivery

Everything above was merged as PR #19 (merge commit `f1d3dcc`, CI green on the PR
and on `main`) and released as v1.3.0. Push rights: the `Gohyang-Matzip` `gh`
account (`gh auth switch -u Gohyang-Matzip`; its token needs the `workflow` scope
to push workflow files, `gh auth refresh -s workflow`). Push with
`git -c credential.helper= -c 'credential.helper=!gh auth git-credential' push`,
because the keychain holds another account. Merge with
`gh pr merge N --merge --match-head-commit <full sha>` after the ten CI jobs
pass, verify the merge commit's CI, then `git merge --ff-only origin/main`.

## Verification (local, Python 3.12, `.venv` in the repository root)

`ruff check *.py tests scripts --select F`, `py_compile`,
`scripts/generate_api_reference.py --check`, and the test list of
`.github/workflows/ci.yml` (23 tests) plus `tests/test_sb_golden.py`,
`tests/test_sb_distribution.py` and `scripts/demo_sideband.py` all pass on the
branch head; CI passed on PR #19 and on its merge commit.

## After v1.3.0 (all under "Unreleased" in `CHANGELOG.md`)

- PR #20: release 1.3.0; tag `v1.3.0` ran `.github/workflows/release.yml` for the
  first time (wheel and sdist on the GitHub release).
- PR #21: `init.mcmc` (`sb_mcmc.py`, manual 8.9) and `init.multistart.global`
  (`sb_analysis.global_start`, manual 8.1). Checks: `tests/test_sb_mcmc.py`,
  `tests/test_sb_analysis.py::check_global_start`, pooled identity in
  `tests/test_sb_parallel.py`.
- PR #22: `Dockerfile` with a CI `container` job and image publishing in the
  release workflow (first publication happens with the next tag; no Docker on
  the development machine, CI is the only test), the MkDocs site
  (`mkdocs.yml`, `.github/workflows/docs.yml`, GitHub Pages enabled with
  workflow builds, https://gohyang-matzip.github.io/sbonest/), and
  `scripts/compare_field_seeds.py` with `results/field_seeds_20261123/`.

## Round 3 proposals: status

Done: 2 (direction labels earlier, multi-seed benchmark in PR #22), 3 simulate,
4 global search, 5 MCMC, 7 release automation and Docker, 9 documentation site.
Dropped: 8 further Ruff rules (the bug-prone families B, ARG, RET, BLE, PLE, PLW
report 37 findings in the root modules, none a defect, and enforcing them would
churn hashed sources).
Not possible here: 6 Bruker fixture from real data. No measured Bruker dataset
exists in this checkout or on this machine (`**/ser`, `**/pdata/` are
gitignored); `tests/test_sb_import.py` still uses the synthetic `write_bruker`.
It needs a small processed pseudo-2D (`pdata/1` with `procs`, `proc2s`, `2rr`)
and its offset list from the user.
PR #16 (`execute/sbonest-user-trust-c`, another agent's branch) still needs a rebase.

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
- The 100-seed study shows small SBONEST biases at 800 MHz in pB (+0.00010 ±
  0.00003) and RF scale (−0.0014 ± 0.0004), about three bias standard errors;
  unexplained (nonlinear-estimator bias, the fixed H rates or chance over six
  tested quantities are all possible). kex is unbiased there.
- `init.mcmc` uses flat priors inside the bounds and needs a positive-definite
  local covariance; with the default automatic H rates (weak, often at bounds)
  users must fix those parameters first.
- The MkDocs site is public; the next release tag publishes the first container image.
- Bug hunting by reading has reached low yield after three rounds; the importer
  and web runner boundaries were the last places with findings.
