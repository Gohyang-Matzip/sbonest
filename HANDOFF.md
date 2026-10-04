# HANDOFF: SBONEST loader debugging and refactor

**Written:** 2026-10-04 (Asia/Seoul)
**Working directory:** `/Users/donghanlee/work/projects/sbonest`
**Repository:** https://github.com/Gohyang-Matzip/sbonest
**Implementation branch:** `codex/refactor-spectrum-loader`
**Target branch:** `main`
**Starting commit:** `010ad74299d09ec1b97de95674849866b0340f31`

## User requests and completion criteria

The user requested, in order: `debugging this repo`, `refactor`, `handoff`,
and `commit push PR merge`. The final request authorizes committing these
changes, pushing the branch, opening a PR, merging after passing CI, verifying
the merged commit's CI, and synchronizing local `main`.

The code work is complete and locally validated. No further feature or research
work is pending. Check live GitHub state before repeating delivery actions; this
file is part of the change being published, not proof that publication finished.

## Changes in this delivery

Every change below is **still applied**; no production fix was reverted.

- **[still applied] `est_data.py`:** reject missing fourth column-header lines,
  observations outside residue blocks, empty files, and empty residue blocks.
  Previously a standalone comment inside a spectrum could silently discard all
  remaining points, and a missing column-header line could swallow the first
  residue. Invalid inputs now raise `ValueError` with the filename and a useful
  explanation. A `#` line still terminates a data block; comments between complete
  residue blocks and blank lines remain accepted.
- **[still applied] `est_data.py` refactor:** module-level regexes share one numeric
  pattern; `_read_conditions()` parses acquisition values; `_read_points()` fills
  one spectrum and returns the next comment/header. `addData()` retains its public
  signature and handles residue format validation and attachment. File length
  decreased from 293 lines after the bug fix to 213. Numerical fitting code,
  accepted numeric grammar, error messages, and random draw order are preserved.
- **[still applied] `test_debugging.py`:** `check_loader_structure()` covers eight
  malformed input cases plus valid data. `check_loader_formats_and_noise()` covers
  full/simple headers, multi-file residue ordering, defaults, inline comments,
  mixed-format rejection, and seeded RF/intensity perturbations. Both are called
  by the executable test runner.
- **[still applied] `SBONEST_MANUAL.md` and `SBONEST_MANUAL.ko.md`:** matching input
  rules and troubleshooting entries for the stricter loader.
- **[still applied] `HANDOFF.md`:** this current briefing. Its previous complete
  version is archived locally at `.archive/handoff-before-loader-20261004-094051.md`.

The initial worktree was clean. Only the five files listed above belong in this
PR. Temporary environments, generated fits, verification logs, and backups are
ignored local artifacts, not files to add to Git.

## Verified results

The eight checks in `session_artifacts/refactor_20261004/verification/summary.json`
all exited 0 under Python 3.12.14:

1. `ruff check *.py --select F`
2. `python -m py_compile *.py`
3. `python test_performance.py`
4. `python test_debugging.py`
5. `python verify_3state.py`
6. `python test_sideband.py`
7. `python demo_sideband.py --out <fresh directory>`
8. The real `run.py` CLI on a copy of `example/sideband_auto_H/two_RF.json` with
   absolute input paths and a fresh output prefix, using `--no-pdf`.

`session_artifacts/refactor_20261004/parser_parity.json` records 272 exact
comparisons (68 tracked spectrum files, each with all four noise-flag
combinations), plus 64 boundary/error comparisons against the pre-refactor
loader. Parsed state, exceptions, partial state after failure, and subsequent
random draws matched. Independent review found no actionable regression.

`session_artifacts/refactor_20261004/fit_parity.json` records exact equality
of all result fields except configuration paths for five demo fits and the
automatic-H fit, compared with `session_artifacts/debug_20261004/`.
The automatic-H fit has 882 observations, 24 parameters, rank 24, 858 degrees
of freedom, and chi2 = 802.8203212916861. This is synthetic evidence.

Before the initial fix, the new structure test failed on all eight malformed
inputs. After the fix, the ONEST CLI, MC CLI, and both Sideband entry points
rejected the interrupted-block reproduction without tracebacks or fit outputs.
Their local logs are in `session_artifacts/debug_20261004/reproductions/`.

## Environment and reproducible commands

Use `session_artifacts/debug_20261004/.venv/bin/python` for Python 3.12 locally.
The root `.venv/bin/python` is Python 3.14.7 and `.venv/bin/ruff` was absent;
that original environment is preserved. Do not mistake it for the CI environment.
A fresh checkout should create a Python 3.12 environment and install
`requirements-sideband.txt` plus `ruff`. CI uses Python 3.12 and PyPI
`optimalcontrol-nmr==0.5.0` was used for these local checks.

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg
session_artifacts/debug_20261004/.venv/bin/ruff check *.py --select F
session_artifacts/debug_20261004/.venv/bin/python -m py_compile *.py
session_artifacts/debug_20261004/.venv/bin/python test_performance.py
session_artifacts/debug_20261004/.venv/bin/python test_debugging.py
session_artifacts/debug_20261004/.venv/bin/python verify_3state.py
session_artifacts/debug_20261004/.venv/bin/python test_sideband.py
```

The pre-refactor sources and diff are preserved in
`session_artifacts/refactor_20261004/before_est_data.py`, `before_test_debugging.py`,
and `before.diff`. Final local demo outputs are in
`session_artifacts/refactor_20261004/demo/`; automatic-H results are
`session_artifacts/refactor_20261004/auto_H_two_RF_result.{json,txt}`.
Never rerun a saved config with its old `Project Name`; select a fresh output
prefix. Dataset paths resolve from the config directory; output prefixes resolve
from the process working directory.

## Delivery: check live state, then resume only missing steps

```bash
git status --short --branch
git fetch origin
gh pr list --repo Gohyang-Matzip/sbonest --state all \
  --head codex/refactor-spectrum-loader \
  --json number,state,url,headRefOid,mergeCommit
gh run list --repo Gohyang-Matzip/sbonest --branch main --limit 5
```

At delivery preparation, local `main` equaled `origin/main` at the starting
commit above, with no open PR. The branch was created from that commit.
GitHub account `Gohyang-Matzip` was verified to have push/admin access;
the globally active account was `dleess`. Use scoped authentication for
`Gohyang-Matzip`, obtained through `gh auth token --user Gohyang-Matzip` into a
subprocess environment. Never print or save the token, and do not change global
credentials. Scoped Git pushes can use `credential.helper=!gh auth git-credential`.

1. If this branch has an open PR, inspect its exact head and require passing CI
   before merging. Do not recreate an existing PR.
2. If merged, verify the push CI for the merge SHA, then fast-forward local `main`
   to `origin/main`. Report the PR URL, merge SHA, and CI outcome.
3. Preserve any unrelated work; do not reset or discard dirty files.
4. When delivery is already complete, stop. This handoff does not authorize new
   features, new experiments, or rerunning historical research studies.

## Constraints and boundaries

Follow `AGENTS.md`. Never destructively delete artifacts; archive them or ignore
large raw data. Preserve public names, parameter order, units, tolerances, and
calculation-time source hashes. The 600/800 MHz study uses fixed H rates and
separate field fits; the current examples fit H rates automatically. Synthetic
fits and local covariance are not experimental validation.

No behavior changed in `sideband.py`, `sbfit.py`, `fit.py`, or `estmodel.py`.
No dependencies were added. No manuscript, historical result, or figure was
modified. The longer field-comparison study was not rerun during this task.
Remote PR/merge CI is live state and must be verified as described above.
