# HANDOFF: SBONEST improvement programme — PR 4 (packaging, Bruker import, web runner)

**Written:** 2026-10-04 (Asia/Seoul)
**Repository:** https://github.com/Gohyang-Matzip/sbonest
**Main checkout:** `/Users/donghanlee/work/projects/sbonest`
**Implementation worktree:** `/Users/donghanlee/work/projects/sbonest/session_artifacts/improve5_20261004/tree`
**Branch:** `codex/packaging-import-web` (stacked on `codex/multifield-and-three-state`, PR 3)
**Programme plan:** `docs/superpowers/plans/2026-10-04-full-improvements.md`

## Goal and current status

The user approved the complete 2026-10-04 improvement list. PR 1 merged as #9;
PR 2 is #10; PR 3 is committed on its branch. This handoff covers PR 4, the last
of the programme: an installable `sbonest` command, a Bruker pseudo-2D
importer and a Sideband web runner. Implementation and local validation are
complete; after PR 3 merges this branch is rebased onto main, pushed, CI-checked
and merged. The previous handoff is archived at
`.archive/HANDOFF.before-packaging-20261004.md`.

## Applied changes (PR 4)

Every item below is **still applied**.

- **[still applied] `pyproject.toml`, `sb_cli.py`:** setuptools build with the
  root modules listed as `py-modules` (flat layout preserved so provenance source
  hashes keep their meaning), dependencies matching the requirement files and
  the console script `sbonest = sb_cli:main`. Subcommands `check`, `fit`,
  `resume`, `report`, `init-demo`, `design`, `compare`, `import-bruker`
  (passthrough to `sb_import.main`), `serve`, `benchmark` (passthrough to the
  new `benchmark.main`) and `version` (package version plus executed-source
  hashes) call the same functions as the scripts. Verified with an editable
  install into a Python 3.13 environment and `sbonest version`/`init-demo` from
  another directory.
- **[still applied] `sb_import.py`:** NumPy-only reader for Bruker processed
  pseudo-2D data (`procs`/`proc2s` JCAMP parameters incl. list values, `2rr`
  submatrix layout, both byte orders, int32/float64, `NC_proc` scaling, proton
  ppm axis from `OFFSET`/`SW_p`/`SF`). `convert` takes explicit offsets (ppm or
  Hz with carrier), peaks with optional starting dw, reference row, noise region,
  saturation time and RF, window statistic (max/sum), header R2a/R2b and excluded
  rows; writes the SBONEST text format and a JSON summary; refuses existing
  outputs and invalid inputs before writing.
- **[still applied] `sb_server.py`:** Flask runner. `POST /jobs` stores the
  configuration (datasets rewritten to uploaded names, prefix `fit`), data and
  optional waveform under `SB_JOBS/<id>/` (`SBONEST_JOBS_DIR`), runs
  `run.py --check --identifiability` and returns the JSON; `POST /jobs/<id>/fit`
  starts a background `run.py` with the chosen workers; `/resume` continues a
  checkpointed job; `/report` regenerates reports with a fresh prefix;
  `GET /jobs/<id>` reports process state, checkpoint records, result summary and
  log tail; downloads are restricted to files inside the job directory. The HTML
  page polls the status. No authentication; local trusted use.
- **[still applied] `benchmark.main(argv)`** extracted from the script guard;
  `sb_design` restricts scenario configurations to active residues (a base
  configuration with switched-off residues previously failed with "Residue not
  found").
- **[still applied] Tests:** `test_sb_import.py` (synthetic Bruker directories in
  three storage variants, parameter parsing, Hz offsets, conversion reloaded
  through `est_data` with recovered profiles and errors, rejections, CLI),
  `test_sb_server.py` (test client: rejections, upload and preflight, background
  fit to completion, double-start conflict, resume, report, PDF download, path
  traversal rejected, job listing), `test_sb_cli.py` (pyproject metadata and
  module list, every subcommand on a small synthetic case, passthroughs, errors).
- **[still applied] Docs/CI:** manual sections 13–15 (both languages), README,
  SIDEBAND.md, AGENTS.md, `.gitignore` (`SB_JOBS/`, egg-info, build), CI checks
  job runs the three new scripts and the workflow job installs the package and
  runs `sbonest version`/`check`.

## Evidence

Verification: `session_artifacts/improve5_20261004/verification_01/summary.json`
(ruff, compile, 22 regression scripts, demo, `git diff --check`). Editable
install log: `session_artifacts/improve2_20261004/venv313` contains the
`sbonest` console script.

## Delivery steps — perform only those still missing

1. After PR 3 merges, rebase onto main, rerun the verification.
2. Push, open the PR with `session_artifacts/improve5_20261004/pr4_body.md`,
   require CI for the exact head, merge, verify merge CI, fast-forward main.
3. The programme is then complete; update the plan's execution record.

## Scientific boundaries

The importer does not judge phasing, baseline, overlap or the reference choice;
converted profiles must be inspected. The web runner executes the same code as
the command line and adds no scientific claims.
