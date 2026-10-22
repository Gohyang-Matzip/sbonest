# manuscript/sideband_30ppm/

Published Sideband CEST manuscript (1.2 GHz), its tracked uncertainty revision,
and the numerical evidence behind Tables 3 and 5-7. Earned its own file: distinct
domain, allowlist tracking, and author/hash rules unlike the rest of the repo.

## Where to look

| Need | Location |
|------|----------|
| What each artifact is | `README.md` (Evidence table) |
| Executed research scripts | `compare_conventional.py`, `check_carbon_configuration.py`, `collect_fit_errors.py`, `revise_fit_errors.py` |
| Source hashes of those scripts | `results/*/summary.json` (`source_sha256`) |
| Revision ledger and table verification | `qa/fit_errors_20261003/` |

## Rules

- Only allowlisted paths are tracked (`.gitignore` `/manuscript/*` block). Adding a
  file here needs a matching `!` line, or it stays local.
- The four scripts are hashed executed sources: never edit or reformat them; a
  needed change becomes a new script with a new output directory.
- Their output paths are fixed and refuse to overwrite. Move the old output to an
  archival path in a separate checkout before regenerating.
- Run them from the repository root with `tests/` on the import path
  (`PYTHONPATH=tests`; `check_carbon_configuration.py` imports `test_sideband`).
- Apply DOCX tracked changes with author `Donghan Lee`.
- Uncertainties are within-model local standard errors under the supplied absolute
  sigma (no reduced chi-square rescaling); keep that labeling. All data is synthetic.
