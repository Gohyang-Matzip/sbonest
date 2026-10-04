# Changelog

All notable changes to SBONEST. Versions follow `pyproject.toml`; dates are Asia/Seoul.

## Unreleased

- Fix `sb_workflow.py compare`: residues present in the data but absent from
  `config.residues` are active (as documented) and are now switched off in the
  per-residue fits and given state-C shifts in `--models` comparisons; previously
  the "individual" fits silently kept them, invalidating AICc/BIC/F statistics, and
  three-state derivation failed on the missing `h_ppm_c`.
- Fix `sb_workflow.py design` with `nitrogen_relaxation.mode = "per_field"` and
  several field groups: grouped truth names (`A1.R2a[1]`) are accepted; the
  design previously always raised "truth must include A1.R2a".
- The web runner rejects configurations whose `init` or `sideband` sections are
  not objects with HTTP 400 instead of a 500 that left an orphan job folder.
- `sbonest check|fit|resume` are documented as Sideband-only (ONEST models run
  through `run.py`).
- Refactor without behaviour change: `sbonest` delegates every command to the
  module command lines (`sb_run.main`, `sb_workflow.main`), `run.py` shares
  `sb_run.build_parser`/`dispatch`, `sideband.profile` and `profile_states` share
  one propagation loop, duplicated helpers (`number_text`, `_finite`,
  `reject_constant`, `validate_config_shape`, `_fsync_directory`, `SOURCES`) are
  single definitions, the unused `sbfit.is_sideband_method` is removed, and the
  tests share `tests/_env.py` for the root path and thread settings. Golden,
  worker-identity and output tests are unchanged.
- Reorganize the repository: tests move to `tests/`, research and maintenance
  scripts to `scripts/`, and documents to `docs/manual/`, `docs/guides/`,
  `docs/reports/` and `docs/images/` (`MANUAL.md` becomes
  `docs/manual/ONEST_MANUAL.md`, `result.md` becomes
  `docs/reports/ONEST_CLUSTER_RESULT.ko.md`). Installed modules stay flat at the
  root; moved tests and scripts insert the root into `sys.path`, and CI, README,
  manuals and the dummy guides use the new paths. No numerical behaviour changes.
- Golden regression (`tests/test_sb_golden.py`): tolerances now reflect what the
  model determines. Linux CI reproduces the archived macOS chi2 to 1e-11 and every
  non-proton value to 1e-6, but standard errors differ by up to 1e-4 and the weakly
  constrained peakwise proton rates by up to 2e-3 (G2.R1H). Standard errors use
  rtol 1e-4 (proton rates 1e-3) and proton-rate values rtol 1e-2; chi2 (1e-7) and
  all other values (1e-6) keep their tolerances.
- Package runtime diagnostics and portable synthetic two-RF demo resources in
  wheel and source distributions; version reporting no longer loads demo data.
- Align installed `compare --models/--h-ppm-c` and `serve --token/--max-age-days`
  with their script entry points.
- Keep browser action failures visible, preserve selection after failed archive,
  ignore obsolete job-status responses, and support previews with or without tokens.
- Add independent runs and correlation direction labels without changing the
  two-sided tests, thresholds, default absolute-sigma errors or numerical model.
  Residual warnings call for inspection of model, noise and acquisition rather
  than identifying a cause; concatenated overall runs are descriptive.

## 1.2.0 — 2026-10-04

- Worker pool pre-warms every process and predicts data blocks in parallel for the
  optimizer's own residual evaluations; matplotlib is imported lazily. The bundled
  882-point fit takes about 10 s with eight workers (31 s serially), bit-identical.
- Residual diagnostics (`sb_diagnostics.py`): reduced chi-square with its expected
  spread, residual-based sigma scale, per-block reduced chi-square, runs tests,
  lag-1 autocorrelation, outlier counts and chi2/dof-rescaled standard errors
  (`stderr_rescaled`), in the result JSON, text report and regenerated reports.
- Optimal design (`optimize` in a design file): backward elimination of whole
  spectrum rows on the Fisher information for kex, pB, D or a parameter criterion,
  with uniform-budget comparison.
- CI runs three test shards on Python 3.12, 3.13 and 3.14; CHANGELOG and CITATION added.
- `validate_uncertainty.py --profile-interval ... --inner-bootstrap B --workers N`:
  coverage of profile-likelihood and bootstrap percentile intervals in the same
  seeded study.
- `compare_joint_fields.py` and `results/field_comparison_joint_20261004/`: joint
  600/800 MHz fits with field-group proton/nitrogen relaxation; proton rates are
  not determined at those fields and need bounds, kex/pB unchanged.
- `sb_workflow.py compare --models`: two-state versus three-state Sideband models
  with multistart-based three-state initialization, AICc/BIC and degeneracy warnings.
- `sb_run.py` holds `check_config`/`run_config` (still importable from `sbfit`);
  every public function has a docstring and `docs/API_REFERENCE.md` is generated
  by `generate_api_reference.py` (checked in CI).
- Bruker import: Bruker frequency lists (`fq1list`) read directly, `--peaks-from-reference`,
  `--qa-pdf` figure with reference row, windows, noise region and extracted profiles.
- Web runner: optional analyses in the upload form, PNG previews, job archiving with
  `--max-age-days`, token access with `SBONEST_TOKEN`/`--token`.

## 1.1.0 — 2026-10-04

- `--workers N` parallel execution (PR #9), `--check --identifiability`, balanced
  report pages, CI on Python 3.12 and 3.13, seeded Monte Carlo smoke test.
- Likelihood-ratio profile intervals (`init.profile_interval`), experimental design
  (`sb_workflow.py design`) and shared-versus-individual model comparison
  (`sb_workflow.py compare`) (PR #10).
- Field groups with field-specific proton and nitrogen relaxation and three-state
  Sideband models (PR #11).
- Installable `sbonest` command, Bruker pseudo-2D import and the Sideband web
  runner (PR #12).

## 1.0.0 — 2026-10-02 to 2026-10-04

- SBONEST sideband fitting with automatic peakwise proton relaxation (PR #1), CEST
  manuscript materials with full fitting uncertainties (PR #2), 600/800 MHz
  comparison (PR #3, #4), stricter spectrum loader (PR #5), diagnostics and
  uncertainty analysis with checkpoints, resume, bootstrap and reports (PR #6–#8).
