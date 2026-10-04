# Changelog

All notable changes to SBONEST. Versions follow `pyproject.toml`; dates are Asia/Seoul.

## Unreleased

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
