# Full improvement programme (approved 2026-10-04)

> **For agentic workers:** the user approved every proposal in the 2026-10-04
> improvement list. Work is delivered as four themed PRs, each with executable
> regressions, synchronized bilingual documentation, local verification and
> remote CI before merge. Do not loosen the physical model, optimizer tolerances
> or numerical parity of existing paths.

**Goal:** make SBONEST faster, better diagnosed, scientifically broader and
easier to deploy, without changing existing results.

**Tech stack:** Python 3.12 (CI also 3.13 where dependencies allow),
NumPy/SciPy/Matplotlib, Flask (already a dependency), standard library only
otherwise, executable `test_*.py` regressions.

## PR 1 — `codex/parallel-and-diagnostics` (perf + diagnostics + usability + CI)

- [x] `sb_parallel.py`: spawn-based worker pool holding one `SidebandModel` per
      process; ordered task mapping; identical numerics to the serial path.
- [x] `fit._block_jacobian(evaluate_many=...)`: evaluate perturbed vectors in a
      batch; serial default unchanged.
- [x] `SidebandModel.fit`: Jacobian columns evaluated through the pool when
      `model.pool` is set.
- [x] `fit_multistart`, `profile_likelihood`, `bootstrap_fit`: optional `pool`;
      attempts/points/replicates run in workers, rows saved in index order so the
      checkpoint ordered-prefix rule holds; serial code path shared per item.
- [x] `run.py --workers N`, `sbfit.main --workers N`; checkpoint identity excludes
      workers (results identical).
- [x] Progress/ETA lines for multistart, profile and bootstrap when verbose.
- [x] `run.py --check --identifiability`: scaled singular values, rank,
      condition, expected SE and relative SE at the initial point, weak
      parameters, strongest correlations.
- [x] Report PDF: section-aware pagination, no orphan pages.
- [x] CI: fast `checks` job (3.12 and 3.13 if constraints install) and a
      `workflow` job for the 882-point fit/demo/guide with `--workers 2`.
- [x] `test_sb_parallel.py`: exact parity serial vs parallel for Jacobian, full
      fit, multistart rows, profile rows, bootstrap rows; ordered checkpoints;
      error propagation. `test_sb_check.py`: identifiability output.
- [x] Docs: README, both manuals, both dummy guides, AGENTS, HANDOFF.

## PR 2 — `codex/design-intervals-comparison` (science tools)

- [x] `sb_design.py` + `sb_workflow.py design`: Fisher-information expected SE,
      correlations, rank/condition for scenario JSON (offset grids, RF levels,
      T, sigma) at supplied truth; JSON/TXT/PDF outputs; no optimization.
- [x] `init.profile_interval`: bracket + Brent root finding on the exact
      nuisance-refit profile for kex/pB/v1n_scale at a chi-square threshold;
      every evaluated point retained and checkpointed; warnings for bracket
      failure or a profile below the base minimum.
- [x] `sb_workflow.py compare`: global shared-exchange fit versus per-residue
      individual fits; chi2, AICc, BIC, nested F-test with caveats; parallel
      sub-fits through the pool.
- [x] Tests, docs, HANDOFF.

## PR 3 — `codex/multifield-and-three-state` (model extensions)

- [x] Field groups by `h_larmor_mhz`; `proton_relaxation.mode = "fit"` with
      several fields gives `<peak>.R1H[g]`/`R2H[g]`; optional
      `sideband.nitrogen_relaxation.mode = "per_field"` gives `R1[g]`/`R2a[g]`/`R2b[g]`.
      Single-field names and layouts unchanged.
- [x] `sideband.profile_states`: n-state Liouvillian (16 n × 16 n) with exchange
      matrix; `init.Method = "Sideband_3st_Linear"` and `"Sideband_3st_Triangle"`
      with `h_ppm_c`, `dwC_ppm`, `R2c`; 2-state results unchanged.
- [x] Tests on synthetic two-field and three-state data; docs; HANDOFF.

## PR 4 — `codex/packaging-import-web` (deployment)

- [x] `pyproject.toml` (flat modules), `sbonest` console entry point with
      `check/fit/resume/report/init-demo/design/compare/import-bruker/serve`.
- [x] `sb_import.py`: pure-NumPy Bruker pseudo-2D (`2rr`, `procs`, `proc2s`)
      reader → SBONEST text input with explicit offsets, peaks, reference row
      and noise region; synthetic Bruker fixture test.
- [x] `sb_server.py`: Flask Sideband UI (upload, check, fit with workers,
      progress from checkpoint records, result downloads); live HTTP test.
- [x] Tests, docs, HANDOFF.

## Constraints

- Preserve public names, config keys, units, ordering and tolerances.
- Serial `--workers 1` must reproduce every existing number exactly.
- Archive, never delete; ignore worktrees, environments and run outputs.
- Update both manual languages and both dummy guides for documented behavior.
- Label synthetic evidence; bootstrap/profile intervals remain conditional.

## Execution record

- 2026-10-04: PR 4 worktree `session_artifacts/improve5_20261004/tree`, branch
  `codex/packaging-import-web`, stacked on PR 3.

- 2026-10-04: PR 3 worktree `session_artifacts/improve4_20261004/tree`, branch
  `codex/multifield-and-three-state`, stacked on PR 2 (`e73df58`).

- 2026-10-04: PR 2 worktree `session_artifacts/improve3_20261004/tree`, branch
  `codex/design-intervals-comparison`, stacked on PR 1 (`675109a`).

- 2026-10-04: PR 1 worktree `session_artifacts/improve2_20261004/tree`, branch
  `codex/parallel-and-diagnostics`, base `791fb1dee0f1d8b2e1a0d34fd43b264044c18a73`.


# Round 2 (approved 2026-10-04, after the first programme)

## PR 5 — `codex/perf2-diagnostics-design`

- [x] Worker pool: pre-start every worker; lazy matplotlib import in `estmodel`;
      block-split residual evaluation (`SidebandModel.errFunc`) for the optimizer's
      own evaluations; 882-point fit 10 s with eight workers, bit-identical.
- [x] `sb_diagnostics.py`: reduced chi2 and spread, sigma scale, per-block runs test,
      lag-1 autocorrelation, outliers, chi2/dof-rescaled errors; result JSON, text
      report, regenerated report; `test_sb_diagnostics.py`.
- [x] Optimal design: `optimize` section in design files, backward elimination of
      spectrum rows with Woodbury downdates (kex/pB/parameter/D criteria),
      optimized and uniform comparison scenarios; tests.
- [x] CI: three shards on Python 3.12/3.13/3.14; `CHANGELOG.md`, `CITATION.cff`,
      version 1.2.0.

## PR 6 — `codex/coverage-fields-models`

- [x] `validate_uncertainty.py`: coverage of profile-likelihood intervals and
      bootstrap percentile intervals in the same seeded study.
- [x] 600/800 MHz joint multi-field fit (automatic proton rates, per-field
      nitrogen relaxation) added to `compare_field_models.py` and documented.
- [x] `sb_workflow.py compare --models`: two-state versus three-state with
      AICc/BIC and a multistart-based three-state initialization.

## PR 7 — `codex/import-qa-web-refactor`

- [ ] Bruker import QA PDF, `fq1list` parsing, peak extraction from a reference 2D.
- [ ] Web runner: token authentication, job expiry/archive, analysis options,
      PDF preview.
- [ ] Split `run_config`/validation out of `sbfit.py`; generated API reference;
      release tag v1.2.0.
- Deferred (needs data): experimental 1.2 GHz import → check → design comparison.
