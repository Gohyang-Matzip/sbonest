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

- [ ] `sb_design.py` + `sb_workflow.py design`: Fisher-information expected SE,
      correlations, rank/condition for scenario JSON (offset grids, RF levels,
      T, sigma) at supplied truth; JSON/TXT/PDF outputs; no optimization.
- [ ] `init.profile_interval`: bracket + Brent root finding on the exact
      nuisance-refit profile for kex/pB/v1n_scale at a chi-square threshold;
      every evaluated point retained and checkpointed; warnings for bracket
      failure or a profile below the base minimum.
- [ ] `sb_workflow.py compare`: global shared-exchange fit versus per-residue
      individual fits; chi2, AICc, BIC, nested F-test with caveats; parallel
      sub-fits through the pool.
- [ ] Tests, docs, HANDOFF.

## PR 3 — `codex/multifield-and-three-state` (model extensions)

- [ ] Field groups by `h_larmor_mhz`; `proton_relaxation.mode = "fit"` with
      several fields gives `<peak>.R1H[g]`/`R2H[g]`; optional
      `sideband.nitrogen_relaxation.mode = "per_field"` gives `R1[g]`/`R2a[g]`/`R2b[g]`.
      Single-field names and layouts unchanged.
- [ ] `sideband.profile_states`: n-state Liouvillian (16 n × 16 n) with exchange
      matrix; `init.Method = "Sideband_3st_Linear"` and `"Sideband_3st_Triangle"`
      with `h_ppm_c`, `dwC_ppm`, `R2c`; 2-state results unchanged.
- [ ] Tests on synthetic two-field and three-state data; docs; HANDOFF.

## PR 4 — `codex/packaging-import-web` (deployment)

- [ ] `pyproject.toml` (flat modules), `sbonest` console entry point with
      `check/fit/resume/report/init-demo/design/compare/import-bruker/serve`.
- [ ] `sb_import.py`: pure-NumPy Bruker pseudo-2D (`2rr`, `procs`, `proc2s`)
      reader → SBONEST text input with explicit offsets, peaks, reference row
      and noise region; synthetic Bruker fixture test.
- [ ] `sb_server.py`: Flask Sideband UI (upload, check, fit with workers,
      progress from checkpoint records, result downloads); live HTTP test.
- [ ] Tests, docs, HANDOFF.

## Constraints

- Preserve public names, config keys, units, ordering and tolerances.
- Serial `--workers 1` must reproduce every existing number exactly.
- Archive, never delete; ignore worktrees, environments and run outputs.
- Update both manual languages and both dummy guides for documented behavior.
- Label synthetic evidence; bootstrap/profile intervals remain conditional.

## Execution record

- 2026-10-04: PR 1 worktree `session_artifacts/improve2_20261004/tree`, branch
  `codex/parallel-and-diagnostics`, base `791fb1dee0f1d8b2e1a0d34fd43b264044c18a73`.
