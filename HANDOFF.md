# HANDOFF: SBONEST improvement programme — PR 3 (field groups and three-state exchange)

**Written:** 2026-10-04 (Asia/Seoul)
**Repository:** https://github.com/Gohyang-Matzip/sbonest
**Main checkout:** `/Users/donghanlee/work/projects/sbonest`
**Implementation worktree:** `/Users/donghanlee/work/projects/sbonest/session_artifacts/improve4_20261004/tree`
**Branch:** `codex/multifield-and-three-state` (stacked on `codex/design-intervals-comparison`, PR #10)
**Programme plan:** `docs/superpowers/plans/2026-10-04-full-improvements.md`

## Goal and current status

The user approved the complete 2026-10-04 improvement list. PR 1 (worker pool,
identifiability, pagination, CI) merged as PR #9 (`ec8f81a`); PR 2 (design,
profile intervals, model comparison) is PR #10. This handoff covers PR 3: field
groups with field-specific relaxation and three-state Sideband models.
Implementation and local validation are complete; after PR #10 merges this
branch is rebased onto main, pushed, CI-checked and merged. The previous handoff
is archived at `.archive/HANDOFF.before-models-20261004.md`.

## Applied changes (PR 3)

Every item below is **still applied**.

- **[still applied] `sideband.profile_states`, `stationary_populations`:** n-site
  NH Liouvillian (16 n dimensions) with an exchange rate matrix, equilibrium
  starting populations from the stationary distribution and the site-A Nz
  observable. Two sites reproduce `profile` to 2e-16; a three-site model with an
  unpopulated third site reproduces the two-site curve to 4e-9. `profile` itself
  is untouched, so existing two-state numerics are unchanged.
- **[still applied] `SidebandModel` layout generalization:** `SIDEBAND_METHODS`
  maps `Sideband`, `Sideband_3st_Linear`, `Sideband_3st_Triangle` to the ONEST
  parameter layouts and site counts; `rate_names`, `shift_keys`, `n_relax_keys`,
  `h_keys`, `n_global`, `physical_bounds`, `exchange_matrix`. Datasets sharing
  `h_larmor_mhz` form field groups (`field_groups`, `dataset_group`); automatic
  proton rates become `R1H[g]`/`R2H[g]` per group and
  `sideband.nitrogen_relaxation.mode = "per_field"` gives `R1[g]`/`R2a[g]`/`R2b[g]`
  (and `R2c[g]`). Suffixes appear only with several groups, so single-field names
  and the parameter order are unchanged. Ungrouped names in `initial`, `bounds`,
  `vary` and multistart starts are aliases for every group (`expand_name`,
  `normalize_init`, `init_config`); `model.config` stays the caller's live object
  and aliases are expanded at every read.
- **[still applied] Three-state fitting:** 48-dimensional propagation in `calc`,
  `h_ppm_c` required, `dwC_ppm`/`R2c` per residue, return rates `kba`/`kcb`
  bounded below by 1e-8, ONEST's default three-state start (no rate grid).
  `diagnostics` adds `model`, `states`, `field_groups`,
  `nitrogen_relaxation_mode` and `exchange` (rates, stationary populations,
  pairwise kex sums); `kex`/`pB`/`derived_se` are null for three states.
  Profiles and profile intervals of kex/pB are rejected for three-state models;
  bootstrap reports percentiles for every rate and omits kex/pB.
- **[still applied] Shared code:** `run.load_config`/`run.main`, `benchmark.py`,
  `sb_workflow.py compare`, `sb_design.py` accept any Sideband method;
  `_profile_jacobian`, `local_jacobian`, `identifiability`, `sb_bootstrap`,
  `_multistart_settings`, `_read_snapshot`, `base_fit_config` use
  `rate_names`/`n_global`/`init_config`; `check_config` reports `model`,
  `states`, `field_groups`; text reports list per-group rates and groups.
- **[still applied] `test_sb_models.py` (new):** two fields (600/800 MHz) with
  field-specific true H and N relaxation recovered by a joint fit with
  `per_field` nitrogen relaxation; expanded-initial layout; alias handling in
  initial/multistart; shared-mode names; linear three-state synthetic truth
  recovered (rates within 25% + 3 s⁻¹, populations within 0.02–0.03, shifts
  within 0.2 ppm); null kex/pB; report; bootstrap; rejection of kex/pB
  profiles/intervals, missing `h_ppm_c`, unknown methods; triangle rate matrix;
  `--check --identifiability` for the triangle model. `test_sideband.py` now
  asserts that two fields with automatic proton rates produce grouped names.
- **[still applied] Docs:** manual section 12 (both languages) and section 9
  wording, README, SIDEBAND.md, AGENTS.md, CI, plan, this handoff.

## Findings

- A first version stored a normalized copy of the configuration on the model;
  `test_sb_output.check_partial_fit_parity` and `test_sb_check` mutate the caller's
  `cfg["init"]` after construction and expect the model to see it. The model now
  keeps the live object and expands aliases on every read.

## Evidence

Bundled 882-point two-state example after the refactor: chi2 802.8203212916861,
every parameter value and standard error identical to
`results/auto_H_refit/fits/two_RF_result.json`
(`session_artifacts/improve4_20261004/tree/session_artifacts/parity_full`).
Verification: `session_artifacts/improve4_20261004/verification_02/summary.json`
(ruff, compile, 19 regression scripts, demo, `git diff --check`, all exit 0;
`verification_01_failed` records the normalized-copy regression).

## Delivery steps — perform only those still missing

1. After PR #10 merges, rebase this branch onto main, rerun the verification.
2. Push, open the PR with `session_artifacts/improve4_20261004/pr3_body.md`,
   require CI for the exact head, merge, verify merge CI, fast-forward main.
3. Continue with PR 4 (`codex/packaging-import-web`).

## Scientific boundaries

Field-specific relaxation and three-state exchange add parameters that the
data may not determine; identifiability diagnostics, restarts and comparisons
with distinct output prefixes are required before interpretation. The 600/800
MHz benchmark was not rerun. All evidence is synthetic.
