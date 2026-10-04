# HANDOFF: SBONEST improvement programme — PR 2 (design, profile intervals, model comparison)

**Written:** 2026-10-04 (Asia/Seoul)
**Repository:** https://github.com/Gohyang-Matzip/sbonest
**Main checkout:** `/Users/donghanlee/work/projects/sbonest`
**Implementation worktree:** `/Users/donghanlee/work/projects/sbonest/session_artifacts/improve3_20261004/tree`
**Branch:** `codex/design-intervals-comparison` (stacked on `codex/parallel-and-diagnostics`, PR #9)
**Programme plan:** `docs/superpowers/plans/2026-10-04-full-improvements.md`

## Goal and current status

The user approved the complete 2026-10-04 improvement list. PR 1 (worker pool,
identifiability preflight, balanced report pages, CI split) is PR #9. This
handoff covers PR 2: experimental design, likelihood-ratio profile intervals and
shared-versus-individual model comparison. Implementation and local validation
are complete; the branch must be rebased onto main after PR #9 merges, then
pushed, CI-checked and merged. The previous handoff is archived at
`.archive/HANDOFF.before-design-20261004.md`.

## Applied changes (PR 2)

Every item below is **still applied**.

- **[still applied] `sb_analysis.profile_intervals`, `_interval_settings`:**
  `init.profile_interval = {parameters, confidence=0.95, max_evaluations=40,
  relative_tolerance=1e-3, max_doublings=8}` for kex, pB and (scale mode)
  v1n_scale. Each side is bracketed outward in doubling steps from z×local SE,
  then `brentq` locates the crossing of the chi-square(1) quantile on the exact
  nuisance-refit profile. Every evaluated point is retained in evaluation order;
  a target cache avoids re-evaluating bracket ends. Open sides, below-base
  profiles and failed refits are reported, not raised. Completed rows are replayed
  by position and a mismatching target is rejected.
- **[still applied] `profile_point(pool=)`, `_profile_jacobian(evaluate_many=)`:**
  the constrained refit's Jacobian columns can run in the worker pool, which
  speeds the sequential root finding.
- **[still applied] `sbfit.run_config`:** `profile_interval-<name>-N` checkpoint
  records; `info["profile_intervals"]`; validation accepts the new key;
  `check_config` lists the analysis. `sb_report`: validation, summary lines and
  a profile-interval plot with the threshold.
- **[still applied] `sb_design.py` (new), `sb_workflow.py design`:** design JSON
  (base config, `truth` or `truth_result`, scenarios with datasets of v1n_hz, T,
  sigma, absolute or relative offset grids, optional v1err_hz/field_mhz/decoupling).
  Writes noise-free synthetic inputs and `design_config.json` per scenario,
  evaluates `identifiability` at the truth, writes `design.json/.txt/.pdf` with
  points, Σ(points×T), rank/condition, expected and relative SE, derived kex/pB
  SE, weak parameters and strong correlations. Rejects existing outputs and
  malformed designs before writing.
- **[still applied] `sb_compare.py` (new), `sb_workflow.py compare`:** global
  (shared kab/kba) fit plus per-residue individual fits through `run_config`
  (own checkpoints), `_residue_config` pruning of per-residue settings, AICc/BIC
  (Gaussian, absolute sigma), nested F-test, preferred model by AICc, warnings
  for failed sub-fits or a global chi2 below the individual sum; JSON/TXT/PDF.
- **[still applied] Tests:** `test_profile_interval.py` (settings, intervals vs
  1.96·SE within 0.7–1.4, crossings at the threshold, replay without refits,
  mismatch rejection, open intervals, failure recording, run_config/resume/report),
  `test_sb_design.py` (validation, SE ordering, exact sigma scaling, noise-free
  inputs reproduce the truth, pooled parity, truth_result, CLI),
  `test_sb_compare.py` (statistics, config pruning, shared vs distinct synthetic
  exchange → global vs individual preferred, CLI with workers).
- **[still applied] Docs:** manual sections 8.5–8.7 (both languages), README,
  both dummy guides, SIDEBAND.md, AGENTS.md, CI checks job, plan, this handoff.

## Evidence

Bundled example at sigma 0.01 (`session_artifacts/improve3_20261004/tree/session_artifacts/design_try/out_01`):
two RF × 147 offsets → expected kex SE 6.73 s⁻¹; two RF × 60 → 10.76 at 41% of the
saturation time; three RF × 60 → 9.98; one RF (100 Hz) × 147 → 116.5. Small
synthetic interval case: kex 95% [298.20, 302.09] around 300.15 with local SE
0.994 (ratio to 1.96·SE ≈ 1.00), 6 evaluations per quantity. Comparison on
synthetic three-residue data: shared truth → global preferred (F-test p = 0.56),
one residue with doubled kex → individual preferred (p ≈ 5e-40).

Verification: `session_artifacts/improve3_20261004/verification_01/summary.json`.

## Delivery steps — perform only those still missing

1. After PR #9 merges, `git rebase origin/main` in this worktree (expect no
   conflicts: PR 2 adds files and appends to shared ones), rerun the verification.
2. Push, open the PR with `session_artifacts/improve3_20261004/pr2_body.md`,
   require CI for the exact head, merge, verify merge CI, fast-forward main.
3. Continue with PR 3 (`codex/multifield-and-three-state`).

## Scientific boundaries

Profile intervals, design errors and model-comparison statistics are
within-model statements under the supplied absolute sigma and the configured
fixed inputs. They rank descriptions and designs; they do not validate a sample
or prove a mechanism. All evidence is synthetic.
