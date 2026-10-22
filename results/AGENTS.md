# results/

Archived numerical evidence. Treat every file as read-only: these outputs are
cited by documents, golden tests, and recorded source hashes. Earned its own
file: distinct domain with hash and archiving rules.

## Where to look

| Directory | Contents | Consumer |
|-----------|----------|----------|
| `auto_H_refit/` | Automatic proton-rate refit; `auto_H_refit/fits/two_RF_result.json` is the canonical golden | `tests/test_sb_golden.py`, two-state bit-identity checks |
| `field_comparison_600_800_20261003_02/` | 600/800 MHz separate-field study with `inputs/` | `scripts/compare_field_models.py --source` |
| `field_comparison_joint_20261004/` | Joint-field study (no checkpoints, no input copies) | `scripts/compare_joint_fields.py` |
| `two_RF_check/`, `random_H_rates/`, `peakwise_H_fit/`, `proton_relaxation_fit/` | Outputs of the matching `scripts/check_*.py` | docs/reports |

## Rules

- Never edit, regenerate in place, or delete anything here. A rerun writes to a new
  dated directory (`results/<study>_repeat_NN` or a fresh `--out`).
- `executed_source/` folders are hashed copies of the code that produced a result;
  never reformat or fix them, even when the live module changed.
- Two-state single-field results must stay bit-identical to `auto_H_refit`; a
  numerical change that moves them is a regression, not a golden update.
- Ignored by `.gitignore`: `auto_H_refit/source_before/`, `auto_H_refit/legacy_demo/`,
  `*_sample.txt`, joint-study `inputs/` and checkpoints. Keep it that way.
- Values quoted in docs must match the JSON here; verify against the file, not memory.
