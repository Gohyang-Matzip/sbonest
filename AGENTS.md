# Repository Guidelines

## Project Structure & Module Organization

Root modules include `run.py` (CLI dispatch), `sbfit.py` (Sideband fitting), `sideband.py` (spin propagation), and `est_data.py`, `estmodel.py`, and `fit.py` (shared infrastructure). `server_run.py`, `prepare.py`, and `mcrun.py` support inherited ONEST workflows.

Root `test_*.py` and `verify_3state.py` contain regression checks. `example/` holds portable synthetic examples; `results/` preserves evidence. `compare_field_models.py` compares 600/800 MHz models. `manuscript/sideband_30ppm/` contains published manuscript materials. Update both `SBONEST_MANUAL` languages when changing documented behavior.

## Build, Test, and Development Commands

Use Python 3.12, matching CI. No separate build step is required.

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-sideband.txt -c constraints-sideband.txt ruff
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg
python run.py example/sideband_auto_H/two_RF.json --no-pdf
```

Use a fresh `Project Name` before rerunning. Dataset paths resolve from the configuration; output prefixes resolve from the working directory.

```bash
ruff check *.py --select F
python -m py_compile *.py
python test_performance.py
python test_debugging.py
python verify_3state.py
python test_sideband.py
python test_grouped_jacobian.py
python test_benchmark.py
python test_sb_analysis.py
python test_sb_output.py
python demo_sideband.py --out session_artifacts/sideband_demo_01
```

These match CI checks; use fresh output directories. Reproduce the field comparison with:

```bash
python compare_field_models.py \
  --source results/field_comparison_600_800_20261003_02/inputs \
  --out results/field_comparison_repeat_01
```

Explicit `--source` is necessary in fresh checkouts. This longer study is separate from CI.

## Coding Style & Naming Conventions

Use four-space indentation, snake_case functions/variables, and PascalCase new classes. Preserve public names, configuration keys, units, and parameter ordering. CI enforces Ruff `F`; no formatter is configured. Avoid bulk reformatting executed research sources referenced by hashes.

## Testing Guidelines

Tests use executable Python, assertions, NumPy comparisons, and `unittest.mock`. Add `check_*` functions to the relevant script and invoke them from its runner. Preserve numerical tolerances and independent references. No coverage percentage is configured. For documentation edits, verify links and table values against source JSON.

## Commit & Pull Request Guidelines

Use history's `feat:`, `perf:`, `test:`, and `docs:` prefixes. PRs should explain changes, link relevant issues, report validation, and show changed figures. Require passing CI; after merging, verify the merged commit's CI and synchronize `main`.

## Research Artifact Safety

Archive artifacts instead of deleting them; gitignore large raw data. Preserve calculation-time hashes. The 600/800 MHz benchmark uses fixed H rates and separate field fits; current default examples fit H rates automatically. Label synthetic evidence and within-model uncertainty. Apply manuscript tracked changes as `Donghan Lee`.
