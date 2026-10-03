# Repository Guidelines

## Project Structure & Module Organization

SBONEST fits nitrogen CEST profiles with proton-decoupling sidebands. Python modules live at the repository root: `run.py` dispatches fits, `sbfit.py` implements Sideband fitting, `sideband.py` propagates spin dynamics, and `est_data.py`, `estmodel.py`, and `fit.py` provide shared data/model infrastructure. `server_run.py`, `prepare.py`, and `mcrun.py` support inherited ONEST workflows.

Root-level `test_*.py` and `verify_3state.py` contain regression checks. `example/` holds synthetic inputs and portable configurations; `results/` preserves numerical evidence and figures. `manuscript/sideband_30ppm/` contains manuscript materials. Consult `SBONEST_MANUAL.md`, its Korean counterpart, and `SIDEBAND.md` for configuration details.

## Build, Test, and Development Commands

Use Python 3.12, matching CI. No separate build step is required.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-sideband.txt ruff
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg
python run.py example/sideband_auto_H/two_RF.json --no-pdf
```

These commands install dependencies, limit numerical-library threads, and run the supported Sideband CLI without PDFs. Before rerunning, give the configuration a fresh `Project Name` output prefix. Dataset paths resolve relative to the configuration; output prefixes resolve relative to the working directory.

```bash
ruff check *.py --select F
python -m py_compile *.py
python test_performance.py
python test_debugging.py
python verify_3state.py
python test_sideband.py
python demo_sideband.py --out session_artifacts/sideband_demo_01
```

These reproduce CI lint, syntax, regression, and demonstration checks. Use a new demo output directory each time.

## Coding Style & Naming Conventions

Use four-space indentation, snake_case for new functions/variables, and PascalCase for new classes. Preserve existing public names and scientific configuration keys. Document units and parameter ordering. CI enforces Ruff's `F` rules; no automatic formatter is configured.

## Testing Guidelines

Tests are executable Python scripts using assertions, NumPy comparisons, and `unittest.mock`. Extend the relevant `test_*.py` script with `check_*` functions and invoke them from its runner. Preserve numerical tolerances and independent reference calculations. No coverage percentage is configured; validate physical behavior and numerical parity for model changes.

## Commit & Pull Request Guidelines

Follow history's prefixes: `feat:`, `perf:`, `test:`, and `docs:`. Keep commits focused. PRs should describe changed behavior, link relevant issues, report validation commands/results, and include comparison figures when scientific outputs change. Require passing CI before merging.

## Research Artifact Safety

Archive artifacts instead of deleting them; gitignore large raw data. Preserve archived provenance. Apply manuscript tracked changes as `Donghan Lee` and ground scientific claims in actual code, data, and figures.
