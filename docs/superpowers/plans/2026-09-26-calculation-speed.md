# Calculation Speed Implementation Plan

> Execution: native, continuing the user's approved performance proposal.

**Goal:** Reduce redundant model evaluations and Python overhead while retaining the physical models, fit tolerances, dense solver, covariance and default reports.
**Architecture:** Group finite differences across independent residues; batch the two-state Matrix kernel; cache fit-invariant inputs and B1 weights. Keep process parallelism at one level and make PDFs optional.
**Tech Stack:** Python standard library, existing NumPy/SciPy/Matplotlib dependencies.
**Spec:** `session_artifacts/speed_review_20260926/review.md` and the accepted proposal in the conversation.

## Constraints and review focus

- Preserve tracked and unrelated untracked files; no destructive deletion.
- No private SciPy APIs in production, no solver/tolerance/model changes.
- Preserve inactive-residue ordering, unequal spectrum lengths and zero lower bounds.
- Keep scalar/array Matrix behavior and B1 zero/nonzero results consistent.
- Cache only for a fit lifetime; repeated fits after data/active-flag changes must refresh.
- MC workers must fit successfully without spawning nested processes.
- Default CLI output remains; `--no-pdf` skips PDFs but retains numeric results.

## Tasks

1. `test_performance.py`: add one runnable assertion-based regression check for numerical equivalence, grouped evaluation count, output equivalence and real CLI/MC paths. Run before implementation to expose missing optimizations.
2. `fit.py`: add a dense callable grouped forward-difference Jacobian, using the existing nonnegative bounds and per-residue parameter layout. Retain the covariance path, reject irregular parameter lengths and use standard differences for a single active residue.
3. `estmodel.py`: cache immutable B1 arrays; prepare residual inputs once per fit; batch two-state Matrix propagation in bounded chunks; vectorize log calculation; suppress pools inside daemon MC workers and count all spectra for pool eligibility.
4. `run.py`, `mcrun.py`, `README.md`: add/document `--no-pdf` with standard CLI parsing; retain default PDF behavior.
5. Run the regression check, existing three-state verification, actual Matrix/MC fits, fresh before/after benchmarks and diff review. Save measurements under `session_artifacts/speed_implementation_20260926/`.

## Completion

All five tasks completed. Numerical, CLI, MC and three-state checks passed; measurements and limitations are recorded in `session_artifacts/speed_implementation_20260926/results.md`.

Subsequent debugging and refactoring added strict fit/input validation, shared parameter layouts and data loading, and common web fit/MC result handling. GitHub Actions runs the regression scripts on pull requests and pushes to `main`.
