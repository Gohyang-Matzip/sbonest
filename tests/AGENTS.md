# tests/

Executable checks; no pytest. Each file runs as `python tests/NAME.py` from the
repository root and exits nonzero on the first failed assertion. Earned its own
file: 26 modules, a local import convention, and CI shard ownership.

## Where to look

| Need | File |
|------|------|
| Repository root on `sys.path`, single-threaded BLAS, `MPLBACKEND=Agg` | `_env.py` |
| Synthetic Sideband datasets reused by other tests | `test_sb_compare.py` (`write_synthetic`), `test_sideband.py` (`reference`) |
| Synthetic Bruker pseudo-2D data | `test_sb_import.py` (`write_bruker`) |
| Historical golden fit and worker identity | `test_sb_golden.py` (reads `results/auto_H_refit/fits/two_RF_result.json`) |
| Installed wheel and sdist acceptance | `test_sb_distribution.py` (needs `build`, `twine`) |
| Web runner (Flask test client, temporary `SBONEST_JOBS_DIR`) | `test_sb_server.py` |

## Conventions

- Start every test with `import _env` or `from _env import ROOT` before importing
  numerical modules; never set thread variables per file.
- Group assertions in `check_*` functions and call each from the file's runner
  (`if __name__ == '__main__':`); a new check that is not called never runs.
- Parallel tests need the `__main__` guard (spawn workers re-import the file).
- `test_sb_golden.py` and `test_sb_distribution.py` take `--out NEW_DIR` and
  refuse to reuse an existing directory.
- Every new `test_*.py` goes into exactly one CI place: a `checks` shard (core,
  analysis, tools) or the `workflow` job in `.github/workflows/ci.yml`, and into
  the root `AGENTS.md` command list.

## Anti-patterns

- Editing a tolerance or a golden value to make a failure pass; golden tolerances
  (stderr 1e-4, proton rates 1e-3) were set from hosted Linux runners on purpose.
- Fixed sleeps; wait on the process (`process.wait(timeout=...)`) or the file.
- Mocks that skip the numerical path a check claims to cover; `unittest.mock` is
  for refits that must not happen (report paths) and failures to inject.
- Writing outside a `tempfile.mkdtemp()` directory or the `--out` directory.
