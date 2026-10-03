For SBONEST sideband fitting, use the [SBONEST manual](SBONEST_MANUAL.md).
This inherited manual covers the original ONEST models.
For the matched 600/800 MHz SBONEST–ONEST comparison, see section 7.1 of the
[SBONEST manual](SBONEST_MANUAL.md) and the
[comparison report](results/field_comparison_600_800_20261003_02/REPORT.txt).
That benchmark fits RF scale with a shared optimizer; its settings differ from
the standard ONEST workflow below.

# ONEST Step-by-Step Manual

ONEST analyzes CEST NMR data to characterize invisible protein excited states via simultaneous multi-field fitting (Baldwin analytical model, plus Matrix and 3-state variants).

## 1. Installation

```bash
git clone https://github.com/jhyeokchoi/ONEST/
cd ONEST
python3 -m venv venv
source venv/bin/activate   # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

Dependencies: `numpy`, `scipy`, `matplotlib` (CLI); `flask`, `werkzeug` (web server only).

## 2. Prepare the input data file

One text file per B0 field / experiment. Format (see `example/syn10.txt`):

```
80.12           # B0 field (MHz)
0.400000        # saturation time T (s)
10.00  0.20     # B1 (v1) in Hz and its error
#offset(ppm)    Intensity    error
# A1 R2a: 10.0 R2b: 10.0 dw: -5.0
 100.000    0.672    0.006
 100.201    0.682    0.006
 ...
```

- **Line 1**: B0 field in MHz. **Line 2**: saturation time in seconds. **Line 3**: B1 field (Hz) and its error.
- **Line 4**: column header comment (skipped).
- **Residue header**: starts each residue block. Two forms:
  - Simple: `# A1`
  - Full (recommended): `# A1 R2a: 10.0 R2b: 10.0 dw: -5.0` — provides initial guesses. Supplying an initial `dw` (with the correct sign) avoids sign-flipped local minima in the fit.
- **Data lines**: `offset(ppm)  intensity  error`, one per point.
- Multiple fields: put each field in its own file and pass all files to `prepare.py`; residue names must match across files.

To generate input files from raw spectra, see [input4onest](https://github.com/jhyeokchoi/input4onest).

## 3. Generate the configuration (prepare.py)

```bash
python prepare.py example/syn10.txt > config.json
# multiple fields:
python prepare.py field1.txt field2.txt > config.json
# choose a calculation method (default: Baldwin):
python prepare.py --method Matrix example/syn10.txt > config.json
```

Methods: `Baldwin` (2-state analytical, fast — default), `Matrix` (2-state numerical), `NoEx` (no exchange), `Matrix_3st_Linear`, `Matrix_3st_Triangle` (3-state).

## 4. Edit config.json (optional)

```json
{
    "Project Name": "default",
    "init": {
        "kex": {"min": 10.0, "max": 400.0, "nsteps": 6},
        "pB":  {"min": 0.01, "max": 0.1,  "nsteps": 6},
        "Method": "Baldwin"
    },
    "datasets": ["example/syn10.txt"],
    "residues": [{"name": "A1", "flag": "on"}, ...]
}
```

- `Project Name` — prefix for all output files.
- `init.kex`, `init.pB` — grid-search ranges for the global exchange rate and excited-state population. If you expect kex > 400 s⁻¹, raise `max` (e.g. 2000).
- `residues[].flag` — set `"off"` to exclude a residue from the fit.

## 5. Run the fit

```bash
python run.py config.json
```

For the Matrix method, force single-threaded BLAS to avoid multiprocessing overhead:

```bash
export OMP_NUM_THREADS=1
python run.py config.json
```

Outputs (with `Project Name` = `default`):

| File | Content |
|------|---------|
| `default_data.pdf` | Raw data plots per residue |
| `default_result.txt` | Fitted parameters, Chi2, dof, per-residue results |
| `default.pdf` | Data with fitted curves |

Check `default_result.txt`: global `kex`, `pB`, per-residue `R2a`, `R2b`, `dw`, and Chi2. Verify each residue's fitted `dw` sign is physically reasonable — a flipped sign with a high per-residue chi2 signals a local minimum (fix by giving the initial `dw` in the full residue header, step 2).

## 6. Error estimation via Monte Carlo (mcrun.py)

```bash
python mcrun.py config.json 100          # 100 MC runs, all CPU cores
python mcrun.py config.json 100 4        # limit to 4 processes
```

Each run re-samples the data within its error bars and refits. Additional outputs: `default_mc.txt` (parameter means ± std dev) and `default_mcmean.pdf`.

Note (macOS): run from a saved script/file, not piped stdin — the `spawn` start method cannot re-import `<stdin>`.

## 7. Web interface (optional)

```bash
python server_run.py
```

Open `http://127.0.0.1:5001`. Upload data, configure, and fit through the browser.

## 8. Benchmarking (optional)

```bash
export OMP_NUM_THREADS=1
python benchmark.py config.json           # timing
python benchmark.py config.json profile   # + cProfile → benchmark_profile.prof
```

## Troubleshooting

- **`Field line parse error` / `V1 line parse error`** — check the 3 header lines: one number, one number, two numbers, in that order.
- **`Residue X not found in dataset`** — residue name in `config.json` doesn't match the data file headers.
- **One residue dominates Chi2** — likely a `dw` sign-flip local minimum; supply the initial `dw` via the full residue header.
- **kex hits the grid boundary** — widen `init.kex.max` and rerun.
