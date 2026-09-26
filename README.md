# ONEST: Optimized Novel Exchange Saturation Transfer

**A Web-Based Platform for the Rapid and Robust Analysis of Protein Excited States through CEST Spectroscopy**

## Overview

ONEST (Optimized Novel Exchange Saturation Transfer) is a user-friendly tool designed to automate and accelerate the analysis of Chemical Exchange Saturation Transfer (CEST) NMR spectroscopy data. It facilitates the characterization of "invisible" protein excited states—transient, low-population conformations that are critical for biological function but undetectable by conventional structural methods.

Unlike traditional methods that rely on computationally intensive numerical integration, ONEST utilizes a simultaneous multi-field fitting algorithm leveraging exact analytical solutions for two-state exchange (the Baldwin model). This approach significantly enhances both the speed and robustness of the analysis.

## Key Features

*   **Rapid Analysis:** Utilizes the Baldwin model for exact analytical solutions, avoiding slow numerical integration.
*   **Robust Fitting:** Simultaneous multi-field fitting algorithm for reliable characterization of exchange parameters.
*   **User-Friendly:** Designed to be accessible for researchers analyzing CEST data.
*   **Web-Based & Local Execution:** Can be run as a local web server or via command-line scripts.

## Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/jhyeokchoi/ONEST/
    cd ONEST
    ```

2.  **Install dependencies:**
    Run the following command to install the required packages:
    ```bash
    pip install -r requirements.txt
    ```
    [SciPy 1.9 or newer](https://docs.scipy.org/doc/scipy/release/1.9.0-notes.html#scipy-linalg-improvements) is required for batched matrix exponentials.
    > **Tip:** It is recommended to use a virtual environment to avoid conflicts with other projects.
    > ```bash
    > python3 -m venv venv
    > source venv/bin/activate  # On Windows: venv\Scripts\activate
    > ```

## Usage

### 1. Running the Web Server
To start the web-based interface:
```bash
python server_run.py
```
Access the interface at `http://127.0.0.1:5001`.

### 2. Command Line Execution
You can run the analysis scripts directly using a configuration file:
```bash
python run.py config.json
```

Relative dataset paths are resolved from the configuration file's directory in
all CLI tools. `prepare.py` writes absolute dataset paths so its output can be
saved in another directory. Invalid inputs and fits that fail to converge stop
with an error; failed Monte Carlo fits are excluded from statistics.

For repeated calculations, skip PDF generation while keeping text results:
```bash
python run.py config.json --no-pdf
python mcrun.py config.json 100 4 --no-pdf
```
PDF reports remain enabled by default. Monte Carlo workers fit independently;
the Matrix model does not start nested process pools inside those workers.

> **Tip (For Matrix Method & Benchmarking):**
> When using the **Matrix method** or running benchmarks, it is recommended to force single-core execution to avoid overhead or ensure consistent timing.
> ```bash
> export OMP_NUM_THREADS=1
> python run.py config.json
> ```

### 3. Benchmarking
To benchmark the performance of the model fitting:
```bash
python benchmark.py config.json [profile]
```
- `config.json`: Path to your configuration file.
- `profile`: (Optional) Add this argument to enable cProfile and generate a `benchmark_profile.prof` file.

Run the numerical and performance regression checks (no extra test dependency):
```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python test_performance.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python test_debugging.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python verify_3state.py
```

## Input Generation

For instructions on how to generate input files for ONEST, please refer to the [input4onest](https://github.com/jhyeokchoi/input4onest) repository.

## Directory Structure

*   **`run.py`**: Standard fitting CLI and shared configuration/data loading for fitting, MC and benchmarks.
*   **`server_run.py`**: Flask-based web server for the graphical user interface.
*   **`mcrun.py`**: Script for Monte Carlo simulations (supports multiprocessing).
*   **`benchmark.py`**: Tool for benchmarking model performance.
*   **`estmodel.py`**: Calculation engines (Baldwin, Matrix, NoEx and three-state models), residuals and reports.
*   **`fit.py`**: Shared parameter layouts, initial guesses, bounds, Jacobian and least-squares solver.
*   **`est_data.py`**: Spectrum file parsing and dataset storage.
*   **`example/`**: Contains example datasets (`syn10.txt`, `syn100.txt`) for testing.
*   **`requirements.txt`**: Python dependencies.

## Citation

If you use ONEST in your research, please cite:

> Choi, J., Lee, SY., Han, K., Carneiro, M. G., Ryu, KS., & Lee, D. "ONEST: A Web-Based Platform for the Rapid and Robust Analysis of Protein Excited States through CEST Spectroscopy." (in preparation)

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.
