#!/usr/bin/env python
# fit.py — fitting logic for estmodel.py
# Original estmodel.py by Donghan Lee 2015; scipy.optimize.least_squares port 2025

import numpy as np
from scipy.optimize import least_squares


VALID_METHODS = [
    "Baldwin",
    "Matrix",
    "NoEx",
    "Matrix_3st_Linear",
    "Matrix_3st_Triangle",
]

# Per-method flat-parameter layout: (global keys, per-residue keys, per-residue defaults)
PARAM_LAYOUT = {
    "NoEx": ([], ["dGs", "r1s", "r2as"], [0.0, 1.0, 10.0]),
    "Baldwin": (
        ["kab", "kba"],
        ["dGs", "dws", "r1s", "r2as", "r2bs"],
        [0.0, 0.1, 1.0, 10.0, 20.0],
    ),
    "Matrix_3st_Linear": (
        ["kab", "kba", "kbc", "kcb"],
        ["dGs", "dws", "dwCs", "r1s", "r2as", "r2bs", "r2cs"],
        [0.0, 0.1, 0.2, 1.0, 10.0, 20.0, 20.0],
    ),
    "Matrix_3st_Triangle": (
        ["kab", "kba", "kbc", "kcb", "kca", "kac"],
        ["dGs", "dws", "dwCs", "r1s", "r2as", "r2bs", "r2cs"],
        [0.0, 0.1, 0.2, 1.0, 10.0, 20.0, 20.0],
    ),
}
PARAM_LAYOUT["Matrix"] = PARAM_LAYOUT["Baldwin"]


def _res_guess(fs):
    """Initial (dG, R1) guess from a residue's first spectrum: dG at the intensity
    minimum, R1 from the decay of the maximum intensity."""
    dG = 0.0
    if len(fs.offset) > 0 and len(fs.int) > 0:
        dG = fs.offset[np.argmin(fs.int)]
    maxint = np.max(fs.int) if len(fs.int) > 0 else 0.0
    r1 = max(0.0, -np.log(maxint) / fs.T) if (maxint > 0 and fs.T > 0) else 1.0
    return dG, r1


def generate_initial_parameters(model_instance, initConf):
    """
    Generates the initial parameter vector p0 based on initConf.
    For 2-state exchange models a (kex, pB) grid search picks the best start;
    NoEx and 3-state models use per-residue heuristics / defaults.
    """
    model_instance.selMethod(initConf)
    active = [r for r in model_instance.dataset.res if r.active]
    gkeys, rkeys, defaults = PARAM_LAYOUT[model_instance.method]

    # Local guesses are shared by all models and independent of the rate grid.
    res_params = []
    have_all_data = True
    for res in active:
        if not res.estSpecs or not res.estSpecs[0].int:
            have_all_data = False
            res_params.extend(defaults)
            if model_instance.method == "NoEx" and model_instance.verbose:
                print(
                    f"Warning (fit.py): No data for {res.label}, using defaults for NoEx."
                )
            continue
        fs = res.estSpecs[0]
        dG, r1 = _res_guess(fs)
        values = {
            "dGs": dG,
            "dws": fs.initdw,
            "dwCs": fs.initdw * 2.0,
            "r1s": r1,
            "r2as": fs.initr2a,
            "r2bs": fs.initr2b,
            "r2cs": fs.initr2b,
        }
        res_params.extend(values[key] for key in rkeys)

    if model_instance.method == "NoEx":
        p_initial_list = res_params
    elif model_instance.method in ("Matrix_3st_Linear", "Matrix_3st_Triangle"):
        # Grid search is too expensive for 4-6 exchange rates; use defaults.
        p_initial_list = [10.0] * len(gkeys) + res_params

    else:  # 2-state exchange (Baldwin, Matrix): grid search over kex and pB
        kex_c = initConf.get("kex", {"min": 10.0, "max": 400.0, "nsteps": 6})
        pb_c = initConf.get("pB", {"min": 0.01, "max": 0.1, "nsteps": 6})
        kex_v = np.linspace(kex_c["min"], kex_c["max"], int(kex_c.get("nsteps", 6)))
        pB_v = np.linspace(pb_c["min"], pb_c["max"], int(pb_c.get("nsteps", 6)))

        minChi2, best_p0_grid = float("inf"), None
        if have_all_data:
            for ikv in kex_v:
                for ipbv in pB_v:
                    kab = ipbv * ikv
                    kba = ikv - kab
                    if kab < 0 or kba < 0:
                        continue
                    curr_p_g = [kab, kba] + res_params
                    try:
                        model_instance.errFunc(np.array(curr_p_g, dtype=float))
                        if model_instance.chi2 < minChi2:
                            minChi2 = model_instance.chi2
                            best_p0_grid = curr_p_g
                    except Exception as e:
                        if model_instance.verbose:
                            print(
                                f"Error in errFunc during grid search for initial params (fit.py): {e}"
                            )

        if best_p0_grid:
            p_initial_list = best_p0_grid
        else:
            if model_instance.verbose:
                print("Grid search failed for initial params (fit.py), using defaults.")
            default_kex = (kex_c["min"] + kex_c["max"]) / 2.0
            default_pB = (pb_c["min"] + pb_c["max"]) / 2.0
            p_initial_list = [
                default_pB * default_kex,
                (1.0 - default_pB) * default_kex,
            ] + res_params

    if not p_initial_list:
        raise ValueError("Failed to generate any initial parameters (fit.py).")
    if model_instance.verbose:
        print("Initial parameter generation (in fit.py) finished.")
    return np.array(p_initial_list, dtype=float)


def _block_jacobian(fun, residue_sizes, n_global, n_local, relative_step=None):
    """Group independent residue parameters while returning a dense Jacobian."""
    rows = np.arange(sum(residue_sizes))
    local_columns = np.repeat(
        n_global + n_local * np.arange(len(residue_sizes)), residue_sizes
    )
    columns = [np.full(len(rows), i) for i in range(n_global)]
    columns.extend(local_columns + i for i in range(n_local))
    groups = [(c, np.unique(c)) for c in columns]

    def jacobian(p):
        base = fun(p)
        jac = np.zeros((len(base), len(p)))
        # Negative parameters are unbounded shifts; bounded rates always step up.
        step = (
            (np.sqrt(np.finfo(float).eps) if relative_step is None else relative_step)
            * np.where(p >= 0, 1.0, -1.0)
            * np.maximum(1.0, np.abs(p))
        )
        for column, changed in groups:
            shifted = p.copy()
            shifted[changed] += step[changed]
            delta = shifted - p
            jac[rows, column] = (fun(shifted) - base) / delta[column]
        return jac

    return jacobian


def perform_least_squares_fit(model_instance, p0_initial):
    """
    Least squares fit via scipy.optimize.least_squares using model_instance.errFunc.
    Returns (optimized params, covariance matrix).
    """
    p0 = np.asarray(p0_initial, dtype=float)
    method = model_instance.method
    num_active = sum(1 for r in model_instance.dataset.res if r.active)

    # Bounds follow the same layout as initialization and unpacking.
    gkeys, rkeys, _ = PARAM_LAYOUT[method]
    n_glob, n_res = len(gkeys), len(rkeys)
    local_bounds = [-np.inf if k in ("dGs", "dws", "dwCs") else 0.0 for k in rkeys]
    bounds_min = [0.0] * n_glob + local_bounds * num_active
    expected_len = len(bounds_min)
    if p0.ndim != 1 or p0.size != expected_len:
        raise ValueError(
            f"Expected {expected_len} parameters for {method}; got shape {p0.shape}."
        )
    scipy_bounds = (np.array(bounds_min), np.full(len(p0), np.inf))
    residue_sizes = [
        sum(len(es.offset) for es in r.estSpecs)
        for r in model_instance.dataset.res
        if r.active
    ]
    if not residue_sizes or any(size == 0 for size in residue_sizes):
        raise ValueError(
            "Select at least one residue; every selected residue must contain data."
        )
    jacobian = (
        _block_jacobian(model_instance.errFunc, residue_sizes, n_glob, n_res)
        if num_active > 1
        else "2-point"
    )

    if model_instance.verbose:
        print(
            f"Initiating fit (in fit.py, method: {method}): Using scipy.optimize.least_squares with {len(p0)} params."
        )

    covariance_matrix = np.full((len(p0), len(p0)), np.nan)

    result = least_squares(
        model_instance.errFunc,
        x0=p0,
        bounds=scipy_bounds,
        jac=jacobian,
        method="trf",
        ftol=1e-9,
        xtol=1e-9,
        gtol=1e-9,
        verbose=2 if model_instance.verbose else 0,
    )
    if not result.success:
        raise RuntimeError(
            f"Fit did not converge (status {result.status}): {result.message}"
        )
    p1_optimized = result.x
    model_instance.errFunc(p1_optimized)
    if result.jac is not None and result.jac.shape[1] == len(p0):
        try:
            jtj = result.jac.T @ result.jac
            if np.linalg.cond(jtj) < 1 / np.finfo(jtj.dtype).eps:
                covariance_matrix = np.linalg.inv(jtj)
            elif model_instance.verbose:
                print(
                    "Jacobian^T * Jacobian is singular or ill-conditioned (fit.py); covariance calculation failed."
                )
        except (np.linalg.LinAlgError, ValueError) as e:
            if model_instance.verbose:
                print(f"Covariance matrix computation failed (fit.py): {e}")
    elif model_instance.verbose:
        details = (
            f"Jacobian shape: {result.jac.shape if result.jac is not None else 'None'}"
        )
        print(
            f"Jacobian not available or has unexpected dimensions for covariance calculation (fit.py). {details}"
        )

    if model_instance.verbose:
        print("Fit (in fit.py) completed.")
        print(f"SciPy Status: {result.status} ({result.message})")
        print(
            f"Final cost: {result.cost:.4e}, NFEV: {result.nfev}, NJEV: {getattr(result, 'njev', 'N/A')}"
        )

    return p1_optimized, covariance_matrix
