#!/usr/bin/env python
# fit.py — fitting logic for estmodel.py
# Original estmodel.py by Donghan Lee 2015; scipy.optimize.least_squares port 2025

import numpy as np
from scipy.optimize import least_squares


def _res_guess(fs):
    """Initial (dG, R1) guess from a residue's first spectrum: dG at the intensity
    minimum, R1 from the decay of the maximum intensity."""
    dG = 0.0
    if len(fs.offset) > 0 and len(fs.int) > 0:
        dG = fs.offset[np.argmin(fs.int)]
    maxint = np.max(fs.int) if len(fs.int) > 0 else 0.0
    r1 = -np.log(maxint) / fs.T if (maxint > 0 and fs.T > 0) else 1.0
    return dG, r1


def generate_initial_parameters(model_instance, initConf):
    """
    Generates the initial parameter vector p0 based on initConf.
    For 2-state exchange models a (kex, pB) grid search picks the best start;
    NoEx and 3-state models use per-residue heuristics / defaults.
    """
    model_instance.selMethod(initConf)
    active = [r for r in model_instance.dataset.res if r.active]

    p_initial_list = []
    if model_instance.method == "NoEx":
        for res_obj in active:
            if not res_obj.estSpecs or not res_obj.estSpecs[0].int:
                if model_instance.verbose:
                    print(
                        f"Warning (fit.py): No data for {res_obj.label}, using defaults for NoEx."
                    )
                p_initial_list.extend([0.0, 1.0, 10.0])
                continue
            fs = res_obj.estSpecs[0]
            dG, r1 = _res_guess(fs)
            p_initial_list.extend([dG, r1, fs.initr2a])

    elif model_instance.method in ["Matrix_3st_Linear", "Matrix_3st_Triangle"]:
        # Grid search is too expensive for 4-6 exchange rates; use defaults.
        n_rates = 4 if model_instance.method == "Matrix_3st_Linear" else 6
        p_initial_list.extend([10.0] * n_rates)
        for r_o in active:
            if not r_o.estSpecs or not r_o.estSpecs[0].int:
                p_initial_list.extend(
                    [0.0, 0.1, 0.2, 1.0, 10.0, 20.0, 20.0]
                )  # dG, dw, dwC, R1, R2a, R2b, R2c
                continue
            fs = r_o.estSpecs[0]
            dG, r1 = _res_guess(fs)
            # dwC guess: initdw*2
            p_initial_list.extend(
                [dG, fs.initdw, fs.initdw * 2.0, r1, fs.initr2a, fs.initr2b, fs.initr2b]
            )

    else:  # 2-state exchange (Baldwin, Matrix): grid search over kex and pB
        kex_c = initConf.get("kex", {"min": 10.0, "max": 400.0, "nsteps": 6})
        pb_c = initConf.get("pB", {"min": 0.01, "max": 0.1, "nsteps": 6})
        kex_v = np.linspace(kex_c["min"], kex_c["max"], int(kex_c.get("nsteps", 6)))
        pB_v = np.linspace(pb_c["min"], pb_c["max"], int(pb_c.get("nsteps", 6)))

        # Per-residue params don't depend on the grid point: compute once.
        have_all_data = all(r.estSpecs and r.estSpecs[0].int for r in active)
        res_params = []
        for r_o in active:
            if not r_o.estSpecs or not r_o.estSpecs[0].int:
                res_params.extend(
                    [0.0, 0.1, 1.0, 10.0, 20.0]
                )  # dG, dw, R1, R2a, R2b defaults
                continue
            fs = r_o.estSpecs[0]
            dG, r1 = _res_guess(fs)
            res_params.extend([dG, fs.initdw, r1, fs.initr2a, fs.initr2b])

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


def perform_least_squares_fit(model_instance, p0_initial):
    """
    Least squares fit via scipy.optimize.least_squares using model_instance.errFunc.
    Returns (optimized params, covariance matrix).
    """
    p0 = np.asarray(p0_initial, dtype=float)
    method = model_instance.method
    num_active = sum(1 for r in model_instance.dataset.res if r.active)

    # Bounds: rate constants >= 0; per residue, leading shift params free, rates >= 0.
    n_glob = {"NoEx": 0, "Matrix_3st_Linear": 4, "Matrix_3st_Triangle": 6}.get(
        method, 2
    )
    n_free = {"NoEx": 1, "Matrix_3st_Linear": 3, "Matrix_3st_Triangle": 3}.get(
        method, 2
    )
    n_res = {"NoEx": 3, "Matrix_3st_Linear": 7, "Matrix_3st_Triangle": 7}.get(method, 5)

    bounds_min = [0.0] * n_glob + (
        [-np.inf] * n_free + [0.0] * (n_res - n_free)
    ) * num_active
    expected_len = len(bounds_min)
    if len(p0) != expected_len:
        if model_instance.verbose:
            print(
                f"Warning (fit.py): p0 length ({len(p0)}) != expected ({expected_len}) for method {method}. Adjusting bounds."
            )
        if len(p0) < expected_len:
            bounds_min = bounds_min[: len(p0)]
        else:
            bounds_min.extend([-np.inf] * (len(p0) - expected_len))
    scipy_bounds = (np.array(bounds_min), np.full(len(p0), np.inf))

    if model_instance.verbose:
        print(
            f"Initiating fit (in fit.py, method: {method}): Using scipy.optimize.least_squares with {len(p0)} params."
        )

    result = None
    p1_optimized = np.copy(p0)
    covariance_matrix = np.full((len(p0), len(p0)), np.nan)

    try:
        result = least_squares(
            model_instance.errFunc,
            x0=p0,
            bounds=scipy_bounds,
            method="trf",
            ftol=1e-9,
            xtol=1e-9,
            gtol=1e-9,
            verbose=2 if model_instance.verbose else 0,
        )
        p1_optimized = result.x
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
            details = f"Jacobian shape: {result.jac.shape if result.jac is not None else 'None'}"
            print(
                f"Jacobian not available or has unexpected dimensions for covariance calculation (fit.py). {details}"
            )
    except Exception as e:
        if model_instance.verbose:
            print(f"Error during least_squares fitting (in fit.py): {e}")

    if model_instance.verbose:
        print("Fit (in fit.py) completed.")
        if result:
            print(f"SciPy Status: {result.status} ({result.message})")
            print(
                f"Final cost: {result.cost:.4e}, NFEV: {result.nfev}, NJEV: {getattr(result, 'njev', 'N/A')}"
            )
        else:
            print("Fitting process did not yield a result object.")

    return p1_optimized, covariance_matrix
