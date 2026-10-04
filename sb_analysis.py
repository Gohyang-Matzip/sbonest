"""Optional multistart fits and constrained, nuisance-refitted likelihood profiles.

Residuals use supplied absolute intensity errors. Profile delta-chi2 values are
within-model likelihood differences, not confidence intervals or proof that the
base fit is a global minimum. All unsuccessful attempts are retained.
"""

import copy
import numbers

import numpy as np
from scipy.optimize import least_squares


class MultiStartError(RuntimeError):
    """All fits failed; ``attempts`` retains their starting vectors and errors."""

    def __init__(self, attempts):
        super().__init__("No multistart attempt converged")
        self.attempts = attempts


_STATE = ("result", "free", "lower", "upper", "initial_parameters", "rank",
          "condition", "npar", "nvar", "dof", "chi2")


def _state(model):
    return {name: copy.copy(getattr(model, name)) if name != "result"
            else model.result for name in _STATE if hasattr(model, name)}


def _restore(model, state):
    for name, value in state.items():
        setattr(model, name, value)


def _chi2(model, p):
    residual = np.asarray(model.errFunc(p), dtype=float)
    if not np.isfinite(residual).all():
        raise ValueError("Non-finite analysis residual")
    value = float(np.sum(residual**2))
    if not np.isfinite(value):
        raise ValueError("Non-finite analysis chi2")
    return value


def _number(value):
    return isinstance(value, numbers.Real) and not isinstance(value, (bool, np.bool_)) and np.isfinite(value)


def _multistart_settings(model, settings):
    if not isinstance(settings, dict) or set(settings) - {"starts", "seed", "random_starts"}:
        raise ValueError("multistart allows starts, seed, and random_starts")
    count = settings.get("random_starts", 0)
    if isinstance(count, bool) or not isinstance(count, int) or count < 0:
        raise ValueError("multistart.random_starts must be a nonnegative integer")
    seed = settings.get("seed")
    if "seed" in settings and (isinstance(seed, bool) or not isinstance(seed, int) or seed < 0):
        raise ValueError("multistart.seed must be a nonnegative integer")
    if count and seed is None:
        raise ValueError("Random multistart requires an explicit seed")
    starts = settings.get("starts", [])
    if not isinstance(starts, list):
        raise ValueError("multistart.starts must be a list of parameter mappings")
    vary = model.config["init"].get("vary", model.parameter_names)
    for entry in starts:
        if not isinstance(entry, dict) or any(name not in model.parameter_names for name in entry):
            raise ValueError("Each multistart start must map known parameter names to values")
        if any(name not in vary for name in entry):
            raise ValueError("Multistart starts cannot change fixed parameters")
        if not all(_number(value) for value in entry.values()):
            raise ValueError("Multistart starting values must be finite numbers")
    return starts, count, seed


def _random_start(model, base, rng):
    """Uniform finite bounds; otherwise local normals, reflected at a lone bound.

    The standard deviation is 0.1 ppm for peak_ppm, max(0.25*abs(value),
    0.1) ppm for dw_ppm, and 0.25*max(abs(value), 1) for other parameters.
    Fixed parameters are never sampled.
    """
    start = base.copy()
    for i in model.free:
        lo, hi, value = model.lower[i], model.upper[i], base[i]
        if np.isfinite(lo) and np.isfinite(hi):
            start[i] = rng.uniform(lo, hi)
            continue
        name = model.parameter_names[i]
        scale = (0.1 if name.endswith(".peak_ppm") else
                 max(0.25*abs(value), 0.1) if name.endswith(".dw_ppm") else
                 0.25*max(abs(value), 1.0))
        sample = rng.normal(value, scale)
        start[i] = (lo + abs(sample - lo) if np.isfinite(lo) else
                    hi - abs(sample - hi) if np.isfinite(hi) else sample)
    return start


def fit_multistart(model, settings):
    """Fit the configured baseline plus starts; restore and return the best fit.

    ``model.fit`` must expose full ``initial_parameters``, ``lower``, and
    ``upper`` arrays before optimization. Arrays in result rows follow
    ``model.parameter_names``. All failed runs raise MultiStartError with the
    complete attempt list, allowing callers to save failure evidence.
    """
    starts, count, seed = _multistart_settings(model, settings)
    cfg = {k: v for k, v in model.config["init"].items() if k not in ("multistart", "profile")}
    attempts, best = [], None

    def attempt(start, source):
        nonlocal best
        row = {"index": len(attempts), "source": source, "start": None,
               "parameters": None, "success": False, "chi2": None,
               "nfev": None, "status": None, "message": "", "selected": False}
        model.result = None
        prior_initial = getattr(model, "initial_parameters", None)
        if start is not None:
            row["start"] = start.tolist()
        try:
            use_cfg = cfg if start is None else {k: v for k, v in cfg.items() if k != "initial"}
            p, covariance = model.fit(p0=start, fitting_config=use_cfg)
            if model.result is not None and not model.result.success:
                raise RuntimeError(str(model.result.message))
            fit_state = _state(model)
            value = _chi2(model, p)
            _restore(model, fit_state)
            row.update(success=True, chi2=value, parameters=np.asarray(p).tolist())
            if best is None or value < best[0]:
                best = (value, np.asarray(p).copy(), np.asarray(covariance).copy(), fit_state, row["index"])
        except Exception as exc:
            row["message"] = f"{type(exc).__name__}: {exc}"
        if (row["start"] is None and hasattr(model, "initial_parameters")
                and (row["success"] or model.initial_parameters is not prior_initial)):
            initial = np.asarray(model.initial_parameters)
            if np.isfinite(initial).all():
                row["start"] = initial.tolist()
        result = model.result
        if result is not None:
            row.update(nfev=int(result.nfev), status=int(result.status))
            if row["success"]:
                row["message"] = str(result.message)
            elif row["start"] is not None:
                last = np.asarray(row["start"]).copy()
                last[model.free] = result.x
                if np.isfinite(last).all():
                    row["parameters"] = last.tolist()
                value = float(np.sum(np.asarray(result.fun)**2))
                if np.isfinite(value):
                    row["chi2"] = value
        attempts.append(row)

    attempt(None, "initial")
    if attempts[0]["start"] is None or not all(hasattr(model, key) for key in ("initial_parameters", "lower", "upper")):
        raise MultiStartError(attempts)
    base = np.asarray(model.initial_parameters).copy()
    for entry in starts:
        start = base.copy()
        for name, value in entry.items():
            start[model.parameter_names.index(name)] = value
        attempt(start, "explicit")
    rng = np.random.default_rng(seed)
    for _ in range(count):
        attempt(_random_start(model, base, rng), "random")
    if best is None:
        raise MultiStartError(attempts)
    _, p, covariance, state, selected = best
    model.errFunc(p)
    _restore(model, state)
    attempts[selected]["selected"] = True
    return p, covariance, attempts


def _profile_settings(settings, names):
    if not isinstance(settings, dict) or set(settings) - {"kex", "pB", "v1n_scale"}:
        raise ValueError("profile allows kex, pB, and v1n_scale grids")
    for name, values in settings.items():
        if name == "v1n_scale" and name not in names:
            raise ValueError("v1n_scale profiling requires scale RF mode")
        if not isinstance(values, list) or not values or not all(_number(v) for v in values):
            raise ValueError(f"profile.{name} must be a nonempty list of finite numbers")


def _constrained_coordinates(model, p, name, target):
    """Return a bounded coordinate map that satisfies the exact constraint."""
    names = model.parameter_names
    lower, upper = np.asarray(model.lower), np.asarray(model.upper)
    free = list(model.free)
    base = p.copy()
    excluded, rate = set(), None
    tolerance = np.zeros(len(p))
    if name == "v1n_scale":
        i = names.index(name)
        if i not in free and target != p[i]:
            raise ValueError("Profile target changes a fixed RF scale")
        base[i] = target
        excluded.add(i)
    else:
        a, b = names.index("kab"), names.index("kba")
        if name == "kex" and target <= 0:
            raise ValueError("kex profile targets must be positive")
        if name == "pB" and not 0 <= target <= 1:
            raise ValueError("pB profile targets must lie between 0 and 1")
        varying = [i for i in (a, b) if i in free]
        excluded.update(varying)
        if len(varying) == 2:
            if name == "kex":
                lo, hi = max(lower[a], target-upper[b]), min(upper[a], target-lower[b])
                rate = (lo, hi, p[a], lambda x: (x, target-x), "kab (kba=kex-kab)")
            else:
                lo, hi = 0., np.inf
                for i, coefficient in ((a, target), (b, 1-target)):
                    if coefficient:
                        lo, hi = max(lo, lower[i]/coefficient), min(hi, upper[i]/coefficient)
                    elif not lower[i] <= 0 <= upper[i]:
                        raise ValueError("Population endpoint violates rate bounds")
                rate = (lo, hi, p[a]+p[b], lambda x: (target*x, (1-target)*x), "kex (fixed pB)")
            lo, hi, initial, transform, label = rate
            # Derived limits can cross by an ULP at an exactly feasible boundary.
            # Preserve the target algebra; allow only floating-point roundoff.
            slack = 16*np.finfo(float).eps*max(np.finfo(float).tiny, abs(lo), abs(hi)) if np.isfinite([lo, hi]).all() else 0.
            if not np.isfinite(lo) or lo > hi + slack:
                raise ValueError("Profile constraint has no feasible rate interval")
            if lo >= hi:
                base[a], base[b] = transform(lo)
                rate = None
            else:
                initial = np.clip(initial, lo, hi)
                base[a], base[b] = transform(initial)
                rate = (lo, hi, initial, transform, label)
        elif len(varying) == 1:
            i = varying[0]
            fixed = b if i == a else a
            if name == "kex":
                base[i] = target-p[fixed]
            else:
                coefficient, other = (1-target, target) if i == a else (target, 1-target)
                if coefficient == 0:
                    if p[fixed] != 0:
                        raise ValueError("Population target conflicts with a fixed rate")
                    excluded.remove(i)  # Endpoint already enforced by the zero fixed rate.
                else:
                    base[i] = other*p[fixed]/coefficient
        else:
            total = p[a]+p[b]
            actual = total if name == "kex" else p[a]/total if total > 0 else np.nan
            if not np.isclose(actual, target, rtol=16*np.finfo(float).eps, atol=0.):
                raise ValueError("Profile target conflicts with fixed exchange rates")
        if base[a] + base[b] <= 0:
            raise ValueError("Exchange population is undefined at zero total rate")
        tolerance[[a, b]] = 16*np.finfo(float).eps*max(np.finfo(float).tiny, abs(base[a]), abs(base[b]))
    if (not np.isfinite(base).all() or np.any(base < lower-tolerance)
            or np.any(base > upper+tolerance)):
        raise ValueError("Profile target violates parameter bounds")

    coordinates, labels, q0, lo, hi = [], [], [], [], []
    rate_added = False
    for i in free:
        if i in excluded:
            if rate is not None and not rate_added:
                coordinates.append(None)
                labels.append(rate[4])
                lo.append(rate[0])
                hi.append(rate[1])
                q0.append(rate[2])
                rate_added = True
        else:
            coordinates.append(i)
            labels.append(names[i])
            lo.append(lower[i])
            hi.append(upper[i])
            q0.append(base[i])

    def expand(q):
        out = base.copy()
        for coordinate, value in zip(coordinates, q):
            if coordinate is None:
                out[a], out[b] = rate[3](value)
            else:
                out[coordinate] = value
        if name != "v1n_scale":
            indices = [a, b]
            rates = out[indices]
            slack = 16*np.finfo(float).eps*max(np.finfo(float).tiny, *np.abs(rates))
            if (np.any(rates < lower[indices]-slack)
                    or np.any(rates > upper[indices]+slack)):
                raise ValueError("Profile constraint violates rate bounds")
            # Repair only ULP-scale reconstruction excursions, including a tiny
            # negative rate; the requested constraint changes by roundoff only.
            for i in excluded:
                out[i] = np.clip(out[i], lower[i], upper[i])
        return out

    return np.asarray(q0), np.asarray(lo), np.asarray(hi), expand, labels


def profile_likelihood(model, p, settings):
    """Nuisance-refit each requested target, preserving the original fitted state.

    Failed targets have null chi2/delta_chi2; negative successful differences
    remain negative and flag a base fit that this scan has improved upon.
    Every grid point starts from the supplied base fit after making its exact
    constraint feasible. No warm start or automatic confidence interval is used.
    """
    _profile_settings(settings, model.parameter_names)
    p = np.asarray(p, dtype=float).copy()
    if p.shape != (len(model.parameter_names),) or not np.isfinite(p).all():
        raise ValueError("Invalid base profile parameter vector")
    state = _state(model)
    old_data = getattr(model, "_fit_data", None)
    out = {"parameter_names": list(model.parameter_names), "warnings": [],
           "interpretation": "Within-model absolute-sigma likelihood differences; no confidence interval or global-optimum guarantee."}
    try:
        if hasattr(model, "_prepare_data") and old_data is None:
            model._fit_data = model._prepare_data()
        baseline = _chi2(model, p)
        out["base_chi2"] = baseline
        for name, values in settings.items():
            rows = out[name] = []
            for target in values:
                row = {"target": float(target), "start": None, "parameters": None,
                       "success": False, "chi2": None, "delta_chi2": None,
                       "below_base_minimum": False, "nuisance_parameters": [],
                       "optimization_performed": False, "nfev": None,
                       "status": None, "message": ""}
                try:
                    q0, lower, upper, expand, labels = _constrained_coordinates(model, p, name, float(target))
                    row["start"] = expand(q0).tolist()
                    row["nuisance_parameters"] = labels
                    if len(q0):
                        def residual(q):
                            r = np.asarray(model.errFunc(expand(q)), dtype=float)
                            if not np.isfinite(r).all():
                                raise ValueError("Non-finite profile residual")
                            return r

                        row["optimization_performed"] = True
                        fit = least_squares(residual, q0, bounds=(lower, upper),
                                            jac="3-point", diff_step=1e-5, x_scale="jac",
                                            ftol=1e-9, xtol=1e-9, gtol=1e-9,
                                            max_nfev=model.config["init"].get("max_nfev", 500))
                        row.update(nfev=int(fit.nfev), status=int(fit.status), message=str(fit.message))
                        if not fit.success:
                            raise RuntimeError(f"Profile fit did not converge: {fit.message}")
                        solution = expand(fit.x)
                    else:
                        solution = expand(q0)
                        row.update(nfev=1, status=0, message="Direct evaluation: no nuisance parameters vary")
                    value = _chi2(model, solution)
                    delta = value-baseline
                    below = delta < -1e-8*max(1., abs(baseline))
                    row.update(success=True, parameters=solution.tolist(), chi2=value,
                               delta_chi2=delta, below_base_minimum=below)
                    if below:
                        out["warnings"].append(f"{name}={target:g} improves on the base chi2; the base fit is not the best solution found.")
                except Exception as exc:
                    row["message"] = f"{type(exc).__name__}: {exc}"
                    out["warnings"].append(f"Profile {name}={target:g} failed: {row['message']}")
                rows.append(row)
    finally:
        try:
            model.errFunc(p)
        finally:
            _restore(model, state)
            if hasattr(model, "_fit_data"):
                model._fit_data = old_data
    return out
