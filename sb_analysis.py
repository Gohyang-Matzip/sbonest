"""Optional multistart fits and constrained, nuisance-refitted likelihood profiles.

Residuals use supplied absolute intensity errors. Profile delta-chi2 values are
within-model likelihood differences, not confidence intervals or proof that the
base fit is a global minimum. All unsuccessful attempts are retained.
"""

import copy
import numbers

import numpy as np
from scipy.optimize import OptimizeResult, least_squares

from fit import _block_jacobian


class MultiStartError(RuntimeError):
    """All fits failed; ``attempts`` retains their starting vectors and errors."""

    def __init__(self, attempts):
        """Record every attempt so callers can save failure evidence."""
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


def _init_config(model):
    """A model's init section; Sideband models expand grouped-name aliases."""
    accessor = getattr(model, "init_config", None)
    return accessor() if callable(accessor) else model.config["init"]


def _json_safe(value):
    if isinstance(value, np.ndarray):
        return _json_safe(value.tolist())
    if isinstance(value, np.generic):
        return _json_safe(value.item())
    if isinstance(value, dict):
        return {key: _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def _initialization(model):
    names = ("initial_parameters", "free", "lower", "upper")
    return (_json_safe({name: getattr(model, name) for name in names})
            if all(hasattr(model, name) for name in names) else None)


def _read_initialization(values, size):
    if (not isinstance(values, dict)
            or not {"initial_parameters", "free", "lower", "upper"} <= set(values)):
        raise ValueError("Missing completed fit initialization")
    try:
        initial = np.asarray(values["initial_parameters"], dtype=float)
        free_values = values["free"]
        if (not isinstance(free_values, list) or not free_values
                or any(isinstance(i, bool) or not isinstance(i, int) for i in free_values)):
            raise ValueError("Invalid completed free parameter indices")
        free = np.asarray(free_values, dtype=int)
        lower = np.array([-np.inf if x is None else x for x in values["lower"]], dtype=float)
        upper = np.array([np.inf if x is None else x for x in values["upper"]], dtype=float)
    except (KeyError, TypeError, OverflowError) as exc:
        raise ValueError("Malformed completed fit initialization") from exc
    if (initial.shape != (size,) or not np.isfinite(initial).all()
            or lower.shape != initial.shape or upper.shape != initial.shape
            or np.isnan(lower).any() or np.isnan(upper).any() or np.any(lower >= upper)
            or len(set(free)) != len(free) or np.any(free < 0) or np.any(free >= size)):
        raise ValueError("Invalid completed fit initialization")
    return dict(initial_parameters=initial, free=free, lower=lower, upper=upper)


def snapshot_fit(model, p, covariance):
    """Return JSON-safe fit state, including the original start and uncertainty.

    Null covariance entries represent unavailable uncertainty; null lower/upper
    bounds represent negative/positive infinity, respectively.
    """
    return _json_safe({"schema_version": 1, "parameter_names": list(model.parameter_names),
                       "parameters": np.asarray(p), "covariance": np.asarray(covariance),
                       "state": _state(model)})


def _read_snapshot(model, snapshot):
    size = len(model.parameter_names)
    if (not isinstance(snapshot, dict) or snapshot.get("schema_version") != 1
            or snapshot.get("parameter_names") != list(model.parameter_names)):
        raise ValueError("Completed fit parameter names or snapshot version do not match")
    try:
        p = np.asarray(snapshot["parameters"], dtype=float)
        covariance = np.asarray(snapshot["covariance"], dtype=float)
        if p.shape != (size,) or not np.isfinite(p).all() or covariance.shape != (size, size):
            raise ValueError("Invalid completed fit parameters or covariance")
        state = copy.deepcopy(snapshot["state"])
        if not isinstance(state, dict) or set(state) != set(_STATE):
            raise ValueError("Completed fit statistics are missing or unknown")
        state.update(_read_initialization(state, size))
        vary = _init_config(model).get("vary", model.parameter_names)
        configured_free = [model.parameter_names.index(name) for name in vary]
        if not np.array_equal(state["free"], configured_free):
            raise ValueError("Completed fit free parameters do not match configured vary order")
        result = state["result"]
        if (not isinstance(result, dict) or result.get("success") is not True
                or not {"x", "fun", "jac", "grad", "active_mask", "nfev", "status", "message"} <= set(result)):
            raise ValueError("Completed fit has no converged optimizer result")
        for key in ("x", "fun", "jac", "grad", "active_mask"):
            if key in result:
                result[key] = np.asarray(result[key], dtype=int if key == "active_mask" else float)
        if (result["x"].shape != (len(state["free"]),)
                or result["fun"].ndim != 1
                or result["jac"].shape != (len(result["fun"]), len(state["free"]))
                or result["grad"].shape != result["x"].shape
                or result["active_mask"].shape != result["x"].shape
                or not all(np.isfinite(result[key]).all() for key in ("x", "fun", "jac"))):
            raise ValueError("Completed optimizer result has invalid dimensions or values")
        if (not np.array_equal(p[state["free"]], result["x"])
                or not _number(state["chi2"]) or state["chi2"] < 0
                or state["npar"] != len(state["free"])
                or state["nvar"] != len(result["fun"])
                or state["dof"] != state["nvar"] - state["npar"]
                or not isinstance(state["rank"], int)
                or not 0 <= state["rank"] <= len(state["free"])):
            raise ValueError("Completed fit statistics do not match its optimizer result")
        state["result"] = OptimizeResult(result)
        if state.get("condition") is None and "condition" in state:
            state["condition"] = np.inf
    except (KeyError, TypeError, OverflowError) as exc:
        raise ValueError("Malformed completed fit snapshot") from exc
    if (p.shape != (size,) or not np.isfinite(p).all() or covariance.shape != (size, size)
            or np.isinf(covariance).any() or np.any(p < state["lower"])
            or np.any(p > state["upper"])):
        raise ValueError("Invalid completed fit parameters or covariance")
    return p, covariance, state


def restore_fit(model, snapshot):
    """Restore predictions, optimizer statistics and uncertainty without fitting."""
    p, covariance, state = _read_snapshot(model, snapshot)
    original = _state(model)
    cache_names = ("_fit_data", "_initializing")
    caches = {name: getattr(model, name) for name in cache_names if hasattr(model, name)}
    accepted = False
    try:
        # Derive Sideband bounds from the actual configuration without a grid
        # search, initial overrides or optimization. Generic prepared models
        # expose their current bounds directly.
        if hasattr(model, "prepare_fit"):
            cfg = {key: value for key, value in _init_config(model).items()
                   if key not in ("initial", "multistart", "profile", "bootstrap", "profile_interval")}
            model.prepare_fit(p0=p, fitting_config=cfg)
        for name in ("lower", "upper"):
            if not hasattr(model, name) or not np.array_equal(state[name], getattr(model, name)):
                raise ValueError("Completed fit bounds do not match current model configuration")
        residual = np.asarray(model.errFunc(p), dtype=float)
        value = float(np.sum(residual**2))
        if (residual.ndim != 1 or not np.isfinite(residual).all()
                or not np.array_equal(residual, state["result"].fun)
                or state["nvar"] != len(residual)
                or state["npar"] != len(state["free"])
                or state["dof"] != len(residual)-len(state["free"])
                or not np.isclose(value, state["chi2"], rtol=1e-12, atol=1e-12)):
            raise ValueError("Completed fit residuals or statistics disagree with current data")
        _restore(model, state)
        accepted = True
    finally:
        if not accepted:
            for name in _STATE:
                if name not in original and hasattr(model, name):
                    delattr(model, name)
            _restore(model, original)
        for name in cache_names:
            if name in caches:
                setattr(model, name, caches[name])
            elif hasattr(model, name):
                delattr(model, name)
    return p, covariance


def _chi2(model, p):
    residual = np.asarray(model.errFunc(p), dtype=float)
    if not np.isfinite(residual).all():
        raise ValueError("Non-finite analysis residual")
    value = float(np.sum(residual**2))
    if not np.isfinite(value):
        raise ValueError("Non-finite analysis chi2")
    return value


def local_jacobian(model, p, *, pool=None):
    """Grouped finite-difference Jacobian of the Sideband residual at full vector p.

    Uses the same stencils, steps and bounds as ``SidebandModel.fit`` for the
    model's current ``free``, ``lower`` and ``upper``; columns follow ``free``.
    """
    p = np.asarray(p, dtype=float)
    free = np.asarray(model.free, dtype=int)
    lower, upper = np.asarray(model.lower), np.asarray(model.upper)
    if p.shape != (len(model.parameter_names),) or not np.isfinite(p).all():
        raise ValueError("Invalid parameter vector for a local Jacobian")
    if np.any(p[free] < lower[free]) or np.any(p[free] > upper[free]):
        raise ValueError("Jacobian point lies outside the parameter bounds")

    def expand(q):
        out = p.copy()
        out[free] = q
        return out

    def residual(q):
        values = np.asarray(model.errFunc(expand(q)), dtype=float)
        if not np.isfinite(values).all():
            raise ValueError("Non-finite residual at the Jacobian point")
        return values

    evaluate_many = None
    if pool is not None:
        def evaluate_many(vectors):
            return pool.evaluate_many([expand(q) for q in vectors])
    sizes = [sum(len(es.offset) for es in r.estSpecs) for r in model.dataset.res if r.active]
    n_global = getattr(model, "n_global", 2 + len(model.rf_names))
    nlocal = len(model.local_keys)
    full_order = np.array_equal(free, np.arange(len(p)))
    jacobian = _block_jacobian(
        residual, sizes, n_global, nlocal,
        relative_step=1e-5 if not full_order or model.proton_mode == "fit" else 1e-6,
        free=free, bounds=(lower[free], upper[free]),
        method="2-point" if full_order else "3-point", evaluate_many=evaluate_many)
    return jacobian(p[free])


def identifiability(model, p, *, pool=None, strong=0.95):
    """Local linear identifiability diagnostics at p with the supplied absolute sigma.

    Reports column-scaled singular values, rank, condition, the expected
    standard errors and relative errors from (J^T J)^-1, zero-sensitivity and
    weakly determined free parameters, strong pairwise correlations, and the
    derived kex/pB errors. These are properties of the point and the design, not
    of a fit; they neither guarantee nor replace convergence.
    """
    from sb_report import derived_errors

    names = list(model.parameter_names)
    free = [int(i) for i in model.free]
    p = np.asarray(p, dtype=float)
    cached = getattr(model, "_fit_data", None)
    try:
        if cached is None and hasattr(model, "_prepare_data"):
            model._fit_data = model._prepare_data()
        jac = local_jacobian(model, p, pool=pool)
        chi2 = _chi2(model, p)
    finally:
        if hasattr(model, "_fit_data"):
            model._fit_data = cached
    norms = np.linalg.norm(jac, axis=0)
    scaled = jac / np.where(norms > 0, norms, 1.0)
    singular = np.linalg.svd(scaled, compute_uv=False)
    rank = int(np.count_nonzero(singular > singular[0] * 1e-8)) if singular.size and singular[0] > 0 else 0
    condition = float(singular[0] / singular[-1]) if singular.size and singular[-1] > 0 else None
    covariance = np.full((len(names), len(names)), np.nan)
    if rank == len(free):
        _, s, vt = np.linalg.svd(jac, full_matrices=False)
        covariance[np.ix_(free, free)] = (vt.T / s**2) @ vt
    std = np.sqrt(np.diag(covariance))
    expected, relative, weak = {}, {}, []
    for i in free:
        value = float(std[i]) if np.isfinite(std[i]) else None
        expected[names[i]] = value
        ratio = None if value is None or p[i] == 0 else value / abs(p[i])
        relative[names[i]] = ratio
        if value is None or (ratio is not None and ratio > 1.0):
            weak.append(names[i])
    zero = [names[i] for k, i in enumerate(free) if norms[k] == 0]
    pairs = []
    for a in range(len(free)):
        for b in range(a + 1, len(free)):
            i, k = free[a], free[b]
            if np.isfinite(std[i]) and np.isfinite(std[k]) and std[i] > 0 and std[k] > 0:
                value = float(covariance[i, k] / (std[i] * std[k]))
                if abs(value) >= strong:
                    pairs.append({"parameters": [names[i], names[k]], "correlation": value})
    pairs.sort(key=lambda row: -abs(row["correlation"]))
    two_state = "kbc" not in names
    derived = (derived_errors(p, covariance) if two_state and rank == len(free) and {0, 1} <= set(free)
               else {"kex": None, "pB": None})
    return {
        "parameters": {name: float(p[i]) for i, name in enumerate(names)},
        "free_parameters": [names[i] for i in free],
        "n_points": int(jac.shape[0]), "n_free": len(free),
        "chi2_at_point": chi2,
        "scaled_singular_values": [float(v) for v in singular],
        "jacobian_rank": rank,
        "scaled_condition": condition,
        "expected_se": expected, "relative_se": relative,
        "zero_sensitivity": zero, "weak_parameters": weak,
        "strong_correlations": pairs, "strong_threshold": strong,
        "derived_se": derived,
        "interpretation": ("Local linear diagnostics at this point with the supplied absolute sigma; "
                           "expected errors assume the point is the solution. Weak means a relative error "
                           "above 100% or no finite error."),
    }


def _number(value):
    return isinstance(value, numbers.Real) and not isinstance(value, (bool, np.bool_)) and np.isfinite(value)


def _multistart_settings(model, settings):
    normalize = getattr(model, "normalize_init", None)
    if callable(normalize) and isinstance(settings, dict):
        settings = normalize({"multistart": settings})["multistart"]
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
    vary = _init_config(model).get("vary", model.parameter_names)
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


def base_fit_config(model):
    """The configured fit settings without optional analyses."""
    return {k: v for k, v in _init_config(model).items()
            if k not in ("multistart", "profile", "bootstrap", "profile_interval")}


def fit_attempt(model, start, source, index, cfg):
    """Fit one start and return its JSON-safe attempt row; never raises for fit failure.

    ``start`` is a full parameter vector or None for the configured initial fit.
    The same function serves the serial loop and pool workers, so rows agree exactly.
    """
    row = {"index": int(index), "source": source, "start": None,
           "parameters": None, "success": False, "chi2": None,
           "nfev": None, "status": None, "message": "", "selected": False,
           "initialization": None, "fit_snapshot": None}
    model.result = None
    prior_initial = getattr(model, "initial_parameters", None)
    if start is not None:
        start = np.asarray(start, dtype=float)
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
        row["fit_snapshot"] = snapshot_fit(model, p, covariance)
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
    if row["start"] is not None:
        row["initialization"] = _initialization(model)
    return row


def fit_multistart(model, settings, *, completed=None, on_complete=None, pool=None):
    """Fit the configured baseline plus starts; restore and return the best fit.

    ``model.fit`` must expose full ``initial_parameters``, ``lower``, and
    ``upper`` arrays before optimization. Arrays in result rows follow
    ``model.parameter_names``. All failed runs raise MultiStartError with the
    complete attempt list, allowing callers to save failure evidence. ``completed``
    is an ordered prefix of prior rows. ``on_complete(row)`` runs after each new
    attempt; persistence errors and interrupts propagate to the caller.
    With ``pool``, explicit and random starts are fitted in worker processes in
    index order; the configured initial fit always runs in this process.
    """
    starts, count, seed = _multistart_settings(model, settings)
    cfg = base_fit_config(model)
    completed = [] if completed is None else copy.deepcopy(completed)
    sources = ["initial"] + ["explicit"] * len(starts) + ["random"] * count
    if not isinstance(completed, list) or len(completed) > len(sources):
        raise ValueError("Completed multistart attempts must be a valid prefix")
    for index, row in enumerate(completed):
        if (not isinstance(row, dict) or row.get("index") != index
                or isinstance(row["index"], bool) or row.get("source") != sources[index]
                or not isinstance(row.get("success"), bool)
                or not isinstance(row.get("message"), str)
                or not {"nfev", "status", "parameters", "selected", "chi2", "fit_snapshot"} <= set(row)):
            raise ValueError("Completed multistart attempt does not match its configured index")
        initial = _read_initialization(row.get("initialization"), len(model.parameter_names))
        if row.get("start") != initial["initial_parameters"].tolist():
            raise ValueError("Completed multistart start does not match its initialization")
        if row["success"]:
            restored_p, _, restored_state = _read_snapshot(model, row["fit_snapshot"])
            if (not _number(row.get("chi2")) or row.get("parameters") != restored_p.tolist()
                    or not np.isclose(row["chi2"], restored_state["chi2"], rtol=1e-12, atol=0.)
                    or row["initialization"] != {key: row["fit_snapshot"]["state"][key]
                                                  for key in ("initial_parameters", "free", "lower", "upper")}):
                raise ValueError("Completed multistart row disagrees with its fit snapshot")
    attempts, best = [], None
    progress = None

    def accept(row):
        nonlocal best
        if row["success"] and (best is None or row["chi2"] < best[0]):
            best = row["chi2"], row["fit_snapshot"], row["index"]
        attempts.append(row)
        if on_complete is not None:
            on_complete(copy.deepcopy(row))
        if progress is not None:
            progress.step()

    def reuse(index, start):
        nonlocal best
        row = completed[index]
        if start is not None and row["start"] != np.asarray(start, dtype=float).tolist():
            raise ValueError("Completed multistart starting vector does not match configuration")
        row["selected"] = False
        _restore(model, _read_initialization(row["initialization"], len(model.parameter_names)))
        if row["success"] and (best is None or row["chi2"] < best[0]):
            best = row["chi2"], row["fit_snapshot"], index
        attempts.append(row)

    if completed:
        reuse(0, None)
    else:
        accept(fit_attempt(model, None, "initial", 0, cfg))
    if attempts[0]["start"] is None or not all(hasattr(model, key) for key in ("initial_parameters", "lower", "upper")):
        raise MultiStartError(attempts)
    base = np.asarray(model.initial_parameters).copy()
    planned = []
    for entry in starts:
        start = base.copy()
        for name, value in entry.items():
            start[model.parameter_names.index(name)] = value
        planned.append((start, "explicit"))
    rng = np.random.default_rng(seed)
    for _ in range(count):
        planned.append((_random_start(model, base, rng), "random"))
    pending = []
    for index, (start, source) in enumerate(planned, 1):
        if index < len(completed):
            reuse(index, start)
        else:
            pending.append((index, start, source))
    if pending:
        from sb_parallel import Progress

        progress = Progress("multistart", len(sources), enabled=bool(getattr(model, "verbose", False)),
                            done=len(attempts))
    if pending and pool is not None:
        from sb_parallel import multistart_task

        payloads = [{"start": start.tolist(), "source": source, "index": index, "cfg": cfg}
                    for index, start, source in pending]
        for row in pool.map(multistart_task, payloads):
            accept(row)
    else:
        for index, start, source in pending:
            accept(fit_attempt(model, start, source, index, cfg))
    if best is None:
        raise MultiStartError(attempts)
    _, snapshot, selected = best
    p, covariance = restore_fit(model, snapshot)
    attempts[selected]["selected"] = True
    return p, covariance, attempts


def _profile_settings(settings, names):
    if not isinstance(settings, dict) or set(settings) - {"kex", "pB", "v1n_scale"}:
        raise ValueError("profile allows kex, pB, and v1n_scale grids")
    for name, values in settings.items():
        if name == "v1n_scale" and name not in names:
            raise ValueError("v1n_scale profiling requires scale RF mode")
        if name in ("kex", "pB") and "kbc" in names:
            raise ValueError("kex and pB profiles are defined for the two-state model only")
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

    # A constrained exchange coordinate remains global even when its label and
    # position differ from either original rate. Locals keep full-vector indices.
    full_indices = [names.index("kab") if i is None else i for i in coordinates]
    return np.asarray(q0), np.asarray(lo), np.asarray(hi), expand, labels, full_indices


def _profile_jacobian(model, residual, full_indices, lower, upper, evaluate_many=None):
    """Group only Sideband's verified global/residue parameter and row layout."""
    from sbfit import SidebandModel

    if type(model) is not SidebandModel:
        return "3-point"
    active = [(i, r) for i, r in enumerate(model.dataset.res) if r.active]
    rate_names = list(getattr(model, "rate_names", ["kab", "kba"]))
    names = [*rate_names, *model.rf_names]
    names.extend(f"{r.label}.{key}" for _, r in active for key in model.local_keys)
    if list(model.parameter_names) != names:
        return "3-point"
    data = model._fit_data if model._fit_data is not None else model._prepare_data()
    sizes = [sum(len(row[2]) for row in data if row[0] == i) for i, _ in active]
    row_order = [row[0] for row in data]
    expected_order = [i for i, _ in active for row in data if row[0] == i]
    if row_order != expected_order or not all(sizes):
        return "3-point"
    grouped = _block_jacobian(residual, sizes, len(rate_names) + len(model.rf_names), len(model.local_keys),
                              relative_step=1e-5, free=full_indices,
                              bounds=(lower, upper), method="3-point",
                              evaluate_many=evaluate_many)
    # SciPy's dense 3-point result is column-major. Matching that layout keeps
    # column norms and gradient reductions identical even for weak H rates.
    return lambda q: np.asfortranarray(grouped(q))


def profile_point(model, p, name, target, baseline, *, pool=None):
    """Constrained nuisance refit for one target; returns its JSON-safe row.

    Requires the model's current bounds and free indices for the base fit and a
    cached ``_fit_data``. Shared by the serial scan and pool workers.
    ``pool`` evaluates the refit's Jacobian columns in worker processes.
    """
    row = {"target": float(target), "start": None, "parameters": None,
           "success": False, "chi2": None, "delta_chi2": None,
           "below_base_minimum": False, "nuisance_parameters": [],
           "optimization_performed": False, "nfev": None,
           "status": None, "message": ""}
    try:
        q0, lower, upper, expand, labels, full_indices = _constrained_coordinates(model, p, name, float(target))
        row["start"] = expand(q0).tolist()
        row["nuisance_parameters"] = labels
        if len(q0):
            memo = {}

            def residual(q):
                key = np.asarray(q, dtype=float).tobytes()
                if memo.get("key") == key:
                    return memo["values"].copy()
                r = np.asarray(model.errFunc(expand(q)), dtype=float)
                if not np.isfinite(r).all():
                    raise ValueError("Non-finite profile residual")
                memo["key"], memo["values"] = key, r.copy()
                return r

            row["optimization_performed"] = True
            evaluate_many = None
            if pool is not None:
                def evaluate_many(vectors):
                    return pool.evaluate_many([expand(q) for q in vectors])
            jacobian = _profile_jacobian(model, residual, full_indices, lower, upper, evaluate_many)
            fit = least_squares(residual, q0, bounds=(lower, upper),
                                jac=jacobian, diff_step=1e-5, x_scale="jac",
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
    except Exception as exc:
        row["message"] = f"{type(exc).__name__}: {exc}"
    return row


def _profile_warning(name, target, row):
    if not row["success"]:
        return f"Profile {name}={target:g} failed: {row['message']}"
    if row["below_base_minimum"]:
        return f"{name}={target:g} improves on the base chi2; the base fit is not the best solution found."
    return None


def profile_likelihood(model, p, settings, *, completed=None, on_complete=None, pool=None):
    """Nuisance-refit each requested target, preserving the original fitted state.

    Failed targets have null chi2/delta_chi2; negative successful differences
    remain negative and flag a base fit that this scan has improved upon.
    Every grid point starts from the supplied base fit after making its exact
    constraint feasible. No warm start or automatic confidence interval is used.
    ``completed`` maps profile names to row prefixes; ``on_complete(name, row)``
    runs after each new point. Callback exceptions propagate after state restore.
    With ``pool``, new points are refitted in worker processes in grid order.
    """
    _profile_settings(settings, model.parameter_names)
    p = np.asarray(p, dtype=float).copy()
    if p.shape != (len(model.parameter_names),) or not np.isfinite(p).all():
        raise ValueError("Invalid base profile parameter vector")
    completed = {} if completed is None else copy.deepcopy(completed)
    if not isinstance(completed, dict) or set(completed) - set(settings):
        raise ValueError("Completed profiles contain unconfigured names")
    for name, rows in completed.items():
        if not isinstance(rows, list) or len(rows) > len(settings[name]):
            raise ValueError("Completed profile rows must be a configured prefix")
        for target, row in zip(settings[name], rows):
            if (not isinstance(row, dict) or not _number(row.get("target"))
                    or row["target"] != target or not isinstance(row.get("success"), bool)
                    or not isinstance(row.get("message"), str)
                    or not isinstance(row.get("below_base_minimum"), bool)
                    or not isinstance(row.get("optimization_performed"), bool)
                    or not isinstance(row.get("nuisance_parameters"), list)
                    or not {"start", "parameters", "chi2", "delta_chi2", "nfev", "status"} <= set(row)):
                raise ValueError("Completed profile target does not match its grid")
            if row["success"]:
                parameters = np.asarray(row.get("parameters"), dtype=float)
                if (parameters.shape != p.shape or not np.isfinite(parameters).all()
                        or not _number(row.get("chi2")) or not _number(row.get("delta_chi2"))):
                    raise ValueError("Invalid completed successful profile row")
    state = _state(model)
    old_data = getattr(model, "_fit_data", None)
    out = {"parameter_names": list(model.parameter_names), "warnings": [],
           "interpretation": "Within-model absolute-sigma likelihood differences; no confidence interval or global-optimum guarantee."}
    try:
        if hasattr(model, "_prepare_data") and old_data is None:
            model._fit_data = model._prepare_data()
        baseline = _chi2(model, p)
        out["base_chi2"] = baseline
        total = sum(len(values) for values in settings.values())
        finished = sum(min(len(values), len(completed.get(name, []))) for name, values in settings.items())
        progress = None
        if finished < total:
            from sb_parallel import Progress

            progress = Progress("profile", total, enabled=bool(getattr(model, "verbose", False)), done=finished)
        for name, values in settings.items():
            rows = out[name] = []
            pending = []
            for index, target in enumerate(values):
                if index < len(completed.get(name, [])):
                    row = completed[name][index]
                    rows.append(row)
                    warning = _profile_warning(name, target, row)
                    if warning:
                        out["warnings"].append(warning)
                else:
                    pending.append(target)

            def accept(target, row):
                warning = _profile_warning(name, target, row)
                if warning:
                    out["warnings"].append(warning)
                rows.append(row)
                if on_complete is not None:
                    on_complete(name, copy.deepcopy(row))
                if progress is not None:
                    progress.step()

            if pending and pool is not None:
                from sb_parallel import profile_task

                cfg = base_fit_config(model)
                payloads = [{"p": np.asarray(p).tolist(), "name": name, "target": float(target),
                             "baseline": baseline, "cfg": cfg} for target in pending]
                for target, row in zip(pending, pool.map(profile_task, payloads)):
                    accept(target, row)
            else:
                for target in pending:
                    accept(target, profile_point(model, p, name, target, baseline))
    finally:
        try:
            model.errFunc(p)
        finally:
            _restore(model, state)
            if hasattr(model, "_fit_data"):
                model._fit_data = old_data
    return out


INTERVAL_NAMES = ("kex", "pB", "v1n_scale")


def _interval_settings(settings, names):
    """Validate init.profile_interval; return normalized settings."""
    allowed = {"parameters", "confidence", "max_evaluations", "relative_tolerance", "max_doublings"}
    if not isinstance(settings, dict) or set(settings) - allowed or "parameters" not in settings:
        raise ValueError("profile_interval requires parameters and allows confidence, max_evaluations, relative_tolerance, max_doublings")
    parameters = settings["parameters"]
    if (not isinstance(parameters, list) or not parameters or len(set(parameters)) != len(parameters)
            or any(name not in INTERVAL_NAMES for name in parameters)):
        raise ValueError("profile_interval.parameters must list unique names among kex, pB, v1n_scale")
    if "v1n_scale" in parameters and "v1n_scale" not in names:
        raise ValueError("v1n_scale intervals require scale RF mode")
    if ("kex" in parameters or "pB" in parameters) and "kbc" in names:
        raise ValueError("kex and pB intervals are defined for the two-state model only")
    confidence = settings.get("confidence", 0.95)
    if isinstance(confidence, bool) or not _number(confidence) or not 0 < confidence < 1:
        raise ValueError("profile_interval.confidence must be between 0 and 1")
    budget = settings.get("max_evaluations", 40)
    if isinstance(budget, bool) or not isinstance(budget, int) or budget < 4:
        raise ValueError("profile_interval.max_evaluations must be an integer of at least 4")
    tolerance = settings.get("relative_tolerance", 1e-3)
    if isinstance(tolerance, bool) or not _number(tolerance) or not 0 < tolerance < 0.5:
        raise ValueError("profile_interval.relative_tolerance must lie in (0, 0.5)")
    doublings = settings.get("max_doublings", 8)
    if isinstance(doublings, bool) or not isinstance(doublings, int) or doublings < 1:
        raise ValueError("profile_interval.max_doublings must be a positive integer")
    return {"parameters": list(parameters), "confidence": float(confidence),
            "max_evaluations": budget, "relative_tolerance": float(tolerance),
            "max_doublings": doublings}


def _interval_value(model, p, name):
    names = list(model.parameter_names)
    if name == "v1n_scale":
        return float(p[names.index(name)])
    total = p[names.index("kab")] + p[names.index("kba")]
    return float(total) if name == "kex" else float(p[names.index("kab")] / total)


def _interval_range(model, p, name):
    """Hard limits of the scanned quantity implied by the parameter bounds."""
    names = list(model.parameter_names)
    lower, upper = np.asarray(model.lower), np.asarray(model.upper)
    if name == "v1n_scale":
        i = names.index(name)
        return float(lower[i]), float(upper[i])
    a, b = names.index("kab"), names.index("kba")
    if name == "kex":
        return float(lower[a] + lower[b]), float(upper[a] + upper[b])
    return 0.0, 1.0


class ProfileIntervalError(RuntimeError):
    """A profile-interval search could not be completed on one side."""
    pass


def profile_intervals(model, p, covariance, settings, *, completed=None, on_complete=None, pool=None):
    """Likelihood-ratio intervals from the exact nuisance-refit profile.

    For each quantity the base value is bracketed outward in doubling steps from a
    local-SE-sized start and the crossing of the chi-square threshold is located
    with Brent's method on the profile delta chi2. Every evaluated point is
    retained in evaluation order. ``completed`` maps names to ordered evaluated
    rows that are replayed instead of refitted; the evaluation sequence is
    deterministic for an unchanged base fit, so a mismatch is reported.
    Intervals are within-model statements under the supplied absolute sigma.
    """
    from scipy.optimize import brentq
    from scipy.stats import chi2 as chi2_dist, norm
    from sb_report import derived_errors

    settings = _interval_settings(settings, model.parameter_names)
    p = np.asarray(p, dtype=float).copy()
    names = list(model.parameter_names)
    completed = {} if completed is None else copy.deepcopy(completed)
    if not isinstance(completed, dict) or set(completed) - set(settings["parameters"]):
        raise ValueError("Completed profile intervals contain unconfigured names")
    threshold = float(chi2_dist.ppf(settings["confidence"], 1))
    z = float(norm.ppf(0.5 + settings["confidence"] / 2))
    state = _state(model)
    old_data = getattr(model, "_fit_data", None)
    out = {"confidence": settings["confidence"], "threshold_delta_chi2": threshold,
           "warnings": [], "parameter_names": names,
           "interpretation": ("Likelihood-ratio intervals where the exact nuisance-refit profile crosses "
                              "the chi-square threshold; within-model statements under the supplied "
                              "absolute sigma, not guarantees of global optimality.")}
    try:
        if hasattr(model, "_prepare_data") and old_data is None:
            model._fit_data = model._prepare_data()
        baseline = _chi2(model, p)
        out["base_chi2"] = baseline
        std = np.sqrt(np.diag(np.asarray(covariance, dtype=float)))
        derived = derived_errors(p, np.asarray(covariance, dtype=float))
        for name in settings["parameters"]:
            estimate = _interval_value(model, p, name)
            local = derived.get(name) if name != "v1n_scale" else (
                float(std[names.index(name)]) if np.isfinite(std[names.index(name)]) else None)
            hard = _interval_range(model, p, name)
            result = {"estimate": estimate, "local_se": local, "lower": None, "upper": None,
                      "lower_message": "", "upper_message": "", "points": [], "success": False,
                      "message": "", "evaluations": 0,
                      "search_range": [float(v) if np.isfinite(v) else None for v in hard]}
            out[name] = result
            rows = list(completed.get(name, []))
            if not isinstance(rows, list):
                raise ValueError("Completed interval rows must be a list")
            cache = {}

            def evaluate(target):
                target = float(target)
                if target in cache:
                    return cache[target]
                index = result["evaluations"]
                if index < len(rows):
                    row = rows[index]
                    if not isinstance(row, dict) or row.get("target") != target:
                        raise ValueError(f"Completed profile_interval rows for {name} do not match the evaluation sequence")
                else:
                    row = profile_point(model, p, name, target, baseline, pool=pool)
                    if on_complete is not None:
                        on_complete(name, copy.deepcopy(row))
                result["evaluations"] += 1
                result["points"].append(row)
                if not row["success"]:
                    raise ProfileIntervalError(f"profile refit failed at {name}={target:g}: {row['message']}")
                if row["below_base_minimum"]:
                    out["warnings"].append(f"{name}={target:g} lies below the base chi2; the base fit is not the best solution found and this interval is unreliable.")
                cache[target] = float(row["delta_chi2"]) - threshold
                return cache[target]

            if local is None or not np.isfinite(local) or local <= 0:
                step = 0.1 * abs(estimate) if estimate else 0.01
                out["warnings"].append(f"No finite local standard error for {name}; bracketing starts from a 10% step.")
            else:
                step = z * local
            try:
                for side, sign in (("lower", -1.0), ("upper", 1.0)):
                    limit = hard[0] if sign < 0 else hard[1]
                    inner = estimate
                    found = None
                    for k in range(settings["max_doublings"] + 1):
                        target = estimate + sign * step * 2**k
                        clipped = False
                        if (sign < 0 and target <= limit) or (sign > 0 and target >= limit):
                            if not np.isfinite(limit):
                                raise ProfileIntervalError(f"{side} bracket for {name} exceeded the search range")
                            target, clipped = limit, True
                            if name == "kex" and target <= 0:
                                target = np.nextafter(0.0, 1.0)
                        if result["evaluations"] >= settings["max_evaluations"]:
                            raise ProfileIntervalError(f"{side} bracket for {name} exceeded max_evaluations")
                        value = evaluate(target)
                        if value >= 0:
                            found = (inner, target) if sign > 0 else (target, inner)
                            break
                        inner = target
                        if clipped:
                            result[f"{side}_message"] = f"profile stays below the threshold up to the bound {limit:g}"
                            break
                    if found is None:
                        if not result[f"{side}_message"]:
                            result[f"{side}_message"] = f"no threshold crossing within {settings['max_doublings']} doublings"
                        continue
                    lo, hi = found
                    remaining = settings["max_evaluations"] - result["evaluations"]
                    if remaining < 2:
                        raise ProfileIntervalError(f"{side} root finding for {name} exceeded max_evaluations")
                    xtol = settings["relative_tolerance"] * max(abs(estimate), np.finfo(float).tiny)
                    root = brentq(evaluate, lo, hi, xtol=xtol, rtol=4 * np.finfo(float).eps,
                                  maxiter=remaining, full_output=False, disp=True)
                    result[side] = float(root)
                    result[f"{side}_message"] = "threshold crossing located"
                result["success"] = result["lower"] is not None and result["upper"] is not None
                if not result["success"]:
                    result["message"] = "; ".join(m for m in (result["lower_message"], result["upper_message"]) if m)
                    out["warnings"].append(f"Profile interval for {name} is one-sided or open: {result['message']}")
                else:
                    result["message"] = "both crossings located"
            except (ProfileIntervalError, RuntimeError, ValueError) as exc:
                if isinstance(exc, ValueError) and "do not match the evaluation sequence" in str(exc):
                    raise
                result["message"] = f"{type(exc).__name__}: {exc}"
                out["warnings"].append(f"Profile interval for {name} failed: {result['message']}")
    finally:
        try:
            model.errFunc(p)
        finally:
            _restore(model, state)
            if hasattr(model, "_fit_data"):
                model._fit_data = old_data
    return out
