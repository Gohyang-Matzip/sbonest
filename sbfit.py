#!/usr/bin/env python3
"""SBONEST: fit all CEST points, including decoupling sidebands, with OC.

Usage: python sbfit.py config.json [--no-pdf]
ONEST text inputs and reports are retained; see SIDEBAND.md for configuration.
"""

from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from estmodel import b1_weights, est_model
from fit import PARAM_LAYOUT, _block_jacobian, generate_initial_parameters
from run import set_residue_flags
from sideband import composite_segments, profile, waveform_segments
from sb_report import derived_errors, prediction_rows
from sb_analysis import (
    _multistart_settings, _number, _profile_settings,
)
from sideband import profile_states, stationary_populations

# Sideband model name -> (ONEST parameter layout, number of exchanging sites).
SIDEBAND_METHODS = {
    "Sideband": ("Matrix", 2),
    "Sideband_3st_Linear": ("Matrix_3st_Linear", 3),
    "Sideband_3st_Triangle": ("Matrix_3st_Triangle", 3),
}
# Rates that must stay positive so every site keeps a return path to site A.
_RETURN_RATES = ("kba", "kcb")


def is_sideband_method(name):
    """True for init.Method values handled by SidebandModel."""
    return name in SIDEBAND_METHODS


def _key_base(name):
    """'A1.R2a[1]' -> 'R2a'; 'kab' -> 'kab'."""
    local = name.rsplit(".", 1)[-1]
    return local.split("[", 1)[0]


def _rate_lower(name):
    return 1e-8 if name in _RETURN_RATES else 0.0


def _known(config, keys, label):
    if not isinstance(config, dict) or set(config) - set(keys):
        raise ValueError(f"Invalid {label} keys; allowed: {', '.join(keys)}")


def _validate_fit_config(model, cfg):
    """Validate settings before the initial grid or an optimizer can run."""
    cfg = model.normalize_init(cfg)
    _known(cfg, ("Method", "kex", "pB", "initial", "bounds", "vary", "max_nfev",
                 "multistart", "profile", "bootstrap", "profile_interval"), "init")
    model.selMethod(cfg)
    for name in ("kex", "pB"):
        if name not in cfg:
            continue
        grid = cfg[name]
        _known(grid, ("min", "max", "nsteps"), f"init.{name}")
        lo, hi, steps = grid.get("min"), grid.get("max"), grid.get("nsteps", 6)
        if (not _number(lo) or not _number(hi) or lo > hi
                or isinstance(steps, bool) or not isinstance(steps, int) or steps <= 0
                or (name == "kex" and lo <= 0)
                or (name == "pB" and (lo < 0 or hi > 1))):
            raise ValueError(f"Invalid init.{name} grid limits or nsteps")
    initial = cfg.get("initial", {})
    _known(initial, model.parameter_names, "init.initial")
    if not all(_number(value) for value in initial.values()):
        raise ValueError("init.initial values must be finite numbers")
    bounds = cfg.get("bounds", {})
    _known(bounds, model.parameter_names, "init.bounds")
    for name, limits in bounds.items():
        if (not isinstance(limits, (list, tuple)) or len(limits) != 2
                or any(not _number(value) and value not in (-np.inf, np.inf)
                       for value in limits)
                or any(isinstance(value, (bool, np.bool_)) for value in limits)
                or limits[0] >= limits[1]):
            raise ValueError(f"Bounds for {name} must be two strictly ordered numbers")
        physical = model.physical_bounds(name)
        if limits[0] < physical[0] or limits[1] > physical[1]:
            raise ValueError(f"bounds for {name} must stay inside physical/v1n bounds")
    vary = cfg.get("vary", model.parameter_names)
    if (not isinstance(vary, list) or not vary
            or not all(isinstance(name, str) for name in vary)
            or len(set(vary)) != len(vary)
            or any(name not in model.parameter_names for name in vary)):
        raise ValueError("init.vary must be a nonempty list of unique known parameter names")
    budget = cfg.get("max_nfev", 500)
    if isinstance(budget, bool) or not isinstance(budget, int) or budget <= 0:
        raise ValueError("init.max_nfev must be a positive integer")
    if "multistart" in cfg:
        _multistart_settings(model, cfg["multistart"])
    if "profile" in cfg:
        _profile_settings(cfg["profile"], model.parameter_names)
    if "profile_interval" in cfg:
        from sb_analysis import _interval_settings

        _interval_settings(cfg["profile_interval"], model.parameter_names)
    if "bootstrap" in cfg:
        from sb_bootstrap import validate_bootstrap

        validate_bootstrap(cfg["bootstrap"])


class SidebandModel(est_model):
    """Reuse ONEST data, residue ordering, residuals and PDF/text reporting."""

    def __init__(self, config, config_dir="."):
        """Build the model from a Sideband configuration; dataset paths resolve from config_dir."""
        super().__init__()
        self.sideband_method = config["init"].get("Method", "Sideband")
        if self.sideband_method not in SIDEBAND_METHODS:
            raise ValueError("init.Method must be Sideband, Sideband_3st_Linear or Sideband_3st_Triangle")
        self.method, self.states = SIDEBAND_METHODS[self.sideband_method]
        self.rate_names = list(PARAM_LAYOUT[self.method][0])
        self.programName = "SBONEST 1.0 — exact NH Sideband model (OC)"
        self.config = config
        self._initializing = False
        self.result = None
        self.pool = None
        self._last_covariance = None
        sb = config["sideband"]
        _known(
            sb,
            ("decoupling", "datasets", "residues", "v1n", "proton_relaxation", "nitrogen_relaxation"),
            "sideband",
        )
        self.rf_config = sb.get("v1n", {"mode": "fixed"})
        _known(self.rf_config, ("mode", "initial", "bounds"), "v1n")
        self.rf_mode = self.rf_config.get("mode", "fixed")
        if self.rf_mode not in ("fixed", "scale", "per_dataset"):
            raise ValueError("v1n.mode must be fixed, scale, or per_dataset")
        for index, path in enumerate(config["datasets"]):
            before = {id(es) for r in self.dataset.res for es in r.estSpecs}
            self.dataset.addData(str(Path(config_dir) / path))
            for r in self.dataset.res:
                new = [es for es in r.estSpecs if id(es) not in before]
                if len(new) > 1:
                    raise ValueError(f"Duplicate residue {r.label} in dataset {index}")
                for es in new:
                    es.dataset_index = index
        error = set_residue_flags(self.dataset, config["residues"])
        if error:
            raise ValueError(error)
        active = [r for r in self.dataset.res if r.active]
        if not active:
            raise ValueError("Select at least one residue")
        self._prepare_data()
        for r in active:
            if not r.estSpecs:
                raise ValueError(f"No spectra for active residue {r.label}")
            for es in r.estSpecs:
                if (
                    es.T <= 0
                    or es.field <= 0
                    or es.v1 <= 0
                    or es.v1err < 0
                    or np.any(np.asarray(es.intstd) <= 0)
                ):
                    raise ValueError(
                        "Sideband fits require T, field, v1, intensity errors > 0; "
                        "v1err >= 0"
                    )
        self.h_shifts = sb["residues"]
        shift_keys = ("h_ppm_a", "h_ppm_b") + (("h_ppm_c",) if self.states == 3 else ())
        for r in active:
            h = self.h_shifts[r.label]
            _known(h, shift_keys, f"sideband.residues.{r.label}")
            if set(shift_keys) - set(h) or not np.isfinite([h[key] for key in shift_keys]).all():
                raise ValueError(
                    f"Supply finite {', '.join(shift_keys)} for every active residue"
                )
        n = len(config["datasets"])
        overrides = sb.get("datasets", [{} for _ in range(n)])
        if not isinstance(overrides, list) or len(overrides) != n:
            raise ValueError("sideband.datasets must align with the datasets list")
        explicit_h = any(
            key in entry
            for entry in [sb["decoupling"], *overrides]
            for key in ("R1H", "R2H")
        )
        h_config = sb.get("proton_relaxation", {})
        _known(h_config, ("mode",), "proton_relaxation")
        self.proton_mode = h_config.get("mode", "fixed" if explicit_h else "fit")
        if self.proton_mode not in ("fit", "fixed"):
            raise ValueError("proton_relaxation.mode must be fit or fixed")
        if self.proton_mode == "fit" and explicit_h:
            raise ValueError(
                "Remove decoupling R1H/R2H for automatic peakwise fitting; "
                "optional starting values belong in init.initial as <peak>.R1H/R2H"
            )
        n_config = sb.get("nitrogen_relaxation", {})
        _known(n_config, ("mode",), "nitrogen_relaxation")
        self.nitrogen_mode = n_config.get("mode", "shared")
        if self.nitrogen_mode not in ("shared", "per_field"):
            raise ValueError("nitrogen_relaxation.mode must be shared or per_field")
        self.decoupling, self.segments = [], []
        for entry in overrides:
            d = {**sb["decoupling"], **entry}
            _known(
                d,
                (
                    "h_larmor_mhz",
                    "h_carrier_ppm",
                    "p90_s",
                    "cycle",
                    "b1_scale",
                    "J_hz",
                    "R1H",
                    "R2H",
                    "waveform_json",
                    "rf_hz",
                ),
                "decoupling",
            )
            d = {"b1_scale": 1.0, "J_hz": 92.0, "R1H": 2.0, "R2H": 25.0, **d}
            required = [
                d["h_larmor_mhz"],
                d["h_carrier_ppm"],
                d["b1_scale"],
                d["J_hz"],
                d["R1H"],
                d["R2H"],
            ]
            if (
                not np.isfinite(required).all()
                or d["h_larmor_mhz"] <= 0
                or d["b1_scale"] <= 0
                or min(d["R1H"], d["R2H"]) < 0
            ):
                raise ValueError("Invalid decoupling field, scale, or relaxation")
            if "waveform_json" in d:
                if "p90_s" in d or "cycle" in d:
                    raise ValueError("Supply either waveform_json or p90_s/cycle")
                seg = waveform_segments(
                    Path(config_dir) / d["waveform_json"], d.get("rf_hz"), d["b1_scale"]
                )
            else:
                if "rf_hz" in d:
                    raise ValueError("rf_hz is only for a.u. OC waveforms")
                seg = composite_segments(
                    d["p90_s"], d.get("cycle", "RR"), d["b1_scale"]
                )
            self.decoupling.append(d)
            self.segments.append(seg)
        # Field groups: datasets sharing one proton Larmor frequency. Group-specific
        # relaxation parameters carry a "[g]" suffix only when several groups exist.
        self.field_groups = sorted({float(d["h_larmor_mhz"]) for d in self.decoupling})
        self.dataset_group = [self.field_groups.index(float(d["h_larmor_mhz"])) for d in self.decoupling]
        groups = len(self.field_groups)
        self.proton_groups = groups if self.proton_mode == "fit" and groups > 1 else 1
        self.nitrogen_groups = groups if self.nitrogen_mode == "per_field" and groups > 1 else 1
        nominal = np.asarray(self.dataset.v1s)
        self.rf_names = (
            ["v1n_scale"]
            if self.rf_mode == "scale"
            else [f"v1n[{i}]" for i in range(n)]
            if self.rf_mode == "per_dataset"
            else []
        )
        default = [1.0] if self.rf_mode == "scale" else nominal if self.rf_names else []
        self.rf_initial = np.atleast_1d(self.rf_config.get("initial", default)).astype(
            float
        )
        limits = (
            [0.5, 1.5]
            if self.rf_mode == "scale"
            else np.column_stack((nominal * 0.5, nominal * 1.5))
            if self.rf_names
            else []
        )
        self.rf_bounds = np.asarray(
            self.rf_config.get("bounds", limits), dtype=float
        ).reshape(-1, 2)
        if (
            self.rf_initial.shape != (len(self.rf_names),)
            or self.rf_bounds.shape != (len(self.rf_names), 2)
            or not np.isfinite(self.rf_initial).all()
            or not np.isfinite(self.rf_bounds).all()
            or np.any(self.rf_bounds[:, 0] <= 0)
            or np.any(self.rf_bounds[:, 0] >= self.rf_bounds[:, 1])
            or np.any(self.rf_initial < self.rf_bounds[:, 0])
            or np.any(self.rf_initial > self.rf_bounds[:, 1])
        ):
            raise ValueError(
                "v1n initial/bounds need one positive bounded value per RF parameter"
            )
        self.shift_keys = ("peak_ppm", "dw_ppm") + (("dwC_ppm",) if self.states == 3 else ())
        self.n_relax_keys = ("R1", "R2a", "R2b") + (("R2c",) if self.states == 3 else ())
        self.h_keys = ("R1H", "R2H")

        def suffix(group, count):
            return f"[{group}]" if count > 1 else ""

        self.local_keys = tuple(self.shift_keys) + tuple(
            f"{key}{suffix(g, self.nitrogen_groups)}"
            for g in range(self.nitrogen_groups) for key in self.n_relax_keys
        )
        if self.proton_mode == "fit":
            self.local_keys += tuple(
                f"{key}{suffix(g, self.proton_groups)}"
                for g in range(self.proton_groups) for key in self.h_keys
            )
        self.n_global = len(self.rate_names) + len(self.rf_names)
        self.parameter_names = (
            self.rate_names
            + self.rf_names
            + [f"{r.label}.{key}" for r in active for key in self.local_keys]
        )
        self.free = np.arange(len(self.parameter_names))

    def init_config(self):
        """The live init section with ungrouped aliases expanded to grouped names."""
        return self.normalize_init(self.config["init"])

    def expand_name(self, name):
        """Actual parameter names addressed by ``name``: itself, or all groups of an alias."""
        if name in self.parameter_names:
            return [name]
        if not isinstance(name, str) or "." not in name or "[" in name:
            return []
        return [actual for actual in self.parameter_names
                if "[" in actual and actual.split("[", 1)[0] == name]

    def normalize_init(self, cfg):
        """Copy of an init section with ungrouped aliases expanded to grouped names."""
        if not isinstance(cfg, dict):
            return cfg
        out = dict(cfg)
        for key in ("initial", "bounds"):
            if isinstance(cfg.get(key), dict):
                mapping = {}
                for name, value in cfg[key].items():
                    for target in self.expand_name(name) or [name]:
                        mapping[target] = value
                out[key] = mapping
        if isinstance(cfg.get("vary"), list):
            out["vary"] = [target for name in cfg["vary"] for target in (self.expand_name(name) or [name])]
        multistart = cfg.get("multistart")
        if isinstance(multistart, dict) and isinstance(multistart.get("starts"), list):
            starts = []
            for start in multistart["starts"]:
                if isinstance(start, dict):
                    expanded = {}
                    for name, value in start.items():
                        for target in self.expand_name(name) or [name]:
                            expanded[target] = value
                    starts.append(expanded)
                else:
                    starts.append(start)
            out["multistart"] = {**multistart, "starts": starts}
        return out

    def physical_bounds(self, name):
        """Hard limits of one parameter: rates, RF, shifts or relaxation."""
        if name in self.rate_names:
            return (_rate_lower(name), np.inf)
        if name in self.rf_names:
            return tuple(self.rf_bounds[self.rf_names.index(name)])
        if _key_base(name) in self.shift_keys:
            return (-np.inf, np.inf)
        return (0.0, np.inf)

    def selMethod(self, initConf):
        """Reject init sections whose Method differs from the model this instance was built for."""
        if initConf.get("Method", "Sideband") != self.sideband_method:
            raise ValueError(f"This model was built for init.Method = {self.sideband_method}")

    def _expand_initial(self, base):
        """ONEST-layout vector (rates, then shifts and one relaxation set per residue)
        to the Sideband layout with RF parameters, per-group relaxation and H rates."""
        base = np.asarray(base, dtype=float)
        n_rates, n_shift, n_relax = len(self.rate_names), len(self.shift_keys), len(self.n_relax_keys)
        local = base[n_rates:].reshape(-1, n_shift + n_relax)
        rows = []
        for row in local:
            expanded = [row[:n_shift]] + [row[n_shift:]] * self.nitrogen_groups
            if self.proton_mode == "fit":
                expanded += [np.array([2.0, 25.0])] * self.proton_groups
            rows.append(np.concatenate(expanded))
        return np.r_[base[:n_rates], self.rf_initial, np.concatenate(rows) if rows else []]

    def _nitrogen_parameters(self, p):
        """Sideband layout to the ONEST layout using group-0 nitrogen relaxation."""
        n_rates, n_shift, n_relax = len(self.rate_names), len(self.shift_keys), len(self.n_relax_keys)
        local = np.asarray(p[self.n_global:]).reshape(-1, len(self.local_keys))
        return np.r_[p[:n_rates], local[:, : n_shift + n_relax].ravel()]

    def seParam(self, p):
        """Unpack a full parameter vector into the ONEST-style dictionary plus RF and relaxation groups."""
        p = np.asarray(p, dtype=float)
        if self._initializing:
            p = self._expand_initial(p)
        if p.shape != (len(self.parameter_names),):
            raise ValueError("Wrong Sideband parameter vector length")
        P = super().seParam(self._nitrogen_parameters(p))
        P["v1n"] = p[len(self.rate_names) : self.n_global]
        local = p[self.n_global :].reshape(-1, len(self.local_keys))
        n_shift, n_relax = len(self.shift_keys), len(self.n_relax_keys)
        active = [r for r in self.dataset.res if r.active]
        P["nitrogen_rates"] = {
            r.label: [row[n_shift + g * n_relax : n_shift + (g + 1) * n_relax]
                      for g in range(self.nitrogen_groups)]
            for r, row in zip(active, local)
        }
        if self.proton_mode == "fit":
            start = n_shift + self.nitrogen_groups * n_relax
            P["proton_rates"] = {
                r.label: [row[start + 2 * g : start + 2 * g + 2] for g in range(self.proton_groups)]
                for r, row in zip(active, local)
            }
        return P

    def exchange_matrix(self, P):
        """Site-to-site rate matrix K[i, j] (i -> j) for the configured model."""
        if self.states == 2:
            return np.array([[0.0, P["kab"]], [P["kba"], 0.0]])
        K = np.array([[0.0, P["kab"], 0.0], [P["kba"], 0.0, P["kbc"]], [0.0, P["kcb"], 0.0]])
        if self.sideband_method == "Sideband_3st_Triangle":
            K[2, 0], K[0, 2] = P["kca"], P["kac"]
        return K

    def rf_values(self, P, es):
        """Actual nitrogen RF amplitude and its spread (Hz) for spectrum es under the RF mode."""
        if self.rf_mode == "scale":
            scale = P["v1n"][0]
        elif self.rf_mode == "per_dataset":
            scale = P["v1n"][es.dataset_index] / es.v1
        else:
            scale = 1.0
        # Keep fractional RF inhomogeneity fixed, independently of fit uncertainty.
        return es.v1 * scale, es.v1err * scale

    def errFunc(self, p_flat):
        """Residuals as in ONEST; with a worker pool the data blocks are predicted in parallel.

        Each block's prediction is the same ``calc`` call a worker would make
        serially, and the residual arithmetic stays in this process, so the
        values are identical to the serial path.
        """
        pool = getattr(self, "pool", None)
        if pool is None or self._initializing:
            return super().errFunc(p_flat)
        try:
            P = self.seParam(p_flat)
        except ValueError:
            return super().errFunc(p_flat)
        data = self._fit_data if self._fit_data is not None else self._prepare_data()
        if len(data) < 2:
            return super().errFunc(p_flat)
        del P
        predictions = pool.predict_blocks(p_flat, [len(row[2]) for row in data])
        residuals_list = []
        for (_, _, _, observed, std), est_calc in zip(data, predictions):
            est_calc = np.nan_to_num(est_calc, nan=1e6, posinf=1e6, neginf=-1e6)
            residuals_list.append((observed - est_calc) / std)
        residuals = np.concatenate(residuals_list)
        self.chi2 = np.sum(residuals**2)
        self.npar = len(p_flat)
        self.nvar = len(residuals)
        self.dof = max(1, self.nvar - self.npar)
        return residuals

    def calc(self, P, i, dRF, es):
        """Predicted normalized intensities of residue i at offsets dRF (ppm) for spectrum es."""
        d = self.decoupling[es.dataset_index]
        label = self.dataset.res[i].label
        h = self.h_shifts[label]
        group = self.dataset_group[es.dataset_index]
        r1h, r2h = (
            P["proton_rates"][label][group if self.proton_groups > 1 else 0]
            if self.proton_mode == "fit"
            else (d["R1H"], d["R2H"])
        )
        rates = P["nitrogen_rates"][label][group if self.nitrogen_groups > 1 else 0]
        fields, weights = b1_weights(*self.rf_values(P, es))
        out = np.zeros_like(np.asarray(dRF, dtype=float))
        offsets = (np.asarray(dRF) - P["dGs"][i]) * es.field
        carrier, larmor = d["h_carrier_ppm"], d["h_larmor_mhz"]
        for nu, weight in zip(fields / (2 * np.pi), weights):
            if self.states == 2:
                out += weight * profile(
                    self.segments[es.dataset_index],
                    offsets,
                    T=es.T,
                    nu=abs(nu),
                    kab=P["kab"],
                    kba=P["kba"],
                    dw=P["dws"][i] * es.field,
                    r1=rates[0],
                    r2a=rates[1],
                    r2b=rates[2],
                    ha=(h["h_ppm_a"] - carrier) * larmor,
                    hb=(h["h_ppm_b"] - carrier) * larmor,
                    J=d["J_hz"],
                    r1h=r1h,
                    r2h=r2h,
                )
            else:
                out += weight * profile_states(
                    self.segments[es.dataset_index],
                    offsets,
                    T=es.T,
                    nu=abs(nu),
                    exchange=self.exchange_matrix(P),
                    shifts=[0.0, P["dws"][i] * es.field, P["dwCs"][i] * es.field],
                    h_shifts=[(h[key] - carrier) * larmor for key in ("h_ppm_a", "h_ppm_b", "h_ppm_c")],
                    r1=rates[0],
                    r2=[rates[1], rates[2], rates[3]],
                    J=d["J_hz"],
                    r1h=r1h,
                    r2h=r2h,
                )
        return out if np.ndim(dRF) else float(out)

    def prepare_fit(self, p0=None, fitting_config=None):
        """Prepare current data, vectors and bounds without running an optimizer.

        Return a full starting vector. Publish its bounds before checking whether
        it lies inside them, so explicit multistart attempts can repair a rejected
        baseline. Each call clears stale fit state and releases the data cache.
        """
        self.result = None
        self._fit_data = None
        self._initializing = False
        for name in ("initial_parameters", "lower", "upper", "rank", "condition"):
            self.__dict__.pop(name, None)
        self.free = np.arange(len(self.parameter_names))
        self.chi2, self.nvar, self.npar, self.dof = 0.0, 0, 0, 0
        cfg = self.normalize_init(fitting_config or self.config["init"])
        _validate_fit_config(self, cfg)
        self._fit_data = self._prepare_data()
        try:
            if p0 is None:
                self._initializing = True
                try:
                    base = generate_initial_parameters(self, cfg)
                finally:
                    self._initializing = False
                p0 = self._expand_initial(base)
            p0 = np.asarray(p0, dtype=float).copy()
            if p0.shape != (len(self.parameter_names),) or not np.isfinite(p0).all():
                raise ValueError("Invalid initial parameter vector")
            nlocal = len(self.local_keys)
            n_rates = len(self.rate_names)
            lower = np.r_[
                [_rate_lower(name) for name in self.rate_names],
                self.rf_bounds[:, 0],
                np.tile(
                    [-np.inf if key in self.shift_keys else 0.0 for key in self.local_keys],
                    (len(p0) - self.n_global) // nlocal,
                ),
            ]
            upper = np.r_[
                np.full(n_rates, np.inf),
                self.rf_bounds[:, 1],
                np.full(len(p0) - self.n_global, np.inf),
            ]
            for name, value in cfg.get("initial", {}).items():
                p0[self.parameter_names.index(name)] = value
            for name, limits in cfg.get("bounds", {}).items():
                k = self.parameter_names.index(name)
                lo, hi = limits
                if lo < lower[k] or hi > upper[k]:
                    raise ValueError(
                        f"bounds for {name} must stay inside physical/v1n bounds"
                    )
                lower[k], upper[k] = lo, hi
            vary = cfg.get("vary", self.parameter_names)
            if not isinstance(vary, list) or not vary or len(set(vary)) != len(vary):
                raise ValueError(
                    "init.vary must be a nonempty list of unique parameter names"
                )
            self.free = np.array([self.parameter_names.index(n) for n in vary])
            if self.rf_mode == "per_dataset":
                used = {es.dataset_index for _, es, *_ in self._fit_data}
                for i in range(len(self.rf_names)):
                    if n_rates + i in self.free and i not in used:
                        raise ValueError(f"No active data for varying v1n[{i}]")
            if (
                np.isnan(lower).any()
                or np.isnan(upper).any()
                or np.any(lower >= upper)
            ):
                raise ValueError("Parameter bounds must be valid and strictly ordered")
            ndata = sum(len(x[2]) for x in self._fit_data)
            if ndata <= len(self.free):
                raise ValueError(
                    "More data points than varying parameters are required"
                )

            self.initial_parameters = p0.copy()
            self.lower, self.upper = lower.copy(), upper.copy()
            if (not np.isfinite(p0).all()
                    or np.any(p0 < lower) or np.any(p0 > upper)):
                raise ValueError("Initial parameters must be finite and inside valid bounds")

            return p0.copy()
        finally:
            self._fit_data = None

    def fit(self, p0=None, fitting_config=None):
        """Least-squares fit from p0 (or the configured grid start); returns (parameters, covariance)."""
        cfg = self.normalize_init(fitting_config or self.config["init"])
        p0 = self.prepare_fit(p0, fitting_config)
        lower, upper = self.lower, self.upper
        nlocal = len(self.local_keys)
        self._fit_data = self._prepare_data()
        try:
            def expand(q):
                p = p0.copy()
                p[self.free] = q
                return p

            memo = {}

            def residual(q):
                # SciPy evaluates the residual at an accepted step and then asks for
                # the Jacobian at the same point; reuse that vector instead of
                # recomputing it. The cached values are exactly the recomputed ones.
                key = np.asarray(q, dtype=float).tobytes()
                if memo.get("key") == key:
                    return memo["values"].copy()
                values = self.errFunc(expand(q))
                if not np.isfinite(values).all():
                    raise ValueError("Non-finite Sideband residual")
                memo["key"], memo["values"] = key, np.asarray(values, dtype=float).copy()
                return values

            sizes = [
                sum(len(es.offset) for es in r.estSpecs)
                for r in self.dataset.res
                if r.active
            ]
            full_order = np.array_equal(self.free, np.arange(len(p0)))
            pool = getattr(self, "pool", None)
            evaluate_many = None
            if pool is not None:
                def evaluate_many(vectors):
                    # Workers hold identical models; results are byte-identical to residual().
                    # A single vector (the Jacobian base point) uses the memo or the
                    # block-parallel residual instead of occupying one worker.
                    if len(vectors) == 1:
                        return [residual(vectors[0])]
                    return pool.evaluate_many([expand(q) for q in vectors])
            jac = _block_jacobian(
                residual,
                sizes,
                self.n_global,
                nlocal,
                relative_step=1e-5 if not full_order or self.proton_mode == "fit" else 1e-6,
                free=self.free,
                bounds=(lower[self.free], upper[self.free]),
                method="2-point" if full_order else "3-point",
                evaluate_many=evaluate_many,
            )
            self.result = least_squares(
                residual,
                p0[self.free],
                jac=jac,
                diff_step=1e-5,
                bounds=(lower[self.free], upper[self.free]),
                x_scale="jac",
                ftol=1e-9,
                xtol=1e-9,
                gtol=1e-9,
                max_nfev=cfg.get("max_nfev", 500),
                verbose=2 if self.verbose else 0,
            )
            if not self.result.success:
                raise RuntimeError(
                    f"Sideband fit did not converge: {self.result.message}"
                )
            p = expand(self.result.x)
            self.errFunc(p)
            self.npar = len(self.free)
            self.dof = self.nvar - self.npar
            # Scale columns for an interpretable rank test across Hz/ppm/rate units.
            j = self.result.jac
            norms = np.linalg.norm(j, axis=0)
            scaled = j / np.where(norms > 0, norms, 1.0)
            singular = np.linalg.svd(scaled, compute_uv=False)
            self.rank = int(np.count_nonzero(singular > singular[0] * 1e-8))
            self.condition = (
                float(singular[0] / singular[-1]) if singular[-1] else np.inf
            )
            covariance = np.zeros((len(p), len(p)))
            cfree = np.full((len(self.free), len(self.free)), np.nan)
            if self.rank == len(self.free):
                _, s, vt = np.linalg.svd(j, full_matrices=False)
                cfree = (vt.T / s**2) @ vt
            covariance[np.ix_(self.free, self.free)] = cfree
            return p, covariance
        finally:
            self._fit_data = None

    def _log_stats(self, p_ref, mc=False):
        lines = super()._log_stats(p_ref, mc)
        self.npar = len(self.free)
        self.dof = self.nvar - self.npar
        return [
            line.replace("Matrix", "Sideband")
            if "Results using" in line
            else f"np:       {self.npar:d}"
            if line.startswith("np:")
            else f"dof:      {self.dof:d}"
            if line.startswith("dof:")
            else f"red_chi2: {self.chi2 / self.dof:.6f}"
            if line.startswith("red_chi2:")
            else line
            for line in lines
        ]

    def _log_params(self, values, stds):
        lines = self._parameter_lines(values, stds)
        try:
            from sb_diagnostics import diagnostics_lines, residual_diagnostics

            report = residual_diagnostics(prediction_rows(self, values), len(self.free),
                                          self._last_covariance)
            lines.extend(["*" * 43, *diagnostics_lines(report, self.parameter_names)])
        except (ValueError, AttributeError):
            pass
        return lines

    def getLogBuffer(self, fit_output_tuple):
        """Text report of a fit, including residual diagnostics and rescaled errors."""
        # Remember the covariance so the text report can list rescaled errors.
        self._last_covariance = np.asarray(fit_output_tuple[1], dtype=float)
        try:
            return super().getLogBuffer(fit_output_tuple)
        finally:
            self._last_covariance = None

    def _parameter_lines(self, values, stds):
        n_rates = len(self.rate_names)
        lines = [
            f"{n}: {values[n_rates + i]:.8g} +/- {stds[n_rates + i]:.6g}"
            for i, n in enumerate(self.rf_names)
        ]
        lines += super()._log_params(
            self._nitrogen_parameters(values), self._nitrogen_parameters(stds)
        )
        if self.nitrogen_groups > 1:
            lines.append("Nitrogen relaxation mode: per_field (the block above lists field group 0)")
            lines.extend(
                f"{name} [s-1]: {values[i]:.8g} +/- {stds[i]:.6g}"
                for i, name in enumerate(self.parameter_names)
                if _key_base(name) in self.n_relax_keys and "[" in name
            )
        lines.append(f"Proton relaxation mode: {self.proton_mode}")
        if self.proton_mode == "fit":
            lines.extend(
                f"{name} [s-1]: {values[i]:.8g} +/- {stds[i]:.6g}"
                for i, name in enumerate(self.parameter_names)
                if _key_base(name) in self.h_keys
            )
        if len(self.field_groups) > 1:
            lines.extend(
                f"field group {g}: {larmor:g} MHz (datasets "
                + ", ".join(str(i) for i, gi in enumerate(self.dataset_group) if gi == g) + ")"
                for g, larmor in enumerate(self.field_groups)
            )
        P = self.seParam(values)
        for i, nominal in enumerate(self.dataset.v1s):
            value = (
                nominal * P["v1n"][0]
                if self.rf_mode == "scale"
                else (P["v1n"][i] if self.rf_mode == "per_dataset" else nominal)
            )
            lines.append(f"dataset {i}: v1n = {value:.8g} Hz (nominal {nominal:g} Hz)")
        lines.append(
            "Fixed parameters: "
            + ", ".join(
                name
                for i, name in enumerate(self.parameter_names)
                if i not in self.free
            )
        )
        return lines

    def diagnostics(self, p, covariance):
        """Result dictionary (schema 2) for a fitted vector and its full covariance."""
        std = np.sqrt(np.diag(covariance))
        n_rates = len(self.rate_names)
        pairs = []
        for i in range(n_rates, self.n_global):
            if std[i] > 0:
                for k in self.free:
                    if i != k and std[k] > 0:
                        pairs.append(
                            {
                                "parameters": [
                                    self.parameter_names[i],
                                    self.parameter_names[k],
                                ],
                                "correlation": float(
                                    covariance[i, k] / (std[i] * std[k])
                                ),
                            }
                        )
        pairs.sort(key=lambda row: -abs(row["correlation"]))
        bound = [
            self.parameter_names[self.free[k]]
            for k in np.flatnonzero(self.result.active_mask)
        ]
        warnings = []
        proton_pairs = []
        weak_proton = []
        for i, name in enumerate(self.parameter_names):
            if _key_base(name) in self.h_keys and i in self.free:
                if not np.isfinite(std[i]) or std[i] >= p[i]:
                    weak_proton.append(name)
            if _key_base(name) == "R1H" and std[i] > 0 and std[i + 1] > 0:
                proton_pairs.append(
                    {
                        "parameters": [name, self.parameter_names[i + 1]],
                        "correlation": float(
                            covariance[i, i + 1] / (std[i] * std[i + 1])
                        ),
                    }
                )
        if weak_proton:
            warnings.append(
                "Weakly constrained proton nuisance parameters: "
                + ", ".join(weak_proton)
            )
        if any(abs(row["correlation"]) > 0.95 for row in proton_pairs):
            warnings.append(
                "Strong peakwise R1H/R2H correlation; individual H rates are poorly separated"
            )
        if self.rank < len(self.free):
            warnings.append(
                "Rank-deficient Jacobian: parameters are not separately identifiable"
            )
        if bound:
            warnings.append("Parameters at bounds: " + ", ".join(bound))
        if pairs and abs(pairs[0]["correlation"]) > 0.95:
            warnings.append(
                "Strong v1n correlation; inspect RF calibration and profile likelihood"
            )
        if self.condition > 1e6:
            warnings.append(
                "Ill-conditioned scaled Jacobian; local covariance may be unstable"
            )
        two_state = self.states == 2
        P = self.seParam(p)
        populations = stationary_populations(self.exchange_matrix(P))
        info = {
            "schema_version": 2,
            "method": "Sideband",
            "model": self.sideband_method,
            "states": self.states,
            "rf_mode": self.rf_mode,
            "proton_relaxation_mode": self.proton_mode,
            "nitrogen_relaxation_mode": self.nitrogen_mode,
            "field_groups": [{"index": g, "h_larmor_mhz": larmor,
                              "datasets": [i for i, gi in enumerate(self.dataset_group) if gi == g]}
                             for g, larmor in enumerate(self.field_groups)],
            "proton_correlations": proton_pairs,
            "config": self.config,
            "v1n_hz": [
                float(
                    nominal * p[n_rates]
                    if self.rf_mode == "scale"
                    else p[n_rates + i]
                    if self.rf_mode == "per_dataset"
                    else nominal
                )
                for i, nominal in enumerate(self.dataset.v1s)
            ],
            "chi2": float(self.chi2),
            "dof": self.nvar - len(self.free),
            "n_points": self.nvar,
            "n_parameters": len(self.free),
            "kex": float(p[0] + p[1]) if two_state else None,
            "pB": float(p[0] / (p[0] + p[1])) if two_state else None,
            "exchange": {
                "rates": {name: float(p[i]) for i, name in enumerate(self.rate_names)},
                "populations": {site: float(value) for site, value in zip("ABC", populations)},
                "kex_AB": float(p[0] + p[1]),
                **({"kex_BC": float(P["kbc"] + P["kcb"])} if not two_state else {}),
                "note": ("kex and pB summarize the two-state model; for three sites the "
                         "populations follow from all rates and kex/pB are null."),
            },
            "parameter_order": self.parameter_names,
            "covariance": [[float(v) if np.isfinite(v) else None for v in row]
                           for row in covariance],
            "derived_se": derived_errors(p, covariance) if two_state else {"kex": None, "pB": None},
            "parameters": {
                n: {
                    "value": float(p[i]),
                    "stderr": float(std[i]) if np.isfinite(std[i]) else None,
                    "vary": i in self.free,
                }
                for i, n in enumerate(self.parameter_names)
            },
            "jacobian_rank": self.rank,
            "scaled_condition": self.condition if np.isfinite(self.condition) else None,
            "v1n_correlations": pairs,
            "at_bounds": bound,
            "warnings": warnings,
            "covariance_note": (
                "Local linear covariance using supplied absolute intensity errors; "
                "not rescaled by chi2/dof"
            ),
            "nfev": self.result.nfev,
            "message": self.result.message,
        }
        return info


# check_config, run_config and the command line live in sb_run.py; they stay
# importable from this module for existing callers and scripts.
def __getattr__(name):
    if name in ("check_config", "run_config", "_output_paths", "main"):
        import sb_run

        return getattr(sb_run, name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


if __name__ == "__main__":
    from sb_run import main

    main()
