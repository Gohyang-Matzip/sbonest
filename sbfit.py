#!/usr/bin/env python3
"""SBONEST: fit all CEST points, including decoupling sidebands, with OC.

Usage: python sbfit.py config.json [--no-pdf]
ONEST text inputs and reports are retained; see SIDEBAND.md for configuration.
"""

import argparse
import json
from pathlib import Path
import time

import numpy as np
from scipy.optimize import least_squares

from estmodel import b1_weights, est_model
from fit import _block_jacobian, generate_initial_parameters
from run import load_config, set_residue_flags
from sideband import composite_segments, profile, waveform_segments
from sb_report import derived_errors, fit_pdf, prediction_rows, provenance, write_predictions
from sb_analysis import MultiStartError, fit_multistart, profile_likelihood


def _known(config, keys, label):
    if not isinstance(config, dict) or set(config) - set(keys):
        raise ValueError(f"Invalid {label} keys; allowed: {', '.join(keys)}")


class SidebandModel(est_model):
    """Reuse ONEST data, residue ordering, residuals and PDF/text reporting."""

    def __init__(self, config, config_dir="."):
        super().__init__()
        self.method = "Matrix"  # ONEST's two-state parameter/report layout
        self.programName = "SBONEST 1.0 — exact NH Sideband model (OC)"
        self.config = config
        self._initializing = False
        self.result = None
        sb = config["sideband"]
        _known(
            sb,
            ("decoupling", "datasets", "residues", "v1n", "proton_relaxation"),
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
        for r in active:
            h = self.h_shifts[r.label]
            _known(h, ("h_ppm_a", "h_ppm_b"), f"sideband.residues.{r.label}")
            if not np.isfinite([h["h_ppm_a"], h["h_ppm_b"]]).all():
                raise ValueError(
                    "Supply finite h_ppm_a and h_ppm_b for every active residue"
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
        self.local_keys = ("peak_ppm", "dw_ppm", "R1", "R2a", "R2b")
        if self.proton_mode == "fit":
            self.local_keys += ("R1H", "R2H")
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
        if (
            self.proton_mode == "fit"
            and len({d["h_larmor_mhz"] for d in self.decoupling}) > 1
        ):
            raise ValueError(
                "Automatic peakwise proton relaxation requires one proton field; "
                "fit fields separately or supply fixed rates for each dataset"
            )
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
        self.parameter_names = (
            ["kab", "kba"]
            + self.rf_names
            + [f"{r.label}.{key}" for r in active for key in self.local_keys]
        )
        self.free = np.arange(len(self.parameter_names))

    def selMethod(self, initConf):
        if initConf.get("Method", "Sideband") != "Sideband":
            raise ValueError("SBONEST supports the two-state Sideband model")

    def _expand_initial(self, base):
        local = np.asarray(base[2:]).reshape(-1, 5)
        if self.proton_mode == "fit":
            local = np.column_stack((local, np.tile([2.0, 25.0], (len(local), 1))))
        return np.r_[base[:2], self.rf_initial, local.ravel()]

    def _nitrogen_parameters(self, p):
        local = np.asarray(p[2 + len(self.rf_names) :]).reshape(
            -1, len(self.local_keys)
        )
        return np.r_[p[:2], local[:, :5].ravel()]

    def seParam(self, p):
        p = np.asarray(p, dtype=float)
        nr = len(self.rf_names)
        if self._initializing:
            p = self._expand_initial(p)
        if p.shape != (len(self.parameter_names),):
            raise ValueError("Wrong Sideband parameter vector length")
        P = super().seParam(self._nitrogen_parameters(p))
        P["v1n"] = p[2 : 2 + nr]
        if self.proton_mode == "fit":
            local = p[2 + nr :].reshape(-1, len(self.local_keys))
            P["proton_rates"] = {
                r.label: row[5:7]
                for r, row in zip((r for r in self.dataset.res if r.active), local)
            }
        return P

    def rf_values(self, P, es):
        if self.rf_mode == "scale":
            scale = P["v1n"][0]
        elif self.rf_mode == "per_dataset":
            scale = P["v1n"][es.dataset_index] / es.v1
        else:
            scale = 1.0
        # Keep fractional RF inhomogeneity fixed, independently of fit uncertainty.
        return es.v1 * scale, es.v1err * scale

    def calc(self, P, i, dRF, es):
        d = self.decoupling[es.dataset_index]
        h = self.h_shifts[self.dataset.res[i].label]
        r1h, r2h = (
            P["proton_rates"][self.dataset.res[i].label]
            if self.proton_mode == "fit"
            else (d["R1H"], d["R2H"])
        )
        fields, weights = b1_weights(*self.rf_values(P, es))
        out = np.zeros_like(np.asarray(dRF, dtype=float))
        for nu, weight in zip(fields / (2 * np.pi), weights):
            out += weight * profile(
                self.segments[es.dataset_index],
                (np.asarray(dRF) - P["dGs"][i]) * es.field,
                T=es.T,
                nu=abs(nu),
                kab=P["kab"],
                kba=P["kba"],
                dw=P["dws"][i] * es.field,
                r1=P["r1s"][i],
                r2a=P["r2as"][i],
                r2b=P["r2bs"][i],
                ha=(h["h_ppm_a"] - d["h_carrier_ppm"]) * d["h_larmor_mhz"],
                hb=(h["h_ppm_b"] - d["h_carrier_ppm"]) * d["h_larmor_mhz"],
                J=d["J_hz"],
                r1h=r1h,
                r2h=r2h,
            )
        return out if np.ndim(dRF) else float(out)

    def fit(self, p0=None, fitting_config=None):
        cfg = fitting_config or self.config["init"]
        _known(
            cfg,
            ("Method", "kex", "pB", "initial", "bounds", "vary", "max_nfev",
             "multistart", "profile"),
            "init",
        )
        self.selMethod(cfg)
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
            nr = len(self.rf_names)
            nlocal = len(self.local_keys)
            lower = np.r_[
                [0.0, 1e-8],
                self.rf_bounds[:, 0],
                np.tile(
                    [-np.inf, -np.inf] + [0.0] * (nlocal - 2),
                    (len(p0) - 2 - nr) // nlocal,
                ),
            ]
            upper = np.r_[
                [np.inf, np.inf],
                self.rf_bounds[:, 1],
                np.full(len(p0) - 2 - nr, np.inf),
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
                    if 2 + i in self.free and i not in used:
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

            def expand(q):
                p = p0.copy()
                p[self.free] = q
                return p

            def residual(q):
                values = self.errFunc(expand(q))
                if not np.isfinite(values).all():
                    raise ValueError("Non-finite Sideband residual")
                return values

            sizes = [
                sum(len(es.offset) for es in r.estSpecs)
                for r in self.dataset.res
                if r.active
            ]
            full_order = np.array_equal(self.free, np.arange(len(p0)))
            jac = _block_jacobian(
                residual,
                sizes,
                2 + nr,
                nlocal,
                relative_step=1e-5 if not full_order or self.proton_mode == "fit" else 1e-6,
                free=self.free,
                bounds=(lower[self.free], upper[self.free]),
                method="2-point" if full_order else "3-point",
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
        lines = [
            f"{n}: {values[2 + i]:.8g} +/- {stds[2 + i]:.6g}"
            for i, n in enumerate(self.rf_names)
        ]
        lines += super()._log_params(
            self._nitrogen_parameters(values), self._nitrogen_parameters(stds)
        )
        lines.append(f"Proton relaxation mode: {self.proton_mode}")
        if self.proton_mode == "fit":
            lines.extend(
                f"{name} [s-1]: {values[i]:.8g} +/- {stds[i]:.6g}"
                for i, name in enumerate(self.parameter_names)
                if name.endswith((".R1H", ".R2H"))
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
        std = np.sqrt(np.diag(covariance))
        pairs = []
        for i in range(2, 2 + len(self.rf_names)):
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
            if name.endswith((".R1H", ".R2H")) and i in self.free:
                if not np.isfinite(std[i]) or std[i] >= p[i]:
                    weak_proton.append(name)
            if name.endswith(".R1H") and std[i] > 0 and std[i + 1] > 0:
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
        return {
            "schema_version": 2,
            "method": "Sideband",
            "rf_mode": self.rf_mode,
            "proton_relaxation_mode": self.proton_mode,
            "proton_correlations": proton_pairs,
            "config": self.config,
            "v1n_hz": [
                float(
                    nominal * p[2]
                    if self.rf_mode == "scale"
                    else p[2 + i]
                    if self.rf_mode == "per_dataset"
                    else nominal
                )
                for i, nominal in enumerate(self.dataset.v1s)
            ],
            "chi2": float(self.chi2),
            "dof": self.nvar - len(self.free),
            "n_points": self.nvar,
            "n_parameters": len(self.free),
            "kex": float(p[0] + p[1]),
            "pB": float(p[0] / (p[0] + p[1])),
            "parameter_order": self.parameter_names,
            "covariance": [[float(v) if np.isfinite(v) else None for v in row]
                           for row in covariance],
            "derived_se": derived_errors(p, covariance),
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


def run_config(config, config_dir=".", no_pdf=False):
    model = SidebandModel(config, config_dir)
    project = Path(config["Project Name"])
    outputs = [
        Path(str(project) + suffix)
        for suffix in ("_result.txt", "_result.json", "_predictions.csv")
    ]
    if not no_pdf:
        outputs += [Path(str(project) + ".pdf"), Path(str(project) + "_data.pdf")]
    if any(path.exists() or path.is_symlink() for path in outputs):
        raise FileExistsError(
            "Output already exists; choose a new Project Name to preserve previous fits"
        )
    project.parent.mkdir(parents=True, exist_ok=True)
    metadata = provenance(config, config_dir)
    began = time.monotonic()

    def save_info(info):
        metadata["elapsed_s"] = time.monotonic() - began
        info["provenance"] = metadata
        with outputs[1].open("x", encoding="utf-8") as stream:
            stream.write(json.dumps(info, indent=2, allow_nan=False) + "\n")

    model.verbose = True
    attempts = None
    if "multistart" in config["init"]:
        try:
            p, covariance, attempts = fit_multistart(model, config["init"]["multistart"])
            fitted = p, covariance
        except MultiStartError as exc:
            save_info({"schema_version": 2, "success": False, "method": "Sideband",
                       "config": config, "parameter_order": model.parameter_names,
                       "multistart": exc.attempts, "message": str(exc)})
            raise
    else:
        fitted = model.fit(fitting_config=config["init"])
    info = model.diagnostics(*fitted)
    info["success"] = True
    if attempts is not None:
        info["multistart"] = attempts
    if "profile" in config["init"]:
        info["profiles"] = profile_likelihood(model, fitted[0], config["init"]["profile"])
        info["warnings"].extend(info["profiles"]["warnings"])
    rows = prediction_rows(model, fitted[0])
    with outputs[0].open("x", encoding="utf-8") as stream:
        stream.write(model.getLogBuffer(fitted))
    save_info(info)
    write_predictions(outputs[2], rows)
    if not no_pdf:
        fit_pdf(str(project) + ".pdf", rows)
        model.datapdf(str(project) + "_data.pdf")
    for warning in info["warnings"]:
        print(f"Warning: {warning}")
    print(f"Saved {outputs[0]} and {outputs[1]}")
    return fitted


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config_file")
    parser.add_argument("--no-pdf", action="store_true")
    args = parser.parse_args()
    try:
        run_config(
            load_config(args.config_file),
            Path(args.config_file).resolve().parent,
            args.no_pdf,
        )
    except (ValueError, KeyError, OSError, RuntimeError) as exc:
        parser.exit(1, f"Error: {exc}\n")


if __name__ == "__main__":
    main()
