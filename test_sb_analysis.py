"""Exact quadratic checks and a small real Sideband analysis integration test."""

# ruff: noqa: E402 -- Limit numerical libraries before importing them.
import os

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(name, "1")

import copy
import json
from pathlib import Path
import tempfile

import numpy as np
from numpy.testing import assert_allclose
from scipy.optimize import least_squares

from sb_analysis import MultiStartError, fit_multistart, profile_likelihood


class QuadraticModel:
    """An independently solvable objective with deliberately reordered names."""

    parameter_names = ["nuisance", "kba", "v1n_scale", "kab"]

    def __init__(self, vary=None, residual=None):
        self.config = {"init": {"initial": {}, "max_nfev": 200}}
        self.free = np.array([3, 0, 1, 2] if vary is None else vary, dtype=int)
        self.config["init"]["vary"] = [self.parameter_names[i] for i in self.free]
        self.lower = np.array([-np.inf, 1e-8, 0.5, 0.0])
        self.upper = np.array([np.inf, np.inf, 1.5, np.inf])
        self.initial_parameters = np.array([15., 8., 1., 2.])
        self.result = None
        self.residual = residual
        self.rank = len(self.free)
        self.condition = 1.0
        self.npar = len(self.free)
        self.nvar = 5
        self.dof = self.nvar - self.npar

    def errFunc(self, p):
        n, b, rf, a = p
        r = np.asarray(self.residual(p) if self.residual else
                       [a - 2., b - 8., n - a - 2*b, rf - 1., 0.2])
        self.chi2 = float(r @ r)
        self.nvar, self.npar = len(r), len(p)
        self.dof = self.nvar - self.npar
        return r

    def fit(self, p0=None, fitting_config=None):
        cfg = self.config["init"] if fitting_config is None else fitting_config
        start = self.initial_parameters.copy() if p0 is None else np.array(p0).copy()
        for name, value in cfg.get("initial", {}).items():
            start[self.parameter_names.index(name)] = value
        self.initial_parameters = start.copy()

        def expand(x):
            p = start.copy()
            p[self.free] = x
            return p

        self.result = least_squares(lambda x: self.errFunc(expand(x)), start[self.free],
                                    bounds=(self.lower[self.free], self.upper[self.free]),
                                    max_nfev=cfg.get("max_nfev", 200),
                                    ftol=1e-12, xtol=1e-12, gtol=1e-12)
        if not self.result.success:
            raise RuntimeError(self.result.message)
        p = expand(self.result.x)
        self.errFunc(p)
        self.npar, self.dof = len(self.free), self.nvar - len(self.free)
        return p, np.eye(len(p))


def check_quadratic_profiles():
    # Failure caught: a scan that fixes nuisance parameters overstates curvature.
    model = QuadraticModel()
    p, _ = model.fit()
    state = model.result, model.free.copy(), model.chi2, model.npar, model.dof
    out = profile_likelihood(model, p, {"kex": [10., 12.], "pB": [0.25], "v1n_scale": [1.2]})
    assert_allclose(out["base_chi2"], 0.04, atol=1e-12)
    assert_allclose(out["kex"][1]["parameters"], [21., 9., 1., 3.], atol=1e-6)
    assert_allclose(out["kex"][1]["delta_chi2"], 2., atol=1e-8)
    assert_allclose(out["pB"][0]["parameters"], [18.2, 7.8, 1., 2.6], atol=1e-6)
    assert_allclose(out["pB"][0]["delta_chi2"], 0.4, atol=1e-8)
    assert_allclose(out["v1n_scale"][0]["delta_chi2"], 0.04, atol=1e-8)
    assert model.result is state[0]
    assert_allclose(model.free, state[1])
    assert_allclose([model.chi2, model.npar, model.dof], state[2:])


def check_fixed_rates_and_bounds():
    # Failure caught: reconstructed fixed rates drifting or user bounds being ignored.
    for vary, setting, wanted in [
        ([3, 0], {"kex": [12.]}, [20., 8., 1., 4.]),
        ([3, 0], {"pB": [0.25]}, [56./3, 8., 1., 8./3]),
        ([1, 0], {"pB": [0.25]}, [14., 6., 1., 2.]),
    ]:
        model = QuadraticModel(vary)
        p, _ = model.fit()
        key = next(iter(setting))
        row = profile_likelihood(model, p, setting)[key][0]
        assert row["success"], row
        assert_allclose(row["parameters"], wanted, atol=1e-6)
        for i in set(range(4)) - set(vary):
            assert row["parameters"][i] == p[i]

    model = QuadraticModel([3])
    model.upper[3] = 4.
    p, _ = model.fit()
    out = profile_likelihood(model, p, {"kex": [12., 13.], "pB": [0., 1.]})
    assert out["kex"][0]["success"]
    assert not out["kex"][0]["optimization_performed"]
    assert out["kex"][0]["parameters"][3] == 4.
    assert not out["kex"][1]["success"] and out["kex"][1]["chi2"] is None
    assert out["warnings"], "Failed scan targets must be visible in warnings"
    assert out["pB"][0]["success"] and not out["pB"][1]["success"]

    model = QuadraticModel([0, 2])
    p, _ = model.fit()
    out = profile_likelihood(model, p, {"kex": [10., 11.], "pB": [0.2, 0.3], "v1n_scale": [0.2]})
    assert out["kex"][0]["success"] and not out["kex"][1]["success"]
    assert out["pB"][0]["success"] and not out["pB"][1]["success"]
    assert not out["v1n_scale"][0]["success"]

    model = QuadraticModel()
    model.lower[[3, 1]], model.upper[[3, 1]] = [2., 3.], [8., 12.]
    p, _ = model.fit()
    out = profile_likelihood(model, p, {"kex": [4., 20.], "pB": [0.2]})
    assert not out["kex"][0]["success"]
    assert_allclose(np.array(out["kex"][1]["parameters"])[[3, 1]], [8., 12.])
    assert out["kex"][1]["success"]
    assert out["pB"][0]["success"]
    assert np.array(out["pB"][0]["parameters"])[3] >= 2.




def check_population_endpoints():
    # Failure caught: dividing by zero or dropping an unconstrained endpoint rate.
    for vary, rates, target in [([3, 0], [2., 0.], 1.), ([1, 0], [0., 8.], 0.)]:
        model = QuadraticModel(vary)
        model.lower[1] = 0.
        model.initial_parameters[[3, 1]] = rates
        p, _ = model.fit()
        out = profile_likelihood(model, p, {"pB": [target, .5]})
        assert out["pB"][0]["success"] and not out["pB"][1]["success"]
        assert out["pB"][0]["optimization_performed"]
        assert len(out["pB"][0]["nuisance_parameters"]) == len(vary)
    model = QuadraticModel()
    model.lower[1] = 0.
    p, _ = model.fit()
    out = profile_likelihood(model, p, {"pB": [0., 1.]})
    assert all(row["success"] for row in out["pB"])
    assert out["pB"][0]["parameters"][3] == 0.
    assert out["pB"][1]["parameters"][1] == 0.


def check_profile_roundoff_boundaries():
    # Failure caught: one-ULP reconstruction errors being called infeasible.
    model = QuadraticModel()
    model.lower[[3, 1]], model.upper[[3, 1]] = [.1, .2], [1., 1.]
    model.initial_parameters[[3, 1]] = [.5, .5]
    p, _ = model.fit()
    row = profile_likelihood(model, p, {"kex": [.3]})["kex"][0]
    assert row["success"], row
    assert_allclose(np.array(row["parameters"])[[3, 1]], [.1, .2], atol=1e-15)

    model = QuadraticModel()
    model.lower[3] = 7.221648081421175
    model.initial_parameters[[3, 1]] = [model.lower[3], .001]
    p = model.initial_parameters.copy()
    row = profile_likelihood(model, p, {"pB": [.22434111607742846]})["pB"][0]
    assert row["success"], row
    rates = np.array(row["parameters"])[[3, 1]]
    assert rates[0] >= model.lower[3] - 1e-13
    assert_allclose(rates[0] / rates.sum(), .22434111607742846, atol=1e-15)
    model = QuadraticModel([3, 1])
    p = model.initial_parameters.copy()
    p[3] = 150.
    row = profile_likelihood(model, p, {"kex": [100.]})["kex"][0]
    assert row["success"], row
    assert np.array(row["start"])[1] >= model.lower[1]
    assert_allclose(np.array(row["start"])[[3, 1]].sum(), 100., atol=1e-13)



def check_infeasible_small_rates():
    # Failure caught: an absolute tolerance floor moving a small-rate pB target.
    model = QuadraticModel([3, 1])
    model.upper[3] = 1e-8
    p = model.initial_parameters.copy()
    p[[3, 1]] = [5e-9, 1e-8]
    row = profile_likelihood(model, p, {"pB": [.50000001]})["pB"][0]
    assert not row["success"] and row["chi2"] is None, row
    model = QuadraticModel([0])
    p = model.initial_parameters.copy()
    row = profile_likelihood(model, p, {"pB": [.2 + 1e-13]})["pB"][0]
    assert not row["success"], "Fixed rates cannot satisfy a different population"


def check_multistart_selection_and_failures():
    # Failure caught: last fit overwriting a better minimum or discarded failures.
    def residual(p):
        x = p[0]
        if x > 9:
            raise ValueError("deliberate bad domain")
        return [x*x-4., 0.25*(x-2.), 0.2]

    model = QuadraticModel([0], residual)
    model.initial_parameters[0] = -3.
    model.config["init"]["initial"] = {"nuisance": -3.}
    p, _, attempts = fit_multistart(model, {"starts": [{"nuisance": 3.}, {"nuisance": -3.}, {"nuisance": 10.}]})
    assert_allclose(p[0], 2., atol=1e-6)
    assert len(attempts) == 4 and [r["success"] for r in attempts] == [True, True, True, False]
    assert [r["selected"] for r in attempts] == [False, True, False, False]
    assert attempts[1]["start"][0] == 3.
    assert attempts[-1]["start"][0] == 10.
    assert_allclose(model.result.x, [2.], atol=1e-6)
    assert_allclose([model.chi2, model.npar, model.dof], [0.04, 1, 2], atol=1e-8)

    model = QuadraticModel([0], residual)
    model.initial_parameters[0] = 10.
    try:
        fit_multistart(model, {"starts": [{"nuisance": 11.}]})
    except MultiStartError as exc:
        assert len(exc.attempts) == 2
        assert all(not row["success"] for row in exc.attempts)
    else:
        raise AssertionError("All failed starts must raise with retained evidence")


def check_random_starts_and_validation():
    # Failure caught: nondeterministic seeds or a random start altering fixed rates.
    records = []
    for _ in range(2):
        model = QuadraticModel([0, 2])
        _, _, rows = fit_multistart(model, {"seed": 317, "random_starts": 3})
        records.append([r["start"] for r in rows])
        assert all(r["start"][1] == 8. and r["start"][3] == 2. for r in rows)
        assert all(0.5 <= r["start"][2] <= 1.5 for r in rows)
    assert_allclose(records[0], records[1], rtol=0, atol=0)
    for settings in ({"random_starts": 1}, {"random_starts": -1}, {"seed": True},
                     {"starts": [{"kab": 5.}]}, {"starts": [{"missing": 1.}]}):
        model = QuadraticModel([0, 2])
        try:
            fit_multistart(model, settings)
        except ValueError:
            assert model.result is None, "Bad settings must fail before fitting"
        else:
            raise AssertionError(f"Invalid settings accepted: {settings}")


def check_lower_profile_and_failed_solver():
    # Failure caught: negative delta-chi2 hidden or nonconvergence called success.
    model = QuadraticModel()
    p, _ = model.fit()
    bad_base = p.copy()
    bad_base[0] += 3.
    out = profile_likelihood(model, bad_base, {"kex": [10.]})
    assert out["kex"][0]["delta_chi2"] < -8.9
    assert out["kex"][0]["below_base_minimum"] and out["warnings"]
    model.config["init"]["max_nfev"] = 1
    out = profile_likelihood(model, bad_base, {"kex": [12.]})
    assert not out["kex"][0]["success"] and out["kex"][0]["chi2"] is None
    assert out["kex"][0]["optimization_performed"]


def check_real_sideband():
    # Small integration: the real model's initial overrides and stats must survive.
    from sbfit import SidebandModel
    from sideband import composite_segments, profile

    folder = Path(tempfile.mkdtemp(prefix="sbonest-analysis-"))
    path = folder / "synthetic.txt"
    offsets = np.array([-120., -60., -10., 0., 30., 80., 160.])
    y = profile(composite_segments(70e-6), offsets, T=0.02, nu=35., kab=15., kba=285.,
                dw=180., r1=1.5, r2a=12., r2b=18., ha=-1380., hb=-1200., J=92., r1h=2., r2h=25.)
    path.write_text("60.0\n0.02\n35.0 0\n# offset intensity error\n# A1 R2a: 12 R2b: 18 dw: 3\n" +
                    "".join(f"{120+x/60:.15g} {v:.15g} 0.001\n" for x, v in zip(offsets, y)))
    cfg = {"Project Name": str(folder / "fit"), "datasets": [str(path)],
           "residues": [{"name": "A1", "flag": "on"}],
           "init": {"Method": "Sideband", "kex": {"min": 300, "max": 300, "nsteps": 1},
                    "pB": {"min": .05, "max": .05, "nsteps": 1}, "vary": ["v1n_scale"],
                    "initial": {"A1.peak_ppm": 120., "A1.dw_ppm": 3., "A1.R1": 1.5,
                                "A1.R2a": 12., "A1.R2b": 18., "v1n_scale": .95}},
           "sideband": {"decoupling": {"h_larmor_mhz": 600., "h_carrier_ppm": 8.5,
                                         "p90_s": 70e-6, "R1H": 2., "R2H": 25.},
                        "residues": {"A1": {"h_ppm_a": 6.2, "h_ppm_b": 6.5}},
                        "v1n": {"mode": "scale", "initial": .95, "bounds": [.8, 1.2]}}}
    model = SidebandModel(cfg)
    p, covariance, attempts = fit_multistart(model, {"starts": [{"v1n_scale": 1.1}]})
    assert_allclose(p[2], 1., atol=1e-6)
    assert attempts[1]["start"][2] == 1.1
    assert_allclose(p[[0, 1, 3, 4, 5, 6, 7]], [15., 285., 120., 3., 1.5, 12., 18.])
    stats = model.result, model.chi2, model.npar, model.dof
    out = profile_likelihood(model, p, {"kex": [300., 320.], "v1n_scale": [1., 1.05]})
    assert out["kex"][0]["success"] and not out["kex"][1]["success"]
    assert out["v1n_scale"][0]["chi2"] < 1e-14
    assert out["v1n_scale"][1]["chi2"] > 1.
    assert model.result is stats[0] and covariance.shape == (8, 8)
    assert_allclose([model.chi2, model.npar, model.dof], stats[1:])
    json.dumps({"attempts": attempts, "profiles": out}, allow_nan=False)

    # Failure caught: rejected baseline preventing valid restarts, or stale
    # initialization on a reused model silently restoring an old fixed R1.
    for reused in (False, True):
        candidate = model if reused else SidebandModel(copy.deepcopy(cfg))
        candidate.config["init"]["bounds"] = {"v1n_scale": [1.05, 1.15]}
        candidate.config["init"]["initial"]["A1.R1"] = 3.
        fitted, _, rows = fit_multistart(candidate, {"starts": [{"v1n_scale": 1.10}],
                                                    "seed": 9, "random_starts": 1})
        assert len(rows) == 3 and not rows[0]["success"]
        assert rows[0]["start"][2] == .95 and rows[0]["start"][5] == 3.
        assert rows[1]["success"] and rows[2]["success"]
        assert fitted[5] == 3. and 1.05 <= fitted[2] <= 1.15
        json.dumps(rows, allow_nan=False)
    print(f"PASS: real Sideband analysis; synthetic inputs retained at {folder}")


if __name__ == "__main__":
    check_quadratic_profiles()
    check_fixed_rates_and_bounds()
    check_population_endpoints()
    check_profile_roundoff_boundaries()
    check_infeasible_small_rates()
    check_multistart_selection_and_failures()
    check_random_starts_and_validation()
    check_lower_profile_and_failed_solver()
    check_real_sideband()
    print("PASS: exact nuisance profiles, fixed rates/bounds, starts/failures/state")
