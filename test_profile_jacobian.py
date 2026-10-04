"""Independent dense checks of constrained, residue-grouped profile derivatives."""

# ruff: noqa: E402 -- Limit numerical libraries before importing them.
import os

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(name, "1")

from pathlib import Path
import tempfile
import time
from unittest.mock import patch

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from scipy.optimize._numdiff import approx_derivative

from sb_analysis import _constrained_coordinates, _profile_jacobian, profile_likelihood
from sbfit import SidebandModel
from sideband import composite_segments, profile
from test_sb_analysis import QuadraticModel


def real_model():
    """Three actual NH propagator blocks, differing point counts, two RF values."""
    folder = Path(tempfile.mkdtemp(prefix="sbonest-profile-jacobian-"))
    residues = [("A1", 120., 3., 6.2, 6.5, 9),
                ("G2", 110., -2., 7.2, 7.5, 10),
                ("S3", 115., 4., 9.8, 9.6, 11)]
    paths, initial = [], {"kab": 15., "kba": 285., "v1n_scale": 1.}
    for nu in (25., 80.):
        path = folder / f"rf{int(nu)}.txt"
        lines = ["60.0", "0.02", f"{nu} 0", "# offset intensity error"]
        for label, peak, dw, ha, hb, size in residues:
            offsets = np.linspace(-240., 240., size)
            y = profile(composite_segments(70e-6), offsets, T=.02, nu=nu,
                        kab=15., kba=285., dw=dw*60., r1=1.5, r2a=12., r2b=18.,
                        ha=(ha-8.5)*600., hb=(hb-8.5)*600., J=92., r1h=2., r2h=25.)
            lines.append(f"# {label} R2a: 12 R2b: 18 dw: {dw}")
            lines.extend(f"{peak+x/60:.15g} {value:.15g} 0.005" for x, value in zip(offsets, y))
            initial.update({f"{label}.{key}": value for key, value in zip(
                ("peak_ppm", "dw_ppm", "R1", "R2a", "R2b", "R1H", "R2H"),
                (peak, dw, 1.5, 12., 18., 2., 25.))})
        path.write_text("\n".join(lines) + "\n")
        paths.append(str(path))
    cfg = {"Project Name": str(folder / "fit"), "datasets": paths,
           "residues": [{"name": r[0], "flag": "on"} for r in residues],
           "init": {"Method": "Sideband", "kex": {"min": 300., "max": 300., "nsteps": 1},
                    "pB": {"min": .05, "max": .05, "nsteps": 1},
                    "initial": initial, "max_nfev": 100},
           "sideband": {"decoupling": {"h_larmor_mhz": 600., "h_carrier_ppm": 8.5,
                                         "p90_s": 70e-6},
                        "residues": {r[0]: {"h_ppm_a": r[3], "h_ppm_b": r[4]} for r in residues},
                        "v1n": {"mode": "scale", "initial": 1., "bounds": [.8, 1.2]}}}
    model = SidebandModel(cfg)
    p = model.prepare_fit()
    model._fit_data = model._prepare_data()
    return model, p, folder


def check_dense_fallback():
    model = QuadraticModel()
    q, lo, hi, expand, _, indices = _constrained_coordinates(
        model, model.initial_parameters, "kex", 10.)
    assert _profile_jacobian(model, lambda q: model.errFunc(expand(q)), indices, lo, hi) == "3-point"


def check_real_derivatives():
    model, p, folder = real_model()
    original_free = model.free.copy()
    original_lower, original_upper = model.lower.copy(), model.upper.copy()
    cases = [
        (original_free, "kex", 300.),
        (original_free[::-1], "pB", .05),
        (np.array([23, 0, 5, 16, 2, 11, 3]), "kex", 300.),
        (np.array([23, 1, 5, 16, 2, 11, 3]), "pB", .05),
        (np.array([23, 5, 3]), "kex", 300.),
        (original_free[::-1], "v1n_scale", 1.),
    ]
    for free, name, target in cases:
        model.free = free
        model.lower, model.upper = original_lower.copy(), original_upper.copy()
        # Same local group needs inward lower/upper stencils at different residues.
        model.lower[7], model.upper[14] = p[7], p[14]
        q, lo, hi, expand, _, indices = _constrained_coordinates(model, p, name, target)
        calls = []

        def residual(q):
            assert np.all(q >= lo) and np.all(q <= hi)
            calls.append(q.copy())
            return model.errFunc(expand(q))

        jacobian = _profile_jacobian(model, residual, indices, lo, hi)
        assert callable(jacobian), (name, free)
        started = time.perf_counter()
        actual = jacobian(q)
        grouped_time, grouped_calls = time.perf_counter() - started, len(calls)
        calls.clear()
        started = time.perf_counter()
        reference = approx_derivative(residual, q, method="3-point", rel_step=1e-5,
                                      bounds=(lo, hi))
        dense_time, dense_calls = time.perf_counter() - started, len(calls)
        indices = np.asarray(indices)
        groups = np.count_nonzero(indices < 3) + len(set((indices[indices >= 3]-3) % 7))
        assert grouped_calls == 1 + 2*groups, (grouped_calls, groups)
        assert dense_calls == 1 + 2*len(q), (dense_calls, len(q))
        assert grouped_calls <= dense_calls
        # Independent finite differences agree on every affected residual block;
        # grouped zero blocks avoid cancellation noise in inward dense stencils.
        assert_allclose(actual, reference, rtol=1e-8, atol=2e-7)
        if len(q) > 20:
            assert grouped_calls < dense_calls / 2
        print(f"PASS: real {name}, {len(q)} nuisance coordinates, calls "
              f"{dense_calls}->{grouped_calls}, seconds {dense_time:.3f}->{grouped_time:.3f}")
    print(f"Synthetic three-residue inputs retained at {folder}")


def check_profile_linear_algebra_order():
    # Equal Jacobian entries in C/F order can give different column reductions,
    # gradients and weak-H optimizer trajectories. Match SciPy's dense layout.
    model, p, _ = real_model()
    q, lo, hi, expand, _, indices = _constrained_coordinates(model, p, "kex", 295.)
    residual = lambda values: model.errFunc(expand(values))
    actual = _profile_jacobian(model, residual, indices, lo, hi)(q)
    reference = approx_derivative(residual, q, method="3-point", rel_step=1e-5, bounds=(lo, hi))
    assert_array_equal(actual, reference)
    assert actual.flags.f_contiguous and reference.flags.f_contiguous
    assert_array_equal(actual.T @ residual(q), reference.T @ residual(q))
    assert_array_equal(np.sum(actual**2, axis=0), np.sum(reference**2, axis=0))


def check_dense_grouped_profile_solutions():
    model, p, _ = real_model()
    # Enough nuisance flexibility for a real refit, with well-identified locals.
    model.free = np.array([17, 0, 10, 1, 3, 2])
    settings = {"kex": [300., 320.], "pB": [.05, .06], "v1n_scale": [1., 1.02]}
    model.errFunc(p)
    grouped_calls = []
    original = model.errFunc

    def counted(values):
        grouped_calls.append(1)
        return original(values)

    with patch.object(model, "errFunc", side_effect=counted):
        started = time.perf_counter()
        grouped = profile_likelihood(model, p, settings)
        grouped_time, count = time.perf_counter() - started, len(grouped_calls)
        grouped_calls.clear()
        with patch("sb_analysis._profile_jacobian", return_value="3-point"):
            started = time.perf_counter()
            dense = profile_likelihood(model, p, settings)
            dense_time, dense_count = time.perf_counter() - started, len(grouped_calls)
    for name in settings:
        for actual, expected in zip(grouped[name], dense[name]):
            assert actual["success"] and expected["success"], (name, actual, expected)
            assert_array_equal(actual["parameters"], expected["parameters"])
            assert actual["chi2"] == expected["chi2"]
            assert actual["nfev"] == expected["nfev"]
            assert_array_equal(actual["start"], expected["start"])
    assert count < dense_count, (count, dense_count)
    print(f"PASS: real profile solutions, residual calls {dense_count}->{count}, "
          f"seconds {dense_time:.3f}->{grouped_time:.3f}")


if __name__ == "__main__":
    check_dense_fallback()
    check_real_derivatives()
    check_profile_linear_algebra_order()
    check_dense_grouped_profile_solutions()
