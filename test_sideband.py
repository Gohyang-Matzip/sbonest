"""Runnable physical checks: OC integration, exact timing and RF-field recovery."""

# ruff: noqa: E402 -- Set BLAS limits before importing NumPy.
import importlib.util
import os

for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(name, "1")

import numpy as np
from numpy.testing import assert_allclose
from scipy.linalg import expm


def reference(
    segments,
    offsets,
    *,
    T,
    nu,
    kab,
    kba,
    dw,
    r1,
    r2a,
    r2b,
    ha,
    hb,
    J=92.0,
    r1h=2.0,
    r2h=25.0,
):
    """Independent complex density matrix, H first; no production helpers."""
    eye = np.eye(2)
    spin = [
        eye,
        np.array([[0, 1], [1, 0]]),
        np.array([[0, -1j], [1j, 0]]),
        np.diag([1, -1]),
    ]
    H = [np.kron(s / 2, eye) for s in spin[1:]]
    N = [np.kron(eye, s / 2) for s in spin[1:]]

    def vec(a):
        return a.ravel(order="F")

    def comm(a):
        return -1j * (np.kron(np.eye(4), a) - np.kron(a.T, np.eye(4)))

    def relaxation(r2):
        out = np.zeros((16, 16), complex)
        for i, hi in enumerate(spin):
            for j, nj in enumerate(spin):
                b = vec(np.kron(hi, nj) / 2)
                rate = [0, r2h, r2h, r1h][i] + [0, r2, r2, r1][j]
                out += rate * np.outer(b, b.conj())
        return out

    pa = kba / (kab + kba)
    initial = np.r_[pa * vec(N[2]), (1 - pa) * vec(N[2])]
    period = sum(row[0] for row in segments)
    count = int(np.floor(T / period + 1e-12))
    values = []
    for off in offsets:
        full, partial = np.eye(32, dtype=complex), np.eye(32, dtype=complex)
        remaining = max(0.0, T - count * period)
        for duration, fx, fy in segments:
            L = np.zeros((32, 32), complex)
            for j, (om, dh, r2, rate) in enumerate(
                [(-off, ha, r2a, kab), (dw - off, hb, r2b, kba)]
            ):
                h = (
                    2
                    * np.pi
                    * (
                        om * N[2]
                        + nu * N[0]
                        + dh * H[2]
                        + fx * H[0]
                        + fy * H[1]
                        + J * N[2] @ H[2]
                    )
                )
                b = slice(16 * j, 16 * (j + 1))
                L[b, b] = comm(h) - relaxation(r2) - rate * np.eye(16)
            L[:16, 16:] = kba * np.eye(16)
            L[16:, :16] = kab * np.eye(16)
            full = expm(L * duration) @ full
            use = min(duration, remaining)
            if use > 1e-14:
                partial = expm(L * use) @ partial
                remaining -= use
        state = partial @ np.linalg.matrix_power(full, count) @ initial
        values.append((vec(N[2]).conj() @ state[:16]).real / pa)
    return np.array(values)


def check_automatic_proton_rates(folder):
    """Catch shared/wrong peak H rates, broken packing, and silently fixed defaults."""
    import copy

    from sbfit import SidebandModel
    from sideband import composite_segments

    seg = composite_segments(70e-6, "RR")
    ppm = np.linspace(105, 135, 147)
    peaks = {
        "A1": [120.0, 3.0, 1.5, 12.0, 15.0, 2.0, 25.0],
        "S3": [121.5, -2.5, 1.7, 16.0, 25.0, 3.0, 35.0],
    }
    h = {"A1": [6.2, 6.5], "S3": [9.8, 9.6]}
    paths = []
    for rf in (25.0, 100.0):
        lines = [f"121.5949416\n0.4\n{rf} 0\n# offset intensity error\n"]
        for name in ("A1", "G2", "S3"):
            key = "A1" if name == "G2" else name
            shift, dw, r1, r2a, r2b, r1h, r2h = peaks[key]
            y = reference(
                seg,
                (ppm - shift) * 121.5949416,
                T=0.4,
                nu=rf * 1.08,
                kab=15.0,
                kba=285.0,
                dw=dw * 121.5949416,
                r1=r1,
                r2a=r2a,
                r2b=r2b,
                ha=(h[key][0] - 8.5) * 1200,
                hb=(h[key][1] - 8.5) * 1200,
                r1h=r1h,
                r2h=r2h,
            )
            lines.append(f"# {name} R2a: {r2a} R2b: {r2b} dw: {dw}\n")
            lines.extend(f"{x:.14g} {v:.14g} 0.001\n" for x, v in zip(ppm, y))
        path = folder / f"auto_h_{rf:g}.txt"
        path.write_text("".join(lines))
        paths.append(str(path))
    cfg = {
        "Project Name": "auto_h_check",
        "datasets": paths,
        "residues": [
            {"name": n, "flag": "off" if n == "G2" else "on"}
            for n in ("A1", "G2", "S3")
        ],
        "init": {
            "Method": "Sideband",
            "kex": {"min": 300, "max": 300, "nsteps": 1},
            "pB": {"min": 0.05, "max": 0.05, "nsteps": 1},
            "initial": {
                f"{n}.{k}": v
                for n, values in peaks.items()
                for k, v in zip(("peak_ppm", "dw_ppm", "R1", "R2a", "R2b"), values)
            },
        },
        "sideband": {
            "decoupling": {
                "h_larmor_mhz": 1200.0,
                "h_carrier_ppm": 8.5,
                "p90_s": 70e-6,
            },
            "residues": {
                n: dict(zip(("h_ppm_a", "h_ppm_b"), shifts)) for n, shifts in h.items()
            },
            "v1n": {"mode": "scale", "initial": 1.03, "bounds": [0.8, 1.2]},
        },
    }
    model = SidebandModel(cfg)
    assert all(
        f"{n}.{k}" in model.parameter_names for n in peaks for k in ("R1H", "R2H")
    ), "Omitting H rates must fit them independently for each active peak"
    p, cov = model.fit()
    expected = np.r_[15.0, 285.0, 1.08, peaks["A1"], peaks["S3"]]
    assert_allclose(p, expected, atol=2e-3, rtol=1e-4)
    assert model.chi2 < 1e-8 and model.rank == 17
    assert model.nvar == 588 and model.dof == 571
    assert np.isfinite(cov).all() and np.all(np.diag(cov) > 0)
    info = model.diagnostics(p, cov)
    assert info["proton_relaxation_mode"] == "fit"
    assert "A1.R1H" in model.getLogBuffer((p, cov))
    assert len(info["proton_correlations"]) == 2
    # Explicit fixed settings retain the original vector layout and forward model.
    fixed = copy.deepcopy(cfg)
    fixed["sideband"]["decoupling"].update(R1H=2.0, R2H=25.0)
    legacy = SidebandModel(fixed)
    assert not any(n.endswith(("R1H", "R2H")) for n in legacy.parameter_names)
    fixed["sideband"]["proton_relaxation"] = {"mode": "fixed"}
    del fixed["sideband"]["decoupling"]["R1H"]
    del fixed["sideband"]["decoupling"]["R2H"]
    explicit = SidebandModel(fixed)
    nitrogen = np.r_[15.0, 285.0, 1.08, peaks["A1"][:5], peaks["S3"][:5]]
    assert_allclose(legacy.errFunc(nitrogen), explicit.errFunc(nitrogen), atol=1e-12)
    for case in ("bad_mode", "conflicting_fixed_rates", "bad_nitrogen_mode"):
        bad = copy.deepcopy(cfg)
        if case == "bad_mode":
            bad["sideband"]["proton_relaxation"] = {"mode": "typo"}
        elif case == "conflicting_fixed_rates":
            bad["sideband"]["proton_relaxation"] = {"mode": "fit"}
            bad["sideband"]["decoupling"]["R1H"] = 2.0
        else:
            bad["sideband"]["nitrogen_relaxation"] = {"mode": "typo"}
        try:
            SidebandModel(bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"Invalid automatic proton settings accepted: {case}")
    # Several proton fields with automatic proton rates now give one R1H/R2H pair
    # per field group; single-field names are unchanged.
    two_fields = copy.deepcopy(cfg)
    two_fields["sideband"]["datasets"] = [{}, {"h_larmor_mhz": 800.0}]
    grouped = SidebandModel(two_fields)
    assert grouped.field_groups == [600.0, 800.0] or len(grouped.field_groups) == 2
    assert grouped.proton_groups == 2 and grouped.nitrogen_groups == 1
    labels = [r.label for r in grouped.dataset.res if r.active]
    assert f"{labels[0]}.R1H[0]" in grouped.parameter_names and f"{labels[0]}.R2H[1]" in grouped.parameter_names
    assert f"{labels[0]}.R1" in grouped.parameter_names and f"{labels[0]}.R1[0]" not in grouped.parameter_names
    print(
        "PASS: automatic peakwise proton fitting, inactive peak, reports and legacy fixed mode"
    )


def main():
    assert importlib.util.find_spec("sideband") is not None, (
        "Sideband engine is missing"
    )
    from estmodel import matrix_calc_worker_chunk
    from sideband import composite_segments, profile

    segments = composite_segments(70e-6, "RR")
    args = dict(
        T=0.003173,
        nu=35.0,
        kab=15.0,
        kba=285.0,
        dw=360.0,
        r1=1.5,
        r2a=12.0,
        r2b=18.0,
        ha=-2700.0,
        hb=-2300.0,
    )
    offsets = np.array([-1700.0, -31.0, 0.0, 360.0, 1530.0, 2250.0])
    # Changing RF phase, J factor, tensor order, exchange direction, or remainder fails.
    for cycle in ("RR", "RRbar", "MLEV4"):
        seg = composite_segments(70e-6, cycle)
        for T in (0.0, 31e-6, seg[:, 0].sum() * 2, 0.003173):
            a = dict(args, T=T)
            assert_allclose(
                profile(seg, offsets, **a),
                reference(seg, offsets, **a),
                atol=2e-11,
                rtol=2e-11,
            )
    # J=0 must reduce to ONEST's independent Bloch-McConnell calculation.
    expected = matrix_calc_worker_chunk(
        (
            15.0,
            285.0,
            120.0,
            3.0,
            1.5,
            12.0,
            18.0,
            120 + offsets / 120,
            0.4,
            35.0,
            0.0,
            120.0,
        )
    )
    assert_allclose(
        profile(segments, offsets, **dict(args, T=0.4, J=0.0)),
        expected,
        atol=2e-10,
        rtol=2e-10,
    )
    print("PASS: independent density matrix, cycle/timing, J=0 ONEST parity")

    assert importlib.util.find_spec("sbfit") is not None, "Sideband fitter is missing"
    import tempfile
    from pathlib import Path

    from sbfit import SidebandModel

    folder = Path(tempfile.mkdtemp(prefix="sbonest-check-"))
    check_automatic_proton_rates(folder)
    paths = []
    for index, nu in enumerate((20.0, 50.0)):
        off = np.r_[
            np.linspace(-70, 70, 13),
            [300.0, 360.0, 420.0],
            np.linspace(1500, 1560, 9),
            [-2800.0, 2800.0],
        ]
        intensity = reference(segments, off, **dict(args, T=0.04, nu=nu * 1.08))
        path = folder / f"v1_{index}.txt"
        path.write_text(
            f"120\n0.04\n{nu} 0\n# offset intensity error\n# A1 R2a: 12 R2b: 18 dw: 3\n"
            + "\n".join(
                f"{120 + o / 120:.12g} {y:.12g} 0.001" for o, y in zip(off, intensity)
            )
            + "\n"
        )
        paths.append(str(path))
    config = {
        "Project Name": "check",
        "datasets": paths,
        "residues": [{"name": "A1", "flag": "on"}],
        "init": {
            "Method": "Sideband",
            "kex": {"min": 300, "max": 300, "nsteps": 1},
            "pB": {"min": 0.05, "max": 0.05, "nsteps": 1},
            "initial": {"A1.peak_ppm": 120.0, "A1.R1": 1.5},
            "vary": ["v1n_scale"],
        },
        "sideband": {
            "decoupling": {
                "h_larmor_mhz": 1200.0,
                "h_carrier_ppm": 8.5,
                "p90_s": 70e-6,
                "cycle": "RR",
            },
            "residues": {"A1": {"h_ppm_a": 6.25, "h_ppm_b": 6.583333333333333}},
            "v1n": {"mode": "scale", "initial": 1.0, "bounds": [0.8, 1.2]},
        },
    }
    m = SidebandModel(config)
    p, cov = m.fit(fitting_config=config["init"])
    assert_allclose(p[2], 1.08, atol=2e-6)
    assert m.chi2 < 1e-10 and cov[2, 2] > 0
    assert "v1n_scale" in m.getLogBuffer((p, cov))
    config["sideband"]["v1n"] = {"mode": "per_dataset"}
    config["init"]["vary"] = ["v1n[0]", "v1n[1]"]
    m = SidebandModel(config)
    p, _ = m.fit(fitting_config=config["init"])
    assert_allclose(p[2:4], [21.6, 54.0], atol=1e-4)
    print(
        "PASS: fitted v1n scale and per-dataset Hz against independent synthetic data"
    )

    # OC's native JSON loader must respect channel order, amplitude units and dwell.
    from optimalcontrol.io import Waveform, export_json

    from sideband import waveform_segments

    waveform = folder / "oc.json"
    w = Waveform(
        ["y", "x"],
        "a.u.",
        np.arange(3) * 20e-6,
        np.array([[0.0, 0.5, 0.0], [1.0, 0.0, 1.0]]),
        {"pulse_dt": 20e-6},
        "test",
    )
    export_json(w, waveform)
    seg = waveform_segments(waveform, rf_hz=4000.0)
    assert_allclose(
        seg, [[20e-6, 4000.0, 0.0], [20e-6, 0.0, 2000.0], [20e-6, 4000.0, 0.0]]
    )
    assert_allclose(
        profile(seg, offsets, **args),
        reference(seg, offsets, **args),
        atol=2e-11,
        rtol=2e-11,
    )
    for invalid in (
        dict(config["sideband"]["v1n"], mode="typo"),
        {"mode": "scale", "bounds": [1.0, 1.0]},
        {"mode": "scale", "initial": float("nan")},
    ):
        import copy

        bad = copy.deepcopy(config)
        bad["sideband"]["v1n"] = invalid
        try:
            SidebandModel(bad)
        except ValueError:
            pass
        else:
            raise AssertionError("Invalid RF configuration was accepted")
    assert m.diagnostics(p, np.zeros((len(p), len(p))))["v1n_hz"] == list(p[2:4])
    print("PASS: OC waveform units/channel order and invalid RF inputs")


if __name__ == "__main__":
    main()
