"""Run with: OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 python3 test_performance.py."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]  # repository root: flat modules live there
sys.path.insert(0, str(ROOT))

import copy
import json
import multiprocessing
import os
import subprocess
import tempfile
from unittest.mock import patch

import numpy as np
from scipy.linalg import expm

import estmodel
import fit
from estmodel import est_model, PARAM_LAYOUT

CONF = {
    "Method": "Baldwin",
    "kex": {"min": 10.0, "max": 400.0, "nsteps": 6},
    "pB": {"min": 0.01, "max": 0.1, "nsteps": 6},
}


def example_model():
    m = est_model()
    for name in ["syn10.txt", "syn100.txt"]:
        m.dataset.addData(str(ROOT / "example" / name))
    return m


def check_initial_parameters():
    """Keep parameter order and defaults stable across all five methods."""
    m = example_model()
    m.dataset.res = m.dataset.res[:2]
    m.dataset.res[1].active = False
    spectrum = m.dataset.res[0].estSpecs[0]
    m.dataset.res[0].estSpecs = [spectrum]
    spectrum.offset, spectrum.int = [3.0, 1.0, 2.0], [0.8, 0.4, 0.6]
    spectrum.intstd, spectrum.T = [0.01] * 3, 0.4
    spectrum.initdw, spectrum.initr2a, spectrum.initr2b = -2.0, 11.0, 22.0
    r1 = -np.log(0.8) / 0.4
    expected = {
        "NoEx": ([1.0, r1, 11.0], [0.0, 1.0, 10.0]),
        "Baldwin": (
            [5.0, 95.0, 1.0, -2.0, r1, 11.0, 22.0],
            [5.0, 95.0, 0.0, 0.1, 1.0, 10.0, 20.0],
        ),
        "Matrix_3st_Linear": (
            [10.0] * 4 + [1.0, -2.0, -4.0, r1, 11.0, 22.0, 22.0],
            [10.0] * 4 + [0.0, 0.1, 0.2, 1.0, 10.0, 20.0, 20.0],
        ),
        "Matrix_3st_Triangle": (
            [10.0] * 6 + [1.0, -2.0, -4.0, r1, 11.0, 22.0, 22.0],
            [10.0] * 6 + [0.0, 0.1, 0.2, 1.0, 10.0, 20.0, 20.0],
        ),
    }
    expected["Matrix"] = expected["Baldwin"]
    for method, (observed, defaults) in expected.items():
        conf = {
            "Method": method,
            "kex": {"min": 100.0, "max": 100.0, "nsteps": 1},
            "pB": {"min": 0.05, "max": 0.05, "nsteps": 1},
        }
        for spectra, target in [([spectrum], observed), ([], defaults)]:
            m.dataset.res[0].estSpecs = spectra
            np.testing.assert_array_equal(
                fit.generate_initial_parameters(m, conf), target
            )


def check_fit_and_jacobian():
    m = example_model()
    with patch.object(m, "errFunc", wraps=m.errFunc) as residual:
        p, covariance = m.fit(fitting_config=CONF)
        assert residual.call_count < 120, residual.call_count
    assert np.isfinite(covariance).all()
    assert abs(np.sum(m.errFunc(p) ** 2) - 2080.693104493508) < 1e-5

    for method, (globals_, locals_, _) in PARAM_LAYOUT.items():
        m = example_model()
        m.method = method
        m.dataset.res[1].active = False
        m.dataset.res[3].active = False
        for i, r in enumerate(m.dataset.res):
            for es in r.estSpecs:
                es.offset, es.int, es.intstd = [
                    v[:: 30 + i] for v in (es.offset, es.int, es.intstd)
                ]
        p = fit.generate_initial_parameters(
            m,
            {
                **CONF,
                "Method": method,
                "kex": {"min": 100.0, "max": 100.0, "nsteps": 1},
                "pB": {"min": 0.05, "max": 0.05, "nsteps": 1},
            },
        )
        # Exercise a free negative shift and a rate exactly at its lower bound.
        p[len(globals_)] = -0.5
        p[-1] = 0.0
        sizes = [
            sum(len(es.offset) for es in r.estSpecs) for r in m.dataset.res if r.active
        ]
        assert hasattr(fit, "_block_jacobian"), "missing grouped derivative"
        jac = fit._block_jacobian(m.errFunc, sizes, len(globals_), len(locals_))
        actual = jac(p)
        base = m.errFunc(p)
        reference = np.empty_like(actual)
        for col in range(len(p)):
            q = p.copy()
            q[col] += (
                np.sqrt(np.finfo(float).eps)
                * (1 if p[col] >= 0 else -1)
                * max(1.0, abs(p[col]))
            )
            reference[:, col] = (m.errFunc(q) - base) / (q[col] - p[col])
        np.testing.assert_allclose(actual, reference, rtol=1e-9, atol=1e-9)


def scalar_matrix(args):
    kab, kba, dG, dw, R1, R2a, R2b, offsets, T, v1, v1err, B0 = args
    pB = kab / (kab + kba) if kab + kba else 0.0
    pA = 1.0 - pB
    wxs = (
        2
        * np.pi
        * (np.linspace(v1 - 2 * v1err, v1 + 2 * v1err, 10) if v1err else np.array([v1]))
    )
    weights = (
        np.exp(-0.5 * ((wxs / (2 * np.pi) - v1) / v1err) ** 2) if v1err else np.ones(1)
    )
    weights /= weights.sum()
    result = []
    for offset in offsets:
        intensity = 0.0
        for wx, w in zip(wxs, weights):
            a = np.diag(
                [-R2a - kab, -R2a - kab, -R1 - kab, -R2b - kba, -R2b - kba, -R1 - kba]
            )
            a[:3, 3:] = kba * np.eye(3)
            a[3:, :3] = kab * np.eye(3)
            for start, shift in [(0, dG), (3, dG + dw)]:
                wa = (shift - offset) * B0 * 2 * np.pi
                a[start, start + 1] = -wa
                a[start + 1, start] = wa
                a[start + 1, start + 2] = -wx
                a[start + 2, start + 1] = wx
            end = expm(a * T) @ np.array([0.0, 0.0, pA, 0.0, 0.0, pB])
            intensity += w * end[2] / (pA or 1.0)
        result.append(max(0.0, intensity))
    return result


def check_matrix_and_weights():
    for rates in [(15.0, 285.0), (0.0, 0.0), (0.0, 100.0), (100.0, 0.0)]:
        for error in [0.0, 0.4]:
            args = (
                *rates,
                0.0,
                -2.0,
                1.0,
                10.0,
                25.0,
                [-8.0, -2.0, 0.0, 0.5, 3.0],
                0.4,
                20.0,
                error,
                80.12,
            )
            np.testing.assert_allclose(
                estmodel.matrix_calc_worker_chunk(args),
                scalar_matrix(args),
                rtol=1e-10,
                atol=1e-11,
            )
    args = (
        15.0,
        285.0,
        0.0,
        -2.0,
        1.0,
        10.0,
        25.0,
        np.linspace(-5, 5, 205),
        0.4,
        20.0,
        0.4,
        80.12,
    )
    with patch.object(estmodel, "expm", wraps=estmodel.expm) as exponential:
        values = estmodel.matrix_calc_worker_chunk(args)
        assert len(values) == 205 and np.isfinite(values).all()
        assert exponential.call_count <= 2, exponential.call_count
    wxs, weights = estmodel.b1_weights(20.0, 0.4)
    assert not wxs.flags.writeable and not weights.flags.writeable
    assert estmodel.b1_weights(20.0, 0.4)[0] is wxs
    assert estmodel.matrix_calc_worker_chunk((*args[:7], [], *args[8:])) == []


def check_log_and_repeated_fit():
    m = example_model()
    p, _ = m.fit(fitting_config=CONF)
    P = m.seParam(p)
    expected = []
    for i, r in enumerate(m.dataset.res):
        for es in r.estSpecs:
            for off, obs, std in zip(es.offset, es.int, es.intstd):
                calc = np.nan_to_num(m.calc(P, i, off, es))
                expected.append(f"{off:8.3f} {obs:12.3f} {std:12.3f} {calc:12.3f}")
    with patch.object(m, "calc", wraps=m.calc) as calc:
        lines = m._log_data_tables(P)
        assert calc.call_count == 10, calc.call_count
    expected_lines = set(expected)
    assert [line for line in lines if line in expected_lines] == expected
    m.dataset.res[1].active = False
    m.dataset.res[0].estSpecs[0].int[0] += 0.01
    m.dataset.res[0].estSpecs[0].intstd[0] = 0.0
    fresh = est_model()
    fresh.dataset = copy.deepcopy(m.dataset)
    p1, c1 = m.fit(fitting_config=CONF)
    p2, c2 = fresh.fit(fitting_config=CONF)
    np.testing.assert_allclose(p1, p2, rtol=1e-10, atol=1e-10)
    np.testing.assert_allclose(c1, c2, rtol=1e-10, atol=1e-10)
    assert m.dataset.res[0].estSpecs[0].intstd[0] == 0.0


def matrix_in_worker(_):
    m = example_model()
    m.dataset.res = m.dataset.res[:2]
    p = np.array(
        [15.0, 285.0, 118.0, -5.0, 1.0, 10.0, 20.0, 110.0, -1.0, 1.0, 10.0, 20.0]
    )
    m.method = "Matrix"
    P = m.seParam(p)
    for i, r in enumerate(m.dataset.res):
        for es in r.estSpecs:
            es.int = m.calc(P, i, np.array(es.offset), es).tolist()
    initial = p.copy()
    initial[0] *= 1.2
    out, _ = m.fit(p0=initial, fitting_config={"Method": "Matrix"})
    return np.sum(m.errFunc(out) ** 2)


def check_mc_worker():
    with multiprocessing.get_context("spawn").Pool(1) as pool:
        chi2 = pool.map(matrix_in_worker, [None])[0]
    assert chi2 < 1e-6, chi2


def check_cli():
    # Retain outputs in a temporary directory for inspection; do not delete artifacts.
    folder = Path(tempfile.mkdtemp(prefix="onest-check-"))
    conf = {
        "Project Name": "smoke",
        "datasets": [str(ROOT / "example/syn10.txt")],
        "residues": [{"name": f"A{i}", "flag": "on"} for i in range(1, 6)],
        "init": CONF,
    }
    (folder / "config.json").write_text(json.dumps(conf))
    env = {**os.environ, "MPLBACKEND": "Agg", "PYTHONUNBUFFERED": "1"}
    # A fixed seed keeps the two noisy refits reproducible; unseeded draws made
    # this smoke test fail intermittently when a refit hit the evaluation limit.
    for script, args in [("run.py", []), ("mcrun.py", ["2", "2", "--seed", "20261004"])]:
        result = subprocess.run(
            [sys.executable, str(ROOT / script), "config.json", *args, "--no-pdf"],
            cwd=folder,
            env=env,
            capture_output=True,
            text=True,
            timeout=90,
        )
        assert result.returncode == 0, result.stdout + result.stderr
    assert (folder / "smoke_result.txt").is_file() and (
        folder / "smoke_mc.txt"
    ).is_file()
    assert not list(folder.glob("*.pdf")), folder
    result = subprocess.run(
        [sys.executable, str(ROOT / "run.py"), "config.json"],
        cwd=folder,
        env=env,
        capture_output=True,
        text=True,
        timeout=90,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    for name in ["smoke.pdf", "smoke_data.pdf"]:
        assert (folder / name).read_bytes().startswith(b"%PDF-")

    conf["Project Name"] = "matrix"
    conf["init"] = {**CONF, "Method": "Matrix"}
    # Two B1 conditions constrain the rates and avoid a slow single-profile fit.
    conf["datasets"] = [
        str(ROOT / "example" / name) for name in ["syn10.txt", "syn100.txt"]
    ]
    (folder / "matrix.json").write_text(json.dumps(conf))
    result = subprocess.run(
        [sys.executable, str(ROOT / "mcrun.py"), "matrix.json", "2", "2", "--no-pdf"],
        cwd=folder,
        env=env,
        capture_output=True,
        text=True,
        timeout=90,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "Completed 2/2 MC runs successfully." in result.stdout
    assert (folder / "matrix_mc.txt").is_file()
    assert not list(folder.glob("matrix*.pdf"))
    print("CLI artifacts:", folder)


if __name__ == "__main__":
    checks = [
        check_initial_parameters,
        check_fit_and_jacobian,
        check_matrix_and_weights,
        check_log_and_repeated_fit,
        check_mc_worker,
        check_cli,
    ]
    failed = []
    for check in checks:
        try:
            check()
            print("PASS", check.__name__, flush=True)
        except AssertionError as exc:
            failed.append(check.__name__)
            print("FAIL", check.__name__, str(exc), flush=True)
        except subprocess.TimeoutExpired as exc:
            print("Timed out:", exc.cmd, flush=True)
            for output in (exc.stdout, exc.stderr):
                if output:
                    print(
                        output.decode(errors="replace")
                        if isinstance(output, bytes)
                        else output,
                        flush=True,
                    )
            raise
    assert not failed, failed
