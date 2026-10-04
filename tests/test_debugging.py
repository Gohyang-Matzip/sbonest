"""Regression checks for rejected fits, configuration paths and web inputs."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]  # repository root: flat modules live there
sys.path.insert(0, str(ROOT))

import json
import io
import multiprocessing
import subprocess
import tempfile
from unittest.mock import patch

import numpy as np
from scipy.optimize import least_squares

import fit
from estmodel import est_model
from mcrun import run_single_mc_iteration
from run import load_config
from test_performance import CONF, ROOT, example_model


def expect_error(call, kinds=(ValueError, RuntimeError)):
    try:
        call()
    except kinds:
        return
    raise AssertionError("Invalid fit/input was accepted")


def check_fit_failures():
    m = example_model()
    p = fit.generate_initial_parameters(m, CONF)
    expect_error(lambda: m.fit(p0=p[:-1], fitting_config=CONF))
    expect_error(lambda: m.fit(p0=np.r_[-1.0, p[1:]], fitting_config=CONF))
    with patch.object(
        fit,
        "least_squares",
        side_effect=lambda *a, **k: least_squares(*a, **k, max_nfev=1),
    ):
        expect_error(lambda: m.fit(p0=p, fitting_config=CONF), (RuntimeError,))
        conf = {
            "datasets": [
                str(ROOT / "example/syn10.txt"),
                str(ROOT / "example/syn100.txt"),
            ],
            "residues": [],
            "init": CONF,
        }
        assert run_single_mc_iteration((conf, p, 0)) is None
    for r in m.dataset.res:
        r.active = False
    expect_error(lambda: m.fit(fitting_config=CONF))
    expect_error(lambda: m.selMethod({"Method": "Baldwn"}))


def check_noisy_initialization_and_stats():
    m = example_model()
    m.dataset.res[0].estSpecs[0].int[0] = 1.01
    p = fit.generate_initial_parameters(m, CONF)
    assert p[4] >= 0, p[4]
    p, c = m.fit(fitting_config=CONF)
    assert np.isfinite(p).all()
    reported = m.chi2
    actual = float(np.sum(m.errFunc(p) ** 2))
    assert reported == actual, (reported, actual)
    m.dataset.res[0].estSpecs[0].intstd[0] = np.inf
    expect_error(lambda: m.fit(fitting_config=CONF))


def check_loader_errors():
    m = est_model()
    expect_error(
        lambda: m.dataset.addData("missing-onest-regression-file.txt"),
        (FileNotFoundError,),
    )
    conf = {
        "datasets": ["missing-onest-regression-file.txt"],
        "residues": [],
        "init": CONF,
    }
    assert run_single_mc_iteration((conf, np.zeros(7), 0)) is None
    with multiprocessing.get_context("spawn").Pool(1) as pool:
        result = pool.apply_async(run_single_mc_iteration, ((conf, np.zeros(7), 0),))
        assert result.get(timeout=10) is None


def check_loader_structure():
    """Never silently discard observations from an incomplete data block."""
    from est_data import EstDataSet

    folder = Path(tempfile.mkdtemp(prefix="onest-loader-check-"))
    header = "120\n0.4\n25 0\n# offset intensity error\n"
    first = "# A1\n119 0.8 0.01\n120 0.7 0.01\n"
    second = "# A2\n121 0.9 0.02\n"
    path = folder / "valid.txt"
    path.write_text(header + "# file comment\n" + first + "\n# separator\n" + second)
    data = EstDataSet()
    data.addData(str(path))
    assert [r.label for r in data.res] == ["A1", "A2"]
    assert data.res[0].estSpecs[0].offset == [119.0, 120.0]
    assert data.res[0].estSpecs[0].int == [0.8, 0.7]
    assert data.res[1].estSpecs[0].intstd == [0.02]

    invalid = {
        "comment_in_block": header + first + "# note\n121 0.9 0.01\n",
        "missing_column_header": "120\n0.4\n25 0\n" + first + second,
        "data_as_column_header": "120\n0.4\n25 0\n119 0.8 0.01\n" + second,
        "malformed_residue": header + first + "# A2 R2a: bad R2b: 15 dw: 3\n121 0.9 0.01\n",
        "orphan_data": header + "119 0.8 0.01\n" + second,
        "unexpected_text": header + "invalid input\n" + first,
        "empty_dataset": header,
        "empty_residue": header + first + "# A2\n",
    }
    accepted = []
    for name, contents in invalid.items():
        path = folder / f"{name}.txt"
        path.write_text(contents)
        try:
            EstDataSet().addData(str(path))
        except ValueError as exc:
            assert str(path) in str(exc), exc
        else:
            accepted.append(name)
    assert not accepted, f"Loader silently accepted incomplete inputs: {accepted}"


def check_loader_formats_and_noise():
    """Preserve header defaults, file ordering, and Monte Carlo random draws."""
    from est_data import EstDataSet

    folder = Path(tempfile.mkdtemp(prefix="onest-loader-formats-"))
    full = folder / "full.txt"
    full.write_text(
        "1.2e2 # MHz\n4e-1 # seconds\n2.5e1 2 # Hz\n# columns\n"
        "# A1 R2a: 12 R2b: 18 dw: -3 # initial values\n"
        "119 0.8 0.01 # point\n\n120 0.7 0.02\n"
        "# between residues\n# G2 R2a: 14 R2b: 20 dw: 2\n"
        "121 0.9 0.03\n# end\n"
    )
    simple = folder / "simple.txt"
    simple.write_text(
        "80\n0.3\n50 0\n# columns\n# G2\n121 0.6 0.01\n"
        "# A1\n120 0.5 0.01\n"
    )
    state = np.random.get_state()
    try:
        np.random.seed(4321)
        expected = np.random.randn(5)
        np.random.seed(4321)
        data = EstDataSet()
        data.addData(str(full), add_error_to_intensity=True, add_error_to_v1=True)
        assert np.random.randn() == expected[4]
    finally:
        np.random.set_state(state)
    data.addData(str(simple))
    assert [r.label for r in data.res] == ["A1", "G2"]
    assert data.fields == [120.0, 80.0] and data.Ts == [0.4, 0.3]
    assert data.v1s == [25 + 2 * expected[0], 50.0]
    assert data.v1errs == [2.0, 0.0] and data.initR2 is True
    a, g = [r.estSpecs[0] for r in data.res]
    assert (a.initr2a, a.initr2b, a.initdw) == (12.0, 18.0, -3.0)
    assert (g.initr2a, g.initr2b, g.initdw) == (14.0, 20.0, 2.0)
    assert a.v1 == g.v1 == data.v1s[0]
    np.testing.assert_array_equal(a.int, [0.8 + 0.01 * expected[1], 0.7 + 0.02 * expected[2]])
    assert g.int == [0.9 + 0.03 * expected[3]]
    for residue, intensity in zip(data.res, (0.5, 0.6)):
        spec = residue.estSpecs[1]
        assert spec.int == [intensity]
        assert (spec.initr2a, spec.initr2b, spec.initdw) == (10.0, 100.0, 0.1)
    for body in (
        "# A1\n119 0.8 0.01\n# G2 R2a: 14 R2b: 20 dw: 2\n121 0.9 0.01\n",
        "# A1 R2a: 12 R2b: 18 dw: -3\n119 0.8 0.01\n# G2\n121 0.9 0.01\n",
    ):
        invalid = folder / "mixed.txt"
        invalid.write_text("120\n0.4\n25 0\n# columns\n" + body)
        expect_error(lambda: EstDataSet().addData(str(invalid)), (ValueError,))


def check_config_and_benchmark():
    folder = Path(tempfile.mkdtemp(prefix="onest-config-check-"))
    (folder / "sample.txt").write_bytes((ROOT / "example/syn10.txt").read_bytes())
    conf = {
        "Project Name": "paths",
        "datasets": ["sample.txt"],
        "residues": [{"name": f"A{i}", "flag": "on"} for i in range(1, 6)],
        "init": CONF,
    }
    path = folder / "config.json"
    path.write_text(json.dumps(conf))
    loaded = load_config(str(path))
    assert Path(loaded["datasets"][0]).resolve() == (folder / "sample.txt").resolve()
    for script in ["run.py", "benchmark.py"]:
        args = ["--no-pdf"] if script == "run.py" else []
        out = subprocess.run(
            [sys.executable, str(ROOT / script), str(path), *args],
            cwd=folder.parent,
            capture_output=True,
            text=True,
            timeout=30,
        )
        assert out.returncode == 0, out.stdout + out.stderr
    conf["residues"] = [{"name": "MISSING", "flag": "on"}]
    path.write_text(json.dumps(conf))
    out = subprocess.run(
        [sys.executable, str(ROOT / "benchmark.py"), str(path)],
        cwd=folder,
        capture_output=True,
        text=True,
        timeout=30,
    )
    assert out.returncode != 0, "benchmark ignored invalid residue"


def check_web_inputs():
    import server_run

    folder = Path(tempfile.mkdtemp(prefix="onest-http-check-"))
    job = folder / "case"
    job.mkdir()
    conf = {
        "Project Name": "web",
        "datasets": [str(ROOT / "example/syn10.txt")],
        "residues": [{"name": f"A{i}", "flag": "on"} for i in range(1, 6)],
        "init": CONF,
    }
    (job / "initial_config.json").write_text(json.dumps(conf))
    (job / "config.json").write_text(json.dumps(conf))
    app = server_run.app
    app.config.update(JOBS_DIR=str(folder), TESTING=True)
    client = app.test_client()
    for route in ["run_fit", "run_mc"]:
        assert client.post(f"/{route}/missing").status_code == 404
    valid = {
        "project_name": "web",
        "method": "Baldwin",
        "residue_flags": conf["residues"],
    }
    for payload in [
        None,
        [],
        {**valid, "method": "typo"},
        {**valid, "residue_flags": [{"name": "A1", "flag": "bad"}]},
        {
            **valid,
            "residue_flags": [{"name": f"A{i}", "flag": "off"} for i in range(1, 6)],
        },
    ]:
        result = client.post(
            "/finalize_config/case",
            data=json.dumps(payload),
            content_type="application/json",
        )
        assert result.status_code == 400, (payload, result.status_code, result.json)
    assert client.post("/finalize_config/case", json=valid).status_code == 200
    for value in [None, [], {}, 1.5, True, 0, -1]:
        result = client.post("/run_mc/case", json={"num_runs": value})
        assert result.status_code == 400, (value, result.status_code)
    for value in [[], {}, 1.5, True, 0, -1]:
        result = client.post(
            "/run_mc/case", json={"num_runs": 1, "num_processes": value}
        )
        assert result.status_code == 400, (value, result.status_code)
    invalid = {
        **valid,
        "residue_flags": [{"name": f"A{i}", "flag": "off"} for i in range(1, 6)],
    }
    conf["residues"] = invalid["residue_flags"]
    (job / "config.json").write_text(json.dumps(conf))
    result = client.post("/run_fit/case")
    assert result.status_code == 200 and result.json["success"] is False, result.json
    assert result.json["output_files"] == []
    result = client.post("/run_mc/case", json={"num_runs": 1})
    assert result.status_code == 200 and result.json["success"] is False, result.json
    assert result.json["output_files"] == []
    del conf["Project Name"]
    (job / "config.json").write_text(json.dumps(conf))
    assert client.post("/run_fit/case").status_code == 500
    assert client.post("/run_mc/case", json={"num_runs": 1}).status_code == 500
    print("Web artifacts:", folder)


def check_web_roundtrip():
    import server_run

    folder = Path(tempfile.mkdtemp(prefix="onest-web-roundtrip-"))
    app = server_run.app
    app.config.update(JOBS_DIR=str(folder), TESTING=True)
    client = app.test_client()
    response = client.post(
        "/prepare_initial",
        data={
            "data_files[]": [
                (io.BytesIO((ROOT / "example" / name).read_bytes()), name)
                for name in ["syn10.txt", "syn100.txt"]
            ]
        },
    )
    assert response.status_code == 200, response.json
    job = response.json["job_id"]
    config = response.json["initial_config_json"]
    response = client.post(
        "/finalize_config/" + job,
        json={
            "project_name": "한글",
            "method": "Baldwin",
            "residue_flags": config["residues"],
        },
    )
    assert response.status_code == 200, response.json
    for route, payload, names in [
        ("run_fit", None, ["한글_data.pdf", "한글_result.txt", "한글.pdf"]),
        (
            "run_mc",
            {"num_runs": 2, "num_processes": 2},
            ["한글_data.pdf", "한글.pdf", "한글_mc.txt", "한글_mcmean.pdf"],
        ),
    ]:
        response = client.post(f"/{route}/{job}", json=payload)
        assert response.json["success"], response.json
        assert [entry["name"] for entry in response.json["output_files"]] == names
        for entry in response.json["output_files"]:
            download = client.get(entry["url"])
            assert download.status_code == 200, (entry, download.status_code)
            assert download.data
    duplicate = client.post(
        "/prepare_initial",
        data={
            "data_files[]": [
                (io.BytesIO((ROOT / "example/syn10.txt").read_bytes()), "same.txt"),
                (io.BytesIO((ROOT / "example/syn100.txt").read_bytes()), "same.txt"),
            ]
        },
    )
    assert duplicate.status_code == 400, duplicate.json
    assert client.get("/download/%2e%2e/README.md").status_code == 400
    assert client.get(f"/download/{job}/%2e%2e/initial_config.json").status_code == 404
    print("Web roundtrip artifacts:", folder)


if __name__ == "__main__":
    failed = []
    for check in [
        check_fit_failures,
        check_noisy_initialization_and_stats,
        check_loader_errors,
        check_loader_structure,
        check_loader_formats_and_noise,
        check_config_and_benchmark,
        check_web_inputs,
        check_web_roundtrip,
    ]:
        try:
            check()
            print("PASS", check.__name__, flush=True)
        except (Exception, SystemExit) as exc:
            failed.append(check.__name__)
            print("FAIL", check.__name__, repr(exc), flush=True)
    assert not failed, failed
