"""Regression checks for rejected fits, configuration paths and web inputs."""

import json
import io
import multiprocessing
from pathlib import Path
import subprocess
import sys
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
