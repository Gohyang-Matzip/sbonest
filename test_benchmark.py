"""Real CLI checks for Sideband/ONEST dispatch and preserved profile artifacts."""

import json
import os
from pathlib import Path
import pstats
import subprocess
import sys
import tempfile

import numpy as np
from optimalcontrol.io import Waveform, export_json

from test_sideband import reference


ROOT = Path(__file__).resolve().parent


def run_cli(folder, config, *args):
    env = dict(os.environ, MPLBACKEND="Agg")
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        env[name] = "1"
    return subprocess.run(
        [sys.executable, str(ROOT / "benchmark.py"), str(config), *args],
        cwd=folder,
        env=env,
        capture_output=True,
        text=True,
        timeout=30,
    )


def write_sideband_config(folder):
    """Tiny independent synthetic spectrum, with data/waveform beside the config."""
    inputs = folder / "inputs"
    inputs.mkdir()
    export_json(
        Waveform(
            ["x", "y"], "Hz", np.array([0.0]),
            np.array([[3500.0], [0.0]]), {"pulse_dt": 70e-6}, "benchmark test",
        ),
        inputs / "waveform.json",
    )
    offsets = np.array([-100.0, -30.0, 0.0, 30.0, 100.0])
    intensity = reference(
        [[70e-6, 3500.0, 0.0]], offsets, T=0.004, nu=25.0,
        kab=15.0, kba=285.0, dw=360.0, r1=1.5, r2a=12.0, r2b=18.0,
        ha=-2700.0, hb=-2300.0,
    )
    (inputs / "data.txt").write_text(
        "120\n0.004\n25 0\n# offset intensity error\n"
        "# A1 R2a: 12 R2b: 18 dw: 3\n"
        + "".join(
            f"{120 + off / 120:.16g} {value:.16g} 0.001\n"
            for off, value in zip(offsets, intensity)
        )
    )
    config = {
        "Project Name": "benchmark",
        "datasets": ["data.txt"],
        "residues": [{"name": "A1", "flag": "on"}],
        "init": {
            "Method": "Sideband",
            "kex": {"min": 300, "max": 300, "nsteps": 1},
            "pB": {"min": 0.05, "max": 0.05, "nsteps": 1},
            "initial": {"A1.peak_ppm": 120.0, "A1.R1": 1.5},
            "vary": ["A1.R1"],
        },
        "sideband": {
            "decoupling": {
                "h_larmor_mhz": 1200.0, "h_carrier_ppm": 8.5,
                "waveform_json": "waveform.json",
            },
            "residues": {"A1": {"h_ppm_a": 6.25, "h_ppm_b": 6.583333333333333}},
        },
    }
    path = inputs / "config.json"
    path.write_text(json.dumps(config))
    return path


def check_sideband_benchmark():
    folder = Path(tempfile.mkdtemp(prefix="sbonest-benchmark-"))
    config = write_sideband_config(folder)
    out = run_cli(folder, config)
    assert out.returncode == 0, out.stdout + out.stderr
    assert "Fit completed" in out.stdout
    assert not (folder / "benchmark_profile.prof").exists()
    output = folder / "sideband.prof"
    out = run_cli(folder, config, "--profile-output", str(output))
    assert out.returncode == 0, out.stdout + out.stderr
    stats = pstats.Stats(str(output))
    assert any(Path(file).name == "sbfit.py" and name == "fit" for file, _, name in stats.stats)
    print("PASS: real Sideband benchmark, relative waveform and custom profile")


def check_onest_benchmark_and_profile_safety():
    folder = Path(tempfile.mkdtemp(prefix="onest-benchmark-"))
    config = folder / "config.json"
    config.write_text(json.dumps({
        "Project Name": "benchmark",
        "datasets": [str(ROOT / "example/syn10.txt")],
        "residues": [{"name": f"A{i}", "flag": "on" if i == 1 else "off"} for i in range(1, 6)],
        "init": {"Method": "Baldwin"},
    }))
    out = run_cli(folder, config, "profile")
    assert out.returncode == 0, out.stdout + out.stderr
    default = folder / "benchmark_profile.prof"
    stats = pstats.Stats(str(default))
    assert any(Path(file).name == "estmodel.py" and name == "fit" for file, _, name in stats.stats)
    for output, args in (
        (default, ("profile",)),
        (folder / "existing.prof", ("--profile-output", "existing.prof")),
    ):
        if not output.exists():
            output.write_bytes(b"preserved research artifact")
        before = output.read_bytes()
        out = run_cli(folder, config, *args)
        assert out.returncode != 0, "Existing profile was overwritten"
        assert output.read_bytes() == before
        assert "Loading config" not in out.stdout, "Profile collision checked after work started"
    print("PASS: real ONEST benchmark and early default/custom profile protection")


if __name__ == "__main__":
    check_sideband_benchmark()
    check_onest_benchmark_and_profile_safety()
