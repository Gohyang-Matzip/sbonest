#!/usr/bin/env python3
"""Run one side of a same-host archived/current numerical-core diagnostic.

Invoke each mode in a separate process with the same interpreter and thread
settings. This is not a reconstruction of the complete original environment:
the available run.py differs from its execution hash, and est_data.py was not
hashed in that record. No archived output or numerical baseline is rewritten.
"""

import argparse
import copy
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import types


ROOT = Path(__file__).resolve().parents[2]
REVISION = "4cdb65fdaea07a944b1aae4e53cb78b5966d6c3c"
ARCHIVE = Path("results/auto_H_refit/fits")
MODULES = ("fit", "est_data", "estmodel", "run", "sideband", "sbfit")
SENTINEL = "SBONEST_CORE_COMPARISON_JSON="


def _sha(data):
    return hashlib.sha256(data).hexdigest()


def _json_sha(value):
    return _sha(json.dumps(value, sort_keys=True, separators=(",", ":")).encode())


def _git(*args):
    return subprocess.check_output(["git", "-C", str(ROOT), *args])


def _sources(mode, metadata):
    sources, identities = {}, {}
    for name in MODULES:
        path = f"{name}.py"
        archive_path = str(ARCHIVE / "executed_source/sbfit.py") if name == "sbfit" else path
        archived = ((ROOT / archive_path).read_bytes() if name == "sbfit"
                    else _git("show", f"{REVISION}:{path}"))
        expected = metadata["execution_source_sha256"].get(archive_path)
        digest = _sha(archived)
        if name in ("fit", "estmodel", "sideband", "sbfit"):
            assert digest == expected, f"Archived execution hash mismatch: {path}"
        data = archived if mode == "archive" else (ROOT / path).read_bytes()
        origin = (f"{ROOT / archive_path}" if name == "sbfit" else f"git:{REVISION}:{path}")
        sources[name] = data
        identities[path] = {
            "loaded_from": origin if mode == "archive" else str(ROOT / path),
            "loaded_sha256": _sha(data),
            "archive_sha256": digest,
            "execution_sha256": expected,
            "archive_matches_execution": digest == expected if expected else None,
        }
    return sources, identities


def _load(mode, sources, identities):
    # Register each dependency before execution; no import may silently choose
    # an installed or current copy of an archived core dependency.
    assert not any(name in sys.modules for name in MODULES)
    if mode == "current":
        sys.path.insert(0, str(ROOT))
    else:
        sys.path[:] = [p for p in sys.path if Path(p or os.getcwd()).resolve() != ROOT]
    for name in MODULES:
        module = types.ModuleType(name)
        module.__file__ = identities[f"{name}.py"]["loaded_from"]
        sys.modules[name] = module
        exec(compile(sources[name], module.__file__, "exec"), module.__dict__)
    assert sys.modules["estmodel"].fit_module is sys.modules["fit"]
    assert sys.modules["estmodel"].EstDataSet is sys.modules["est_data"].EstDataSet
    assert sys.modules["sbfit"].est_model is sys.modules["estmodel"].est_model
    assert sys.modules["sbfit"].profile is sys.modules["sideband"].profile
    assert sys.modules["sbfit"].generate_initial_parameters is sys.modules["fit"].generate_initial_parameters
    assert sys.modules["sbfit"].set_residue_flags is sys.modules["run"].set_residue_flags
    return sys.modules["sbfit"]


def _config(out, metadata):
    original = json.loads((ROOT / ARCHIVE / "two_RF.json").read_text())
    historical = json.loads((ROOT / ARCHIVE / "two_RF_result.json").read_text())
    assert original == historical["config"], "Config differs from historical result"
    example_path = ROOT / "example/sideband_auto_H/two_RF.json"
    current = json.loads(example_path.read_text())
    path_keys = {"datasets", "Project Name"}
    assert {k: v for k, v in original.items() if k not in path_keys} == {
        k: v for k, v in current.items() if k not in path_keys
    }, "Current example differs from archived non-path configuration"
    config = copy.deepcopy(original)
    inputs = {}
    for index, filename in enumerate(original["datasets"]):
        relative = Path("results") / filename.split("/results/", 1)[1]
        digest = _sha((ROOT / relative).read_bytes())
        expected = metadata["execution_source_sha256"][str(relative)]
        assert digest == expected, f"Input hash mismatch: {relative}"
        current_input = (example_path.parent / current["datasets"][index]).resolve()
        assert current_input == (ROOT / relative).resolve(), "Current dataset path differs"
        assert _sha(current_input.read_bytes()) == digest, "Current dataset bytes differ"
        config["datasets"][index] = str(ROOT / relative)
        inputs[str(relative)] = digest
    config["Project Name"] = str(out.with_suffix(""))
    restored = copy.deepcopy(config)
    restored["datasets"] = original["datasets"]
    restored["Project Name"] = original["Project Name"]
    assert restored == original, "Rebasing changed non-path configuration"
    assert len(current["datasets"]) == len(original["datasets"])
    return config, historical, inputs, _json_sha(original)


def main():
    """Write an exclusive JSON diagnostic after a real, unmodified model fit."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", required=True, choices=("archive", "current"))
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    out = args.out.resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    # Fail before expensive work if this name has already been used.
    with out.open("x") as stream:
        metadata = json.loads((ROOT / ARCHIVE / "metadata.json").read_text())
        config, historical, inputs, config_sha = _config(out, metadata)
        sources, identities = _sources(args.mode, metadata)
        sbfit = _load(args.mode, sources, identities)
        import numpy as np
        import scipy

        optimalcontrol_version = importlib.metadata.version("optimalcontrol-nmr")
        model = sbfit.SidebandModel(config, ROOT / ARCHIVE)
        model.verbose = True
        names = model.parameter_names
        assert names == list(historical["parameters"]), "Parameter ordering differs"
        historical_vector = np.array([historical["parameters"][n]["value"] for n in names])
        fixed_residual = np.asarray(model.errFunc(historical_vector)).copy()
        observed = {}
        real_least_squares = sbfit.least_squares

        def observe_least_squares(fun, x0, **kwargs):
            # Forward every argument and return the actual SciPy result. Capture
            # the supplied x0, not a manually recreated or historical start.
            assert not observed, "Expected exactly one optimizer invocation"
            assert np.array_equal(model.free, np.arange(len(names)))
            observed["initial_parameters"] = np.asarray(x0).tolist()
            observed["optimizer_settings"] = {
                key: kwargs[key] for key in
                ("diff_step", "x_scale", "ftol", "xtol", "gtol", "max_nfev")
            }
            observed["jacobian_callable"] = callable(kwargs["jac"])
            print("SBONEST_CORE_OPTIMIZER_START=" + json.dumps(observed), flush=True)
            return real_least_squares(fun, x0, **kwargs)

        sbfit.least_squares = observe_least_squares
        started = time.monotonic()
        try:
            values, covariance = model.fit(fitting_config=config["init"])
        finally:
            sbfit.least_squares = real_least_squares
        elapsed = time.monotonic() - started
        assert observed, "Optimizer entry was not observed"
        result = model.result
        residual = np.asarray(model.errFunc(values))
        assert result.success, f"Optimizer failed: {result.message}"
        assert (residual.size, len(values), model.dof, model.rank) == (882, 24, 858, 24)
        np.testing.assert_array_equal(residual, result.fun)
        jacobian = np.ascontiguousarray(result.jac, dtype="<f8")
        loaded_local = {}
        for name, module in list(sys.modules.items()):
            filename = getattr(module, "__file__", None)
            if filename and Path(filename).is_absolute():
                path = Path(filename).resolve()
                if path.parent == ROOT and path.suffix == ".py":
                    loaded_local[name] = {"path": str(path), "sha256": _sha(path.read_bytes())}
        if args.mode == "archive":
            assert not loaded_local, f"Current root modules leaked into archive: {loaded_local}"
        report = {
            "mode": args.mode,
            "scope": "Same-host numerical-core comparison, not complete original reconstruction",
            "limitations": [
                "Available archived run.py differs from execution_source_sha256; its helpers are used, not its CLI.",
                "Archived est_data.py has no execution_source_sha256 record; its identity is not execution-verified.",
                "Both modes use present-host dependencies, not a recovered historical environment.",
            ],
            "archive_revision": REVISION,
            "checkout_revision": _git("rev-parse", "HEAD").decode().strip(),
            "source_identities": identities,
            "additional_loaded_current_sources": loaded_local,
            "inputs_sha256": inputs,
            "canonical_config_sha256": config_sha,
            "rebased_config": config,
            "config_semantically_equivalent": True,
            "runtime": {
                "python": sys.version, "executable": sys.executable,
                "platform": platform.platform(), "numpy": np.__version__,
                "scipy": scipy.__version__,
                "optimalcontrol-nmr": optimalcontrol_version,
                "threads": {key: os.environ.get(key) for key in
                            ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS")},
            },
            "parameter_names": names,
            **observed,
            "historical_vector": historical_vector.tolist(),
            "historical_stderr": [historical["parameters"][n]["stderr"] for n in names],
            "historical_chi2": historical["chi2"],
            "fixed_historical_vector_residual": fixed_residual.tolist(),
            "fixed_historical_vector_chi2": float(fixed_residual @ fixed_residual),
            "fitted_values": values.tolist(),
            "stderr": np.sqrt(np.diag(covariance)).tolist(),
            "covariance": covariance.tolist(),
            "residuals": residual.tolist(),
            "chi2": float(residual @ residual),
            "jacobian": jacobian.tolist(),
            "jacobian_sha256": _sha(jacobian.tobytes()),
            "jacobian_shape": list(jacobian.shape),
            "jacobian_encoding": "C-contiguous little-endian float64",
            "status": int(result.status), "success": bool(result.success),
            "n_points": int(residual.size), "n_parameters": len(values), "dof": int(model.dof),
            "message": result.message, "nfev": int(result.nfev), "njev": int(result.njev),
            "rank": int(model.rank), "scaled_condition": float(model.condition),
            "elapsed_seconds": elapsed,
        }
        text = json.dumps(report, allow_nan=False, sort_keys=True)
        stream.write(text + "\n")
        stream.flush()
        print(SENTINEL + text, flush=True)


if __name__ == "__main__":
    main()
