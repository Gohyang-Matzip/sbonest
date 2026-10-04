"""Sideband run orchestration: preflight checks, checkpointed fits, analyses and exports.

The model lives in sbfit.py; this module owns the configuration-level workflow
(check_config, run_config) and the command line shared by run.py and sbfit.py.
"""
import argparse
import contextlib
import json
from pathlib import Path
import tempfile
import time

import numpy as np

from run import load_config
from sb_analysis import (
    MultiStartError, _json_safe, fit_multistart, profile_intervals, profile_likelihood,
    restore_fit, snapshot_fit,
)
from sb_report import fit_pdf, prediction_rows, provenance, write_predictions
from sbfit import SidebandModel, _validate_fit_config


def _output_paths(config, no_pdf=False):
    project = Path(config["Project Name"])
    outputs = [
        Path(str(project) + suffix)
        for suffix in ("_result.txt", "_result.json", "_predictions.csv")
    ]
    if not no_pdf:
        outputs += [Path(str(project) + ".pdf"), Path(str(project) + "_data.pdf")]
    return outputs


def check_config(config, config_dir=".", *, no_pdf=False, identifiability=False, workers=1):
    """Return a JSON-safe preflight summary without optimization or file writes.

    The existing initialization grid is evaluated to check the actual baseline.
    Invalid settings and output conflicts are reported in ``errors``. An invalid
    baseline is reported even when later multistart attempts could recover it.
    With ``identifiability``, one grouped Jacobian at the initial point adds local
    expected errors, rank, condition, weak parameters and strong correlations.
    """
    from sb_parallel import WorkerPool, validate_workers

    workers = validate_workers(workers)
    summary = {"valid": False, "method": "Sideband", "errors": [], "warnings": [],
               "datasets": [], "waveforms": [], "outputs": [], "output_conflicts": []}
    try:
        if (not isinstance(config, dict)
                or not isinstance(config.get("Project Name"), str)
                or not config["Project Name"].strip()
                or not isinstance(config.get("datasets"), list) or not config["datasets"]
                or not all(isinstance(path, str) and path for path in config["datasets"])
                or not isinstance(config.get("residues"), list)
                or not isinstance(config.get("init"), dict)):
            raise ValueError("Invalid config value types or missing required sections")
        outputs = _output_paths(config, no_pdf)
        outputs.append(Path(config["Project Name"] + "_checkpoint"))
        summary["project"] = str(Path(config["Project Name"]).absolute())
        summary["outputs"] = [str(path.absolute()) for path in outputs]
        summary["output_conflicts"] = [str(path.absolute()) for path in outputs
                                       if path.exists() or path.is_symlink()]
        if summary["output_conflicts"]:
            summary["errors"].append("Output already exists; choose a new Project Name")
        for parent in outputs[0].parents:
            if parent.exists() and not parent.is_dir():
                summary["errors"].append(f"Output parent is not a directory: {parent.absolute()}")
                break
        model = SidebandModel(config, config_dir)
        _validate_fit_config(model, config["init"])
        summary["model"] = model.sideband_method
        summary["states"] = model.states
        summary["field_groups"] = [{"index": g, "h_larmor_mhz": larmor,
                                    "datasets": [i for i, gi in enumerate(model.dataset_group) if gi == g]}
                                   for g, larmor in enumerate(model.field_groups)]
        data = model._prepare_data()
        for i, path in enumerate(config["datasets"]):
            spectra = [item for item in data if item[1].dataset_index == i]
            summary["datasets"].append({
                "index": i, "path": str((Path(config_dir) / path).resolve()),
                "n_points": sum(len(item[2]) for item in spectra),
                "n_spectra": len(spectra),
                "residues": [model.dataset.res[item[0]].label for item in spectra],
                "nitrogen_fields_mhz": sorted({float(item[1].field) for item in spectra}),
                "proton_field_mhz": float(model.decoupling[i]["h_larmor_mhz"]),
                "nominal_v1n_hz": float(model.dataset.v1s[i]),
            })
        summary["waveforms"] = sorted({str((Path(config_dir) / d["waveform_json"]).resolve())
                                       for d in model.decoupling if "waveform_json" in d})
        vary = model.init_config().get("vary", model.parameter_names)
        summary.update(
            n_points=sum(len(item[2]) for item in data),
            n_residues=sum(r.active for r in model.dataset.res),
            n_parameters=len(vary), parameter_order=list(model.parameter_names),
            free_parameters=list(vary),
            fixed_parameters=[name for name in model.parameter_names if name not in vary],
            rf_mode=model.rf_mode, proton_relaxation_mode=model.proton_mode,
            analyses=[name for name in ("multistart", "profile", "profile_interval", "bootstrap")
                      if name in config["init"]],
        )
        initial = model.prepare_fit()
        summary["initial_parameters"] = dict(zip(model.parameter_names, initial.tolist()))
        summary["bounds"] = {
            name: [float(value) if np.isfinite(value) else None
                   for value in (model.lower[i], model.upper[i])]
            for i, name in enumerate(model.parameter_names)
        }
        summary["bounds_note"] = "null bounds denote unbounded directions"
        if identifiability:
            from sb_analysis import identifiability as _identifiability

            pool_context = (WorkerPool(workers, config, config_dir) if workers > 1
                            else contextlib.nullcontext())
            with pool_context as pool:
                summary["identifiability"] = _identifiability(model, initial, pool=pool)
            report = summary["identifiability"]
            if report["jacobian_rank"] < report["n_free"]:
                summary["warnings"].append(
                    "Rank-deficient Jacobian at the initial point: parameters are not separately identifiable there.")
            if report["weak_parameters"]:
                summary["warnings"].append(
                    "Weakly determined at the initial point: " + ", ".join(report["weak_parameters"]))
            if report["strong_correlations"]:
                summary["warnings"].append(
                    "Strongly correlated at the initial point: " + ", ".join(
                        "/".join(row["parameters"]) for row in report["strong_correlations"][:5]))
        summary["warnings"].append(
            "Preflight validates configuration and initialization; convergence and identifiability require fitting."
        )
    except (ValueError, KeyError, TypeError, OSError, RuntimeError, IndexError) as exc:
        summary["errors"].append(f"{type(exc).__name__}: {exc}")
    summary["valid"] = not summary["errors"]
    return summary


def run_config(config, config_dir=".", no_pdf=False, *, resume=False, workers=1):
    """Fit a Sideband configuration with checkpoints, optional analyses and exclusive exports; returns (parameters, covariance)."""
    from sb_bootstrap import bootstrap_fit
    from sb_checkpoint import Checkpoint, execution_identity, output_plan, publish_outputs
    from sb_parallel import WorkerPool, validate_workers

    workers = validate_workers(workers)
    model = SidebandModel(config, config_dir)
    _validate_fit_config(model, config["init"])
    project = Path(config["Project Name"])
    outputs = _output_paths(config, no_pdf)
    if not resume and any(path.exists() or path.is_symlink() for path in outputs):
        raise FileExistsError("Output already exists; choose a new Project Name to preserve previous fits")
    metadata = provenance(config, config_dir)
    began = time.monotonic()
    checkpoint = Path(str(project) + "_checkpoint")
    # Worker count never enters the checkpoint identity: results do not depend on it.
    pool_context = WorkerPool(workers, config, config_dir) if workers > 1 else contextlib.nullcontext()
    with Checkpoint(checkpoint, execution_identity(metadata, no_pdf=no_pdf), resume=resume) as journal, \
            pool_context as pool:
        model.pool = pool
        if not resume:
            journal.save("provenance", metadata)
        metadata = journal.read("provenance")
        if not isinstance(metadata, dict):
            raise ValueError("Missing checkpoint provenance")

        def completed(prefix):
            rows = []
            while (row := journal.read(f"{prefix}-{len(rows)}")) is not None:
                rows.append(row)
            # A hole is corruption, never permission to silently redo/drop work.
            if len(list(checkpoint.glob(prefix + "-*.json"))) != len(rows):
                raise ValueError(f"Checkpoint records are not an ordered prefix: {prefix}")
            return rows

        def publish_failure():
            plan = journal.read("failure-exports")
            if plan is None:
                stage = Path(tempfile.mkdtemp(prefix="failure-export-", dir=checkpoint)) / outputs[1].name
                stage.write_text(json.dumps(journal.read("failure"), indent=2, allow_nan=False) + "\n",
                                 encoding="utf-8")
                plan = output_plan(journal, [stage], [outputs[1]])
                journal.save("failure-exports", plan)
            publish_outputs(journal, plan, [outputs[1]])

        def save_failure(exc):
            info = {"schema_version": 2, "success": False, "method": "Sideband",
                    "config": config, "parameter_order": model.parameter_names,
                    "message": f"{type(exc).__name__}: {exc}",
                    "provenance": {**metadata, "elapsed_s": time.monotonic() - began}}
            if isinstance(exc, MultiStartError):
                info["multistart"] = exc.attempts
            if model.result is not None:
                info["optimizer"] = _json_safe(dict(model.result))
            record = journal.read("failure")
            if record is None:
                journal.save("failure", info)
            publish_failure()

        model.verbose = True
        baseline = journal.read("baseline")
        prior_failure = journal.read("failure")
        if baseline is None and prior_failure is not None:
            publish_failure()
            raise RuntimeError("This checkpoint records a completed failed fit; inspect its result and use a fresh Project Name after correcting settings")
        if baseline is None:
            try:
                attempts = None
                if "multistart" in config["init"]:
                    p, covariance, attempts = fit_multistart(
                        model, config["init"]["multistart"], completed=completed("attempt"),
                        on_complete=lambda row: journal.save(f"attempt-{row['index']}", row),
                        pool=pool)
                    fitted = p, covariance
                else:
                    fitted = model.fit(fitting_config=config["init"])
            except (ValueError, RuntimeError) as exc:
                save_failure(exc)
                raise
            baseline = {"snapshot": snapshot_fit(model, *fitted), "multistart": attempts}
            journal.save("baseline", baseline)
        else:
            fitted = restore_fit(model, baseline["snapshot"])
        if journal.read("exports") is not None:
            publish_outputs(journal, journal.read("exports"), outputs)
            journal.save("complete", {"success": True})
            print(f"Resumed completed outputs for {project}")
            return fitted
        info = model.diagnostics(*fitted)
        info["success"] = True
        rows = prediction_rows(model, fitted[0])
        from sb_diagnostics import residual_diagnostics

        info["residual_diagnostics"] = residual_diagnostics(rows, info["n_parameters"], fitted[1])
        for name, value in zip(model.parameter_names, info["residual_diagnostics"]["rescaled_stderr"]):
            info["parameters"][name]["stderr_rescaled"] = value if info["parameters"][name]["vary"] else None
        info["warnings"].extend(info["residual_diagnostics"]["warnings"])
        if baseline["multistart"] is not None:
            info["multistart"] = [
                {key: value for key, value in row.items() if key not in ("initialization", "fit_snapshot")}
                for row in baseline["multistart"]
            ]
        if "profile" in config["init"]:
            existing = {name: completed(f"profile-{name}") for name in config["init"]["profile"]}
            counts = {name: len(rows) for name, rows in existing.items()}

            def save_point(name, row):
                journal.save(f"profile-{name}-{counts[name]}", row)
                counts[name] += 1

            info["profiles"] = profile_likelihood(model, fitted[0], config["init"]["profile"],
                                                  completed=existing, on_complete=save_point,
                                                  pool=pool)
            info["warnings"].extend(info["profiles"]["warnings"])
        if "profile_interval" in config["init"]:
            names = config["init"]["profile_interval"]["parameters"]
            existing = {name: completed(f"profile_interval-{name}") for name in names}
            counts = {name: len(rows) for name, rows in existing.items()}

            def save_interval_point(name, row):
                journal.save(f"profile_interval-{name}-{counts[name]}", row)
                counts[name] += 1

            info["profile_intervals"] = profile_intervals(
                model, fitted[0], fitted[1], config["init"]["profile_interval"],
                completed=existing, on_complete=save_interval_point, pool=pool)
            info["warnings"].extend(info["profile_intervals"]["warnings"])
        if "bootstrap" in config["init"]:
            info["bootstrap"] = bootstrap_fit(
                model, fitted[0], config["init"]["bootstrap"], completed=completed("bootstrap"),
                on_complete=lambda row: journal.save(f"bootstrap-{row['index']}", row),
                pool=pool)
            info["warnings"].extend(info["bootstrap"]["warnings"])
        metadata = {**metadata, "elapsed_s": time.monotonic() - began,
                    "elapsed_note": "Seconds in the final invocation; excludes earlier interrupted invocations.",
                    "resumed": bool(resume), "workers": workers}
        info["provenance"] = metadata
        info["checkpoint"] = str(checkpoint.absolute())
        stage = Path(tempfile.mkdtemp(prefix="export-", dir=checkpoint))
        staged = [stage / path.name for path in outputs]
        staged[0].write_text(model.getLogBuffer(fitted), encoding="utf-8")
        staged[1].write_text(json.dumps(info, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        write_predictions(staged[2], rows)
        if not no_pdf:
            fit_pdf(staged[3], rows)
            model.datapdf(str(staged[4]))
        plan = output_plan(journal, staged, outputs)
        journal.save("exports", plan)
        publish_outputs(journal, plan, outputs)
        journal.save("complete", {"success": True})
    for warning in info["warnings"]:
        print(f"Warning: {warning}")
    print(f"Saved {outputs[0]} and {outputs[1]}; checkpoint: {checkpoint}")
    return fitted


def main():
    """Command line shared by run.py and sbfit.py for Sideband configurations."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config_file")
    parser.add_argument("--no-pdf", action="store_true")
    parser.add_argument("--check", action="store_true", help="Validate without fitting or writing outputs")
    parser.add_argument("--identifiability", action="store_true",
                        help="With --check: add local Jacobian diagnostics at the initial point")
    parser.add_argument("--resume", action="store_true", help="Resume a matching Sideband checkpoint")
    parser.add_argument("--workers", type=int, default=1,
                        help="Worker processes for Jacobian columns, restarts, profile points and bootstrap replicates (default 1)")
    args = parser.parse_args()
    if args.check and args.resume:
        parser.error("--check and --resume cannot be combined")
    if args.identifiability and not args.check:
        parser.error("--identifiability requires --check")
    if args.workers < 1:
        parser.error("--workers must be a positive integer")
    try:
        config = load_config(args.config_file)
        config_dir = Path(args.config_file).resolve().parent
        if args.check:
            summary = check_config(config, config_dir, no_pdf=args.no_pdf,
                                   identifiability=args.identifiability, workers=args.workers)
            print(json.dumps(summary, indent=2, allow_nan=False))
            parser.exit(0 if summary["valid"] else 1)
        run_config(
            config,
            config_dir,
            args.no_pdf,
            resume=args.resume,
            workers=args.workers,
        )
    except (ValueError, KeyError, OSError, RuntimeError) as exc:
        parser.exit(1, f"Error: {exc}\n")


if __name__ == "__main__":
    main()
