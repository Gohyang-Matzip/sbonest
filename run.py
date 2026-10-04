#!/usr/bin/env python

"""
Original written by Donghan Lee 2025
Python 3 / scipy.optimize.least_squares port: 2025
"""

import sys
import json
import argparse
from pathlib import Path
from estmodel import est_model


def load_config(config_file_path):
    """설정 파일을 로드하고 검증합니다."""
    try:
        with open(config_file_path, "r") as config_file:
            config = json.load(config_file)
            if not isinstance(config, dict) or not all(
                key in config
                for key in ["Project Name", "datasets", "residues", "init"]
            ):
                raise ValueError("Invalid config file structure.")
            if (
                not isinstance(config["Project Name"], str)
                or not config["Project Name"].strip()
                or not isinstance(config["datasets"], list)
                or not config["datasets"]
                or not all(isinstance(p, str) and p for p in config["datasets"])
                or not isinstance(config["residues"], list)
                or not isinstance(config["init"], dict)
            ):
                raise ValueError("Invalid config value types or empty dataset list.")
            if ("sideband" in config) != (config["init"].get("Method") == "Sideband"):
                raise ValueError(
                    'Use init.Method = "Sideband" together with a sideband section.'
                )
            config_dir = Path(config_file_path).resolve().parent
            config["datasets"] = [
                str((config_dir / p).resolve()) for p in config["datasets"]
            ]
            return config
    except FileNotFoundError:
        sys.stderr.write(f"Error: Config file not found: {config_file_path}\n")
        sys.exit(1)
    except json.JSONDecodeError:
        sys.stderr.write(
            f"Error: Invalid JSON format in config file: {config_file_path}\n"
        )
        sys.exit(1)
    except ValueError as e:
        sys.stderr.write(f"Error: {e}\n")
        sys.exit(1)


def set_residue_flags(dataset, residues_config):
    """config의 on/off 플래그를 데이터셋 잔기에 적용. 성공 시 None, 실패 시 오류 메시지 반환."""
    for entry in residues_config:
        if not isinstance(entry, dict) or not isinstance(entry.get("name"), str):
            return "Each residue must have a name and an on/off flag."
        name, flag = entry["name"], entry.get("flag")
        res = next((r for r in dataset.res if r.label == name), None)
        if res is None:
            return f"Residue {name} not found in dataset."
        if flag not in ("on", "off"):
            return f"Wrong flag '{flag}' for residue {name}"
        res.active = flag == "on"
    return None


def load_datasets(model, config, *, with_error=False):
    """Load configured spectra and apply residue selection for CLI and MC fits."""
    add_data = model.dataset.addDataWithError if with_error else model.dataset.addData
    for filename in config["datasets"]:
        if model.verbose:
            print(f"Loading dataset: {filename}")
        add_data(filename)
    error = set_residue_flags(model.dataset, config["residues"])
    if error:
        raise ValueError(error)


def main():
    parser = argparse.ArgumentParser(description="Fit CEST data.")
    parser.add_argument("config_file")
    parser.add_argument(
        "--no-pdf",
        action="store_true",
        help="Skip PDF reports; save numeric results only.",
    )
    parser.add_argument("--check", action="store_true", help="Validate Sideband inputs without fitting or writing outputs")
    parser.add_argument("--resume", action="store_true", help="Resume a matching Sideband checkpoint")
    args = parser.parse_args()
    if args.check and args.resume:
        parser.error("--check and --resume cannot be combined")
    config = load_config(args.config_file)

    if config["init"].get("Method") == "Sideband":
        from sbfit import check_config, run_config

        try:
            if args.check:
                summary = check_config(config, Path(args.config_file).resolve().parent,
                                       no_pdf=args.no_pdf)
                print(json.dumps(summary, indent=2, allow_nan=False))
                parser.exit(0 if summary["valid"] else 1)
            run_config(config, Path(args.config_file).resolve().parent, args.no_pdf, resume=args.resume)
        except (ValueError, KeyError, OSError, RuntimeError) as exc:
            parser.exit(1, f"Error: {exc}\n")
        return

    if args.check or args.resume:
        parser.exit(1, "Error: --check/--resume currently support init.Method = Sideband\n")

    model = est_model()
    model.verbose = True

    print("**************")
    print(model.programName)
    print("**************")

    try:
        load_datasets(model, config)
    except ValueError as exc:
        sys.stderr.write(f"Error: {exc}\n")
        sys.exit(1)

    project_name = config["Project Name"]
    if not args.no_pdf:
        if model.verbose:
            print(f"Generating data PDF: {project_name}_data.pdf")
        model.datapdf(f"{project_name}_data.pdf")

    if model.verbose:
        print("Starting model fitting...")
    fit_result_tuple = model.fit(fitting_config=config["init"])
    optimized_params = fit_result_tuple[0]

    if model.verbose:
        print("Generating log buffer...")
    log_buffer_content = model.getLogBuffer(fit_result_tuple)

    result_file_name = f"{project_name}_result.txt"
    if model.verbose:
        print(f"Saving results to: {result_file_name}")
    with open(result_file_name, "w") as result_file:
        result_file.write(log_buffer_content)

    if not args.no_pdf:
        if model.verbose:
            print(f"Generating results PDF: {project_name}.pdf")
        model.pdf(optimized_params, f"{project_name}.pdf")

    print("########\nRun script finished.\n")


if __name__ == "__main__":
    main()
