#!/usr/bin/env python

"""
Original written by Donghan Lee 2025
Python 3 / scipy.optimize.least_squares port: 2025
"""

import sys
import json
from estmodel import est_model


def load_config(config_file_path):
    """설정 파일을 로드하고 검증합니다."""
    try:
        with open(config_file_path, "r") as config_file:
            config = json.load(config_file)
            if not all(
                key in config
                for key in ["Project Name", "datasets", "residues", "init"]
            ):
                raise ValueError("Invalid config file structure.")
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
        name, flag = entry["name"], entry["flag"]
        res = next((r for r in dataset.res if r.label == name), None)
        if res is None:
            return f"Residue {name} not found in dataset."
        if flag not in ("on", "off"):
            return f"Wrong flag '{flag}' for residue {name}"
        res.active = flag == "on"
    return None


def main():
    if len(sys.argv) < 2:
        sys.stderr.write(f"Usage: {sys.argv[0]} config_file.json\n")
        sys.exit(1)

    config = load_config(sys.argv[1])

    model = est_model()
    model.verbose = True

    print("**************")
    print(model.programName)
    print("**************")

    for dataset_name in config["datasets"]:
        if model.verbose:
            print(f"Loading dataset: {dataset_name}")
        model.dataset.addData(dataset_name)

    err = set_residue_flags(model.dataset, config["residues"])
    if err:
        sys.stderr.write(f"Error: {err}\n")
        sys.exit(1)

    project_name = config["Project Name"]
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

    if model.verbose:
        print(f"Generating results PDF: {project_name}.pdf")
    model.pdf(optimized_params, f"{project_name}.pdf")

    print("########\nRun script finished.\n")


if __name__ == "__main__":
    main()
