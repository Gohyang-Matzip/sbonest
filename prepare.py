#!/usr/bin/env python

import argparse
import json
import sys
from time import ctime
from pathlib import Path
from estmodel import est_model, VALID_METHODS


def process_data(input_file_paths_str, method="Baldwin"):
    """CEST/DEST 데이터를 처리하고 JSON 결과를 생성합니다."""
    m2 = est_model()

    new_dict = {
        "Header": [m2.programName, f"Time: {ctime()}"],
        "Project Name": "default",
        "init": {
            "kex": {"min": 10.0, "max": 400.0, "nsteps": 6},
            "pB": {"min": 0.01, "max": 0.1, "nsteps": 6},
            "Method": method,
        },
        "datasets": [str(Path(p).resolve()) for p in input_file_paths_str],
    }

    for file_path_str in input_file_paths_str:
        try:
            m2.dataset.addData(file_path_str)
        except FileNotFoundError:
            print(f"Error: Input file '{file_path_str}' not found.", file=sys.stderr)
            sys.exit(1)
        except Exception as e:
            print(f"Error processing file '{file_path_str}': {e}", file=sys.stderr)
            sys.exit(1)

    new_dict["residues"] = m2.dataset.getResidues()
    return new_dict


def main():
    """명령행 인수를 처리하고 데이터를 처리합니다."""
    parser = argparse.ArgumentParser(
        description="Process CEST/DEST experiment data and generate JSON config."
    )
    parser.add_argument(
        "input_files", nargs="+", help="Path(s) to the input data file(s)"
    )
    parser.add_argument(
        "--method",
        type=str,
        default="Baldwin",
        choices=VALID_METHODS,
        help="Calculation method to use (default: Baldwin)",
    )
    args = parser.parse_args()

    result_dict = process_data(args.input_files, method=args.method)
    print(json.dumps(result_dict, indent=4, sort_keys=True))

    print("Done preparing JSON configuration.", file=sys.stderr)


if __name__ == "__main__":
    main()
