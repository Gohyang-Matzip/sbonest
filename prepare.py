#!/usr/bin/env python

import argparse
import json
import sys
from time import ctime
from pathlib import Path # Using pathlib for path manipulation
from estmodel import est_model # Assuming estmodel.py is correctly named and in path

def process_data(input_file_paths_str, method='Baldwin'): # Takes list of string paths
    """CEST/DEST 데이터를 처리하고 JSON 결과를 생성합니다."""
    m2 = est_model()
    # m2.verbose = True # Set if debugging est_model's internal prints is needed

    new_dict = {
        'Header': [m2.programName, f'Time: {ctime()}'], # Python 3 f-string
        'Project Name': 'default',
        'init': {
            'kex': {'min': 10.0, 'max': 400.0, 'nsteps': 6},
            'pB': {'min': 0.01, 'max': 0.1, 'nsteps': 6},
            'Method': method, # Use selected method
        },
        'datasets': input_file_paths_str, # Store as list of strings as per original
    }

    for file_path_str in input_file_paths_str:
        try:
            m2.dataset.addData(file_path_str) # addData expects a string filename
        except FileNotFoundError:
            print(f"Error: Input file '{file_path_str}' not found.", file=sys.stderr)
            sys.exit(1)
        except Exception as e: # Catch other potential errors from addData
            print(f"Error processing file '{file_path_str}': {e}", file=sys.stderr)
            sys.exit(1)

    new_dict['residues'] = m2.dataset.getResidues()
    return new_dict

def main():
    """명령행 인수를 처리하고 데이터를 처리합니다."""
    parser = argparse.ArgumentParser(description="Process CEST/DEST experiment data and generate JSON config.")
    parser.add_argument("input_files", nargs="+", type=str, help="Path(s) to the input data file(s)") # Keep as str for now
    parser.add_argument("--method", type=str, default="Baldwin", 
                        choices=['Baldwin', 'Matrix', 'NoEx', 'Matrix_3st_Linear', 'Matrix_3st_Triangle'],
                        help="Calculation method to use (default: Baldwin)")
    args = parser.parse_args()

    # input_files from args is already a list of strings
    result_dict = process_data(args.input_files, method=args.method)
    
    # Custom float formatting for json.dumps using the default argument
    # This replaces the old json.encoder.FLOAT_REPR method
    json_output = json.dumps(result_dict, indent=4, sort_keys=True, 
                             default=lambda o: format(o, '.4f') if isinstance(o, float) else o)
    print(json_output)

    print("Done preparing JSON configuration.", file=sys.stderr)

if __name__ == "__main__":
    main()
