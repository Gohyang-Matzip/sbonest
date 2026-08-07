#!/usr/bin/env python

"""
Original by Donghan Lee
Python 3 port, scipy.optimize.least_squares integration, parallelization: 2025
"""

import sys
import numpy as np
from estmodel import est_model
from run import load_config, set_residue_flags
import multiprocessing
import os


def run_single_mc_iteration(args_tuple):
    """One MC iteration: reload data with noise, refit from the optimized params."""
    conf, pp_optimized_params_for_mc, worker_id = args_tuple
    mmc = est_model()
    mmc.verbose = False
    try:
        for dataset_name_mc in conf["datasets"]:
            mmc.dataset.addDataWithError(dataset_name_mc)
        if set_residue_flags(mmc.dataset, conf["residues"]) is not None:
            return None
        out_mc_fit = mmc.fit(p0=pp_optimized_params_for_mc, fitting_config=conf["init"])
        return out_mc_fit[0]
    except Exception:
        return None


def main():
    if len(sys.argv) < 3:
        sys.stderr.write(
            f"Usage: {sys.argv[0]} config_file.json number_of_mc_runs [num_processes]\n"
        )
        sys.stderr.write(
            "  [num_processes] is optional, defaults to number of CPU cores.\n"
        )
        sys.exit(1)

    try:
        nrun = int(sys.argv[2])
        if nrun <= 0:
            raise ValueError("Number of MC runs must be positive.")
    except ValueError as e:
        sys.stderr.write(f"Error: Invalid number of MC runs '{sys.argv[2]}'. {e}\n")
        sys.exit(1)

    num_processes = None
    if len(sys.argv) > 3:
        try:
            num_processes = int(sys.argv[3])
            if num_processes <= 0:
                raise ValueError("Number of processes must be positive.")
        except ValueError as e:
            sys.stderr.write(
                f"Error: Invalid num_processes '{sys.argv[3]}'. {e}\nUsing default.\n"
            )
            num_processes = None

    conf = load_config(sys.argv[1])

    m2 = est_model()
    m2.verbose = True
    print(
        "**************\n"
        + m2.programName
        + " - MC Run Setup (Parallelized)\n**************"
    )
    project_name = conf["Project Name"]
    for dataset_name in conf["datasets"]:
        m2.dataset.addData(dataset_name)
    err = set_residue_flags(m2.dataset, conf["residues"])
    if err:
        sys.stderr.write(f"Error: {err}\n")
        sys.exit(1)

    if m2.verbose:
        print(f"Generating data PDF: {project_name}_data.pdf")
    m2.datapdf(f"{project_name}_data.pdf")
    if m2.verbose:
        print("Performing initial fit...")
    out_initial_fit = m2.fit(fitting_config=conf["init"])
    pp_optimized_params = out_initial_fit[0]
    if m2.verbose:
        print(f"Generating PDF for initial fit: {project_name}.pdf")
    m2.pdf(pp_optimized_params, f"{project_name}.pdf")

    tasks_args = [(conf, pp_optimized_params, i) for i in range(nrun)]
    if num_processes is None:
        num_processes = os.cpu_count()
    print(
        f"\nStarting {nrun} Monte Carlo runs in parallel using {num_processes} CPU cores..."
    )

    with multiprocessing.Pool(processes=num_processes) as pool:
        results_from_pool = pool.map(run_single_mc_iteration, tasks_args)

    all_mc_fitted_params_list = [res for res in results_from_pool if res is not None]
    successful_runs = len(all_mc_fitted_params_list)
    failed_runs = nrun - successful_runs
    print(f"Completed {successful_runs}/{nrun} MC runs successfully.")
    if failed_runs > 0:
        print(f"Warning: {failed_runs} MC runs failed and were excluded.")
    if not all_mc_fitted_params_list:
        print("No MC runs successful. Exiting.")
        sys.exit(1)

    if m2.verbose:
        print("\nGenerating MC log buffer...")
    log_buffer_content_mc, mean_mc_params = m2.getLogBufferMC(
        pp_optimized_params, all_mc_fitted_params_list
    )

    mc_result_filename = f"{project_name}_mc.txt"
    if m2.verbose:
        print(f"Saving MC results to: {mc_result_filename}")
    with open(mc_result_filename, "w") as mc_file:
        mc_file.write(log_buffer_content_mc)

    if m2.verbose:
        print(f"Generating PDF for mean MC parameters: {project_name}_mcmean.pdf")
    m2.pdf(mean_mc_params, f"{project_name}_mcmean.pdf")

    if m2.verbose:
        std_dev_mc_params = np.std(np.array(all_mc_fitted_params_list), axis=0)
        print("\n--- Monte Carlo Parameter Statistics (from mcrun.py) ---")
        print(f"Number of successful MC runs: {successful_runs}")
        print("Mean Parameters (from MC) +/- StdDev (from MC):")
        for i in range(min(len(mean_mc_params), 15)):
            print(
                f"  Param {i}: {mean_mc_params[i]:.4f} +/- {std_dev_mc_params[i]:.4f}"
            )

    print("########\nMC Run script finished (Parallelized).\n")


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
