import sys
import os
import time
import json
import cProfile
import pstats
from estmodel import est_model

# Add current directory to sys.path
sys.path.append(os.getcwd())


def run_benchmark(config_path, profile=False):
    print(f"Loading config from {config_path}")
    with open(config_path, "r") as f:
        config = json.load(f)

    # Adjust dataset paths to be relative to the config file or absolute
    config_dir = os.path.dirname(config_path)
    new_datasets = []
    for ds in config["datasets"]:
        new_datasets.append(os.path.join(config_dir, ds))
    config["datasets"] = new_datasets

    model = est_model()
    model.verbose = True  # Keep verbose to see progress

    print("Loading datasets...")
    for dataset_name in config["datasets"]:
        model.dataset.addData(dataset_name)

    from run import set_residue_flags

    set_residue_flags(model.dataset, config["residues"])

    print("Starting fit...")
    start_time = time.time()

    if profile:
        profiler = cProfile.Profile()
        profiler.enable()

    # Run fit
    # We use the config['init'] for fitting config
    model.fit(fitting_config=config["init"])

    if profile:
        profiler.disable()
        stats = pstats.Stats(profiler).sort_stats("cumtime")
        stats.print_stats(20)
        stats.dump_stats("benchmark_profile.prof")

    end_time = time.time()
    print(f"Fit completed in {end_time - start_time:.4f} seconds.")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python benchmark.py <config_file> [profile]")
        sys.exit(1)

    config_file = sys.argv[1]
    do_profile = False
    if len(sys.argv) > 2 and sys.argv[2] == "profile":
        do_profile = True

    run_benchmark(config_file, do_profile)
