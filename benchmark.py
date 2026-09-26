import sys
import time
import cProfile
import pstats
from estmodel import est_model
from run import load_config, load_datasets


def run_benchmark(config_path, profile=False):
    print(f"Loading config from {config_path}")
    config = load_config(config_path)

    model = est_model()
    model.verbose = True

    print("Loading datasets...")
    load_datasets(model, config)

    print("Starting fit...")
    start_time = time.time()

    if profile:
        profiler = cProfile.Profile()
        profiler.enable()

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
    do_profile = len(sys.argv) > 2 and sys.argv[2] == "profile"
    run_benchmark(config_file, do_profile)
