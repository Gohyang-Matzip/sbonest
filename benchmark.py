import argparse
import time
import cProfile
import marshal
from pathlib import Path
import pstats
from estmodel import est_model
from run import load_config, load_datasets


def run_benchmark(config_path, profile=False, profile_output=None):
    profile = profile or profile_output is not None
    profile_path = Path(profile_output or "benchmark_profile.prof")
    if profile and (profile_path.exists() or profile_path.is_symlink()):
        raise FileExistsError(f"Profile already exists; choose a fresh --profile-output: {profile_path}")

    print(f"Loading config from {config_path}")
    config = load_config(config_path)

    print("Loading datasets...")
    if str(config["init"].get("Method", "")).startswith("Sideband"):
        from sbfit import SidebandModel

        model = SidebandModel(config, Path(config_path).resolve().parent)
    else:
        model = est_model()
        model.verbose = True
        load_datasets(model, config)
    model.verbose = True

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
        # Exclusive creation also protects profiles created while this fit ran.
        with profile_path.open("xb") as profile_file:
            marshal.dump(stats.stats, profile_file)
        print(f"Profile saved to {profile_path}")

    end_time = time.time()
    print(f"Fit completed in {end_time - start_time:.4f} seconds.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Benchmark a Sideband or ONEST fit.")
    parser.add_argument("config_file")
    parser.add_argument("mode", nargs="?", choices=["profile"])
    parser.add_argument(
        "--profile-output",
        help="Enable profiling and save to this new file (default: benchmark_profile.prof).",
    )
    args = parser.parse_args()
    try:
        run_benchmark(args.config_file, args.mode == "profile", args.profile_output)
    except (ValueError, KeyError, OSError, RuntimeError) as exc:
        parser.exit(1, f"Error: {exc}\n")
