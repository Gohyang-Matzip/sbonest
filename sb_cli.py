"""Unified SBONEST command line: ``sbonest <command> ...``.

Every command delegates to the existing module functions, so behaviour,
outputs and provenance are identical to running the scripts directly.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

COMMANDS = ("check", "fit", "resume", "report", "init-demo", "design", "compare",
            "import-bruker", "serve", "benchmark", "version")


def _config(path):
    from run import load_config

    return load_config(path), Path(path).resolve().parent


def _fit_arguments(parser):
    parser.add_argument("config_file")
    parser.add_argument("--no-pdf", action="store_true", help="Skip PDF figures")
    parser.add_argument("--workers", type=int, default=1, help="Worker processes (default 1)")


def build_parser():
    """Argument parser of the sbonest command."""
    parser = argparse.ArgumentParser(prog="sbonest", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)

    check = commands.add_parser("check", help="Validate a configuration without fitting or writing outputs")
    _fit_arguments(check)
    check.add_argument("--identifiability", action="store_true", help="Add local Jacobian diagnostics")

    fit = commands.add_parser("fit", help="Fit a Sideband configuration")
    _fit_arguments(fit)
    resume = commands.add_parser("resume", help="Resume a matching checkpoint")
    _fit_arguments(resume)

    report = commands.add_parser("report", help="Report a saved result without refitting")
    report.add_argument("result_json")
    report.add_argument("--out", required=True)
    report.add_argument("--predictions")

    demo = commands.add_parser("init-demo", help="Copy the bundled synthetic example into a new folder")
    demo.add_argument("--out", required=True)

    design = commands.add_parser("design", help="Expected errors of planned acquisitions at a known truth")
    design.add_argument("design_json")
    design.add_argument("--out", required=True)
    design.add_argument("--workers", type=int, default=1)

    compare = commands.add_parser("compare", help="Shared-exchange versus per-residue fits")
    compare.add_argument("config_file")
    compare.add_argument("--out", required=True)
    compare.add_argument("--workers", type=int, default=1)
    compare.add_argument("--pdf", action="store_true")
    compare.add_argument("--models", nargs="+", metavar="MODEL",
                         help="Compare Sideband models on the same data instead of residues")
    compare.add_argument("--h-ppm-c", nargs="+", metavar="LABEL=ppm",
                         help="State-C proton shifts (default: the state-B shift)")

    commands.add_parser("import-bruker", help="Convert Bruker pseudo-2D data (see sb_import.py --help)",
                        add_help=False)
    serve = commands.add_parser("serve", help="Start the Sideband web runner")
    serve.add_argument("--host", default="127.0.0.1")
    serve.add_argument("--port", type=int, default=5050)
    serve.add_argument("--token", help="Access token (default: SBONEST_TOKEN environment variable)")
    serve.add_argument("--max-age-days", type=float, help="Archive idle jobs older than this many days")

    commands.add_parser("benchmark", help="Time or profile one fit (see benchmark.py --help)", add_help=False)
    commands.add_parser("version", help="Print the package version and executed-source hashes")
    return parser


def main(argv=None):
    """Entry point of the sbonest command; returns the exit status."""
    argv = list(sys.argv[1:] if argv is None else argv)
    passthrough = {"import-bruker": "sb_import", "benchmark": "benchmark"}
    if argv and argv[0] in passthrough:
        module = __import__(passthrough[argv[0]])
        return module.main(argv[1:])
    parser = build_parser()
    args = parser.parse_args(argv)
    for name in ("workers",):
        if getattr(args, name, 1) < 1:
            parser.error(f"--{name} must be a positive integer")
    try:
        if args.command in ("check", "fit", "resume"):
            from sbfit import check_config, run_config

            config, config_dir = _config(args.config_file)
            if args.command == "check":
                summary = check_config(config, config_dir, no_pdf=args.no_pdf,
                                       identifiability=args.identifiability, workers=args.workers)
                print(json.dumps(summary, indent=2, allow_nan=False))
                return 0 if summary["valid"] else 1
            run_config(config, config_dir, args.no_pdf, resume=args.command == "resume", workers=args.workers)
            return 0
        if args.command == "report":
            from sb_report import regenerate_report

            paths = regenerate_report(args.result_json, args.out, args.predictions)
        elif args.command == "init-demo":
            from sb_workflow import init_demo

            path = init_demo(args.out)
            paths = {"config": str(path)}
        elif args.command == "design":
            from sb_design import run_design

            paths = run_design(args.design_json, args.out, workers=args.workers)
        elif args.command == "compare":
            from sb_compare import compare_models, run_comparison

            config, config_dir = _config(args.config_file)
            if not str(config['init'].get('Method', '')).startswith('Sideband'):
                raise ValueError('compare requires a Sideband init.Method')
            if args.models:
                shifts = None
                if args.h_ppm_c:
                    shifts = {}
                    for item in args.h_ppm_c:
                        label, _, value = item.partition('=')
                        if not label or not value:
                            raise ValueError('--h-ppm-c entries must be LABEL=ppm')
                        shifts[label] = float(value)
                paths = compare_models(config, config_dir, args.out, models=tuple(args.models),
                                       h_ppm_c=shifts, no_pdf=not args.pdf, workers=args.workers)
            else:
                paths = run_comparison(config, config_dir, args.out, no_pdf=not args.pdf, workers=args.workers)
        elif args.command == "serve":
            import sb_server

            server_args = ["--host", args.host, "--port", str(args.port)]
            if args.token is not None:
                server_args.extend(["--token", args.token])
            if args.max_age_days is not None:
                server_args.extend(["--max-age-days", str(args.max_age_days)])
            return sb_server.main(server_args)
        else:
            from importlib.metadata import PackageNotFoundError, version

            try:
                print(f"sbonest {version('sbonest')}")
            except PackageNotFoundError:
                print("sbonest (not installed as a package; running from source)")
            root = Path(__file__).resolve().parent
            sources = ('run.py', 'sbfit.py', 'sideband.py', 'fit.py', 'estmodel.py',
                       'est_data.py', 'sb_report.py', 'sb_analysis.py', 'sb_workflow.py',
                       'sb_checkpoint.py', 'sb_bootstrap.py', 'sb_run.py', 'sb_diagnostics.py')
            for name in sources:
                print(f"{name}: {hashlib.sha256((root / name).read_bytes()).hexdigest()}")
            return 0
        for label, path in paths.items():
            print(f"{label}: {path}")
        return 0
    except (ValueError, KeyError, OSError, RuntimeError) as exc:
        parser.exit(1, f"Error: {exc}\n")


if __name__ == "__main__":
    sys.exit(main())
