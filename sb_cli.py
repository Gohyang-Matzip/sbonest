"""Unified SBONEST command line: ``sbonest <command> ...``.

Every command delegates to the existing module command lines, so behaviour,
outputs and provenance are identical to running the scripts directly.
``check``, ``fit`` and ``resume`` accept Sideband configurations only; the
inherited ONEST models run through ``run.py``.
"""
import argparse
import hashlib
from pathlib import Path
import sys

COMMANDS = ("check", "fit", "resume", "report", "init-demo", "design", "compare",
            "import-bruker", "serve", "benchmark", "version")
# Commands whose arguments are handed to a module command line unchanged
# (``sbonest check CONFIG`` becomes ``sb_run CONFIG --check``).
_DELEGATED = {"check": ("sb_run", ["--check"]), "fit": ("sb_run", []), "resume": ("sb_run", ["--resume"]),
              "report": ("sb_workflow", ["report"]), "init-demo": ("sb_workflow", ["init-demo"]),
              "design": ("sb_workflow", ["design"]), "compare": ("sb_workflow", ["compare"]),
              "import-bruker": ("sb_import", []), "benchmark": ("benchmark", [])}
_HELP = {"check": "Validate a configuration without fitting or writing outputs (see run.py --help)",
         "fit": "Fit a Sideband configuration (see run.py --help)",
         "resume": "Resume a matching checkpoint (see run.py --help)",
         "report": "Report a saved result without refitting (see sb_workflow.py report --help)",
         "init-demo": "Copy the bundled synthetic example into a new folder",
         "design": "Expected errors of planned acquisitions at a known truth (see sb_workflow.py design --help)",
         "compare": "Shared-exchange versus per-residue or two- versus three-state fits (see sb_workflow.py compare --help)",
         "import-bruker": "Convert Bruker pseudo-2D data (see sb_import.py --help)",
         "benchmark": "Time or profile one fit (see benchmark.py --help)"}


def build_parser():
    """Argument parser of the sbonest command (delegated commands keep their module's help)."""
    parser = argparse.ArgumentParser(prog="sbonest", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in COMMANDS:
        if name in _DELEGATED:
            commands.add_parser(name, help=_HELP[name], add_help=False)
    serve = commands.add_parser("serve", help="Start the Sideband web runner")
    serve.add_argument("--host", default="127.0.0.1")
    serve.add_argument("--port", type=int, default=5050)
    serve.add_argument("--token", help="Access token (default: SBONEST_TOKEN environment variable)")
    serve.add_argument("--max-age-days", type=float, help="Archive idle jobs older than this many days")
    commands.add_parser("version", help="Print the package version and executed-source hashes")
    return parser


def main(argv=None):
    """Entry point of the sbonest command; returns the exit status."""
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] in _DELEGATED:
        module_name, prefix = _DELEGATED[argv[0]]
        module = __import__(module_name)
        if module_name == "sb_run":
            return module.main(argv[1:] + prefix)
        return module.main(prefix + argv[1:])
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.command == "serve":
        import sb_server

        server_args = ["--host", args.host, "--port", str(args.port)]
        if args.token is not None:
            server_args.extend(["--token", args.token])
        if args.max_age_days is not None:
            server_args.extend(["--max-age-days", str(args.max_age_days)])
        return sb_server.main(server_args)
    from importlib.metadata import PackageNotFoundError, version
    from sb_report import SOURCES

    try:
        print(f"sbonest {version('sbonest')}")
    except PackageNotFoundError:
        print("sbonest (not installed as a package; running from source)")
    root = Path(__file__).resolve().parent
    for name in SOURCES:
        print(f"{name}: {hashlib.sha256((root / name).read_bytes()).hexdigest()}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
