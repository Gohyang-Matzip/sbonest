"""Sideband workflow helpers: demo setup, saved-result reports, design and model comparison."""
import argparse
from importlib import resources
import json
from pathlib import Path


def init_demo(out):
    """Copy bundled synthetic two-RF inputs into a new, exclusively owned folder."""
    out = Path(out).expanduser().absolute()
    if out.exists() or out.is_symlink():
        raise FileExistsError(f'Demo target already exists: {out}')
    root = resources.files('sbonest_data').joinpath('sideband_auto_H')
    config = json.loads(root.joinpath('two_RF.json').read_bytes())
    sources = [root.joinpath(name) for name in config['datasets']]
    # Read every bundled input before creating any output.
    contents = [(source.name, source.read_bytes()) for source in sources]
    config['datasets'] = [f'data/{name}' for name, _ in contents]
    config['Project Name'] = str(out / 'fit')
    out.mkdir(parents=True)
    (out / 'data').mkdir()
    for name, content in contents:
        with (out / 'data' / name).open('xb') as stream:
            stream.write(content)
    path = out / 'fit.json'
    with path.open('x', encoding='utf-8') as stream:
        stream.write(json.dumps(config, indent=2, allow_nan=False) + '\n')
    return path


def main():
    """Command line: init-demo, report, design and compare."""
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    demo = commands.add_parser('init-demo', help='Copy bundled synthetic two-RF data')
    demo.add_argument('--out', required=True, help='New demo directory (must not exist)')
    report = commands.add_parser('report', help='Report saved JSON/CSV without refitting')
    report.add_argument('result_json')
    report.add_argument('--out', required=True, help='New report output prefix')
    report.add_argument('--predictions', help='Explicit predictions CSV, if renamed')
    design = commands.add_parser('design', help='Expected errors of planned acquisitions at a known truth (no fitting)')
    design.add_argument('design_json')
    design.add_argument('--out', required=True, help='New output directory')
    design.add_argument('--workers', type=int, default=1, help='Worker processes for the Jacobian columns')
    compare = commands.add_parser('compare', help='Fit shared-exchange and per-residue models and compare them')
    compare.add_argument('config_file')
    compare.add_argument('--out', required=True, help='New output directory')
    compare.add_argument('--workers', type=int, default=1, help='Worker processes for every sub-fit')
    compare.add_argument('--pdf', action='store_true', help='Also write fit PDFs for every sub-fit')
    compare.add_argument('--models', nargs='+', metavar='MODEL',
                         help='Compare Sideband models on the same data instead of residues, e.g. Sideband Sideband_3st_Linear')
    compare.add_argument('--h-ppm-c', nargs='+', metavar='LABEL=ppm',
                         help='State-C proton shifts for three-state models (default: the state-B shift)')
    args = parser.parse_args()
    if getattr(args, 'workers', 1) < 1:
        parser.error('--workers must be a positive integer')
    try:
        if args.command == 'init-demo':
            path = init_demo(args.out)
            print(f'Config: {path}')
            config = json.loads(path.read_text(encoding='utf-8'))
            for name in config['datasets']:
                print(f'Data: {path.parent / name}')
            print(f'Fit output prefix: {config["Project Name"]}')
        elif args.command == 'report':
            from sb_report import regenerate_report
            paths = regenerate_report(args.result_json, args.out, args.predictions)
            for label, path in paths.items():
                print(f'{label}: {path}')
        elif args.command == 'design':
            from sb_design import run_design
            paths = run_design(args.design_json, args.out, workers=args.workers)
            for label, path in paths.items():
                print(f'{label}: {path}')
        else:
            from run import load_config
            from sb_compare import compare_models, run_comparison
            config = load_config(args.config_file)
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
                paths = compare_models(config, Path(args.config_file).resolve().parent, args.out,
                                       models=tuple(args.models), h_ppm_c=shifts, no_pdf=not args.pdf,
                                       workers=args.workers)
            else:
                paths = run_comparison(config, Path(args.config_file).resolve().parent, args.out,
                                       no_pdf=not args.pdf, workers=args.workers)
            for label, path in paths.items():
                print(f'{label}: {path}')
    except (ValueError, KeyError, OSError, RuntimeError) as exc:
        parser.exit(1, f'Error: {exc}\n')


if __name__ == '__main__':
    main()
