"""Prepare a portable Sideband demo or report an existing fit without refitting."""
import argparse
import json
from pathlib import Path


def init_demo(out):
    """Copy bundled synthetic two-RF inputs into a new, exclusively owned folder."""
    out = Path(out).expanduser().absolute()
    if out.exists() or out.is_symlink():
        raise FileExistsError(f'Demo target already exists: {out}')
    root = Path(__file__).resolve().parent
    template = root / 'example/sideband_auto_H/two_RF.json'
    config = json.loads(template.read_text(encoding='utf-8'))
    sources = [(template.parent / name).resolve() for name in config['datasets']]
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
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    demo = commands.add_parser('init-demo', help='Copy bundled synthetic two-RF data')
    demo.add_argument('--out', required=True, help='New demo directory (must not exist)')
    report = commands.add_parser('report', help='Report saved JSON/CSV without refitting')
    report.add_argument('result_json')
    report.add_argument('--out', required=True, help='New report output prefix')
    report.add_argument('--predictions', help='Explicit predictions CSV, if renamed')
    args = parser.parse_args()
    try:
        if args.command == 'init-demo':
            path = init_demo(args.out)
            print(f'Config: {path}')
            config = json.loads(path.read_text(encoding='utf-8'))
            for name in config['datasets']:
                print(f'Data: {path.parent / name}')
            print(f'Fit output prefix: {config["Project Name"]}')
        else:
            from sb_report import regenerate_report
            paths = regenerate_report(args.result_json, args.out, args.predictions)
            for label, path in paths.items():
                print(f'{label}: {path}')
    except (ValueError, KeyError, OSError, RuntimeError) as exc:
        parser.exit(1, f'Error: {exc}\n')


if __name__ == '__main__':
    main()
