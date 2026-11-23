"""Sideband workflow helpers: demo setup, saved-result reports, design, model comparison and simulation."""
import argparse
from importlib import resources
import json
import sys
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


def simulate(config, config_dir, out, *, seed, truth=None):
    """Write seeded synthetic data at the offsets, sigmas and conditions of a configuration's datasets.

    The generating parameters are the configuration's starting point with ``truth``
    entries (parameter names, ungrouped aliases, or ``kex`` together with ``pB`` for
    two-state models) replacing them. ``out`` must not exist; it receives
    ``data/dataset_N.txt``, ``fit.json`` (the configuration pointed at the new data,
    output prefix ``out/fit``) and ``truth.json``. Returns the three kinds of paths.
    """
    import numpy as np

    from sb_bootstrap import bootstrap_draws
    from sbfit import SidebandModel, _validate_fit_config

    out = Path(out).expanduser().absolute()
    if out.exists() or out.is_symlink():
        raise FileExistsError(f'Simulation target already exists: {out}')
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError('seed must be a nonnegative integer')
    model = SidebandModel(config, config_dir)
    _validate_fit_config(model, config['init'])
    vector = np.array(model.prepare_fit(), dtype=float)
    names = list(model.parameter_names)
    truth = dict(truth or {})
    if 'kex' in truth or 'pB' in truth:
        if 'kbc' in names or not {'kex', 'pB'} <= set(truth) or {'kab', 'kba'} & set(truth):
            raise ValueError('truth kex and pB go together, without kab/kba, for two-state models only')
        kex, pb = truth.pop('kex'), truth.pop('pB')
        truth.update(kab=kex * pb, kba=kex * (1 - pb))
    for name, value in truth.items():
        targets = model.expand_name(name)
        if not targets:
            raise ValueError(f'Unknown truth parameter: {name}')
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not np.isfinite(value):
            raise ValueError(f'truth {name} must be a finite number')
        for target in targets:
            vector[names.index(target)] = value
    data = model._prepare_data()
    draws = bootstrap_draws(model, vector, {'seed': seed}, 0)
    files = []
    for index in range(len(config['datasets'])):
        lines = [f'{model.dataset.fields[index]:.15g}', f'{model.dataset.Ts[index]:.15g}',
                 f'{model.dataset.v1s[index]:.15g} {model.dataset.v1errs[index]:.15g}',
                 f'# offset(ppm) intensity error  (synthetic, seed {seed})']
        for (_, es, offsets, _, sigma), values in zip(data, draws):
            if es.dataset_index != index:
                continue
            label = next(r.label for r in model.dataset.res if any(spec is es for spec in r.estSpecs))
            lines.append(f'# {label} R2a: {es.initr2a:.15g} R2b: {es.initr2b:.15g} dw: {es.initdw:.15g}')
            lines.extend(f'{o:.15g} {v:.15g} {e:.15g}' for o, v, e in zip(offsets, values, sigma))
        files.append((f'data/dataset_{index}.txt', '\n'.join(lines) + '\n'))
    new = json.loads(json.dumps(config))
    new['datasets'] = [name for name, _ in files]
    new['Project Name'] = str(out / 'fit')
    # Only active residues are simulated; waveforms stay where they are.
    active = {r.label for r in model.dataset.res if r.active}
    new['residues'] = [row for row in new.get('residues', []) if row.get('name') in active]
    sideband = new.get('sideband', {})
    for entry in [sideband.get('decoupling', {}), *sideband.get('datasets', [])]:
        if 'waveform_json' in entry:
            entry['waveform_json'] = str((Path(config_dir) / entry['waveform_json']).resolve())
    record = {'synthetic': True, 'seed': seed, 'parameters': dict(zip(names, vector.tolist())),
              'note': 'Gaussian noise with the template sigma added to the model at these parameters.'}
    out.mkdir(parents=True)
    (out / 'data').mkdir()
    for name, text in files:
        (out / name).write_text(text, encoding='utf-8')
    (out / 'fit.json').write_text(json.dumps(new, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    (out / 'truth.json').write_text(json.dumps(record, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    return {'config': str(out / 'fit.json'), 'truth': str(out / 'truth.json'),
            **{f'data_{i}': str(out / name) for i, (name, _) in enumerate(files)}}


def main(argv=None):
    """Command line: init-demo, report, design, compare and simulate; returns the exit status."""
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
    sim = commands.add_parser('simulate', help='Write seeded synthetic data shaped like a configuration\'s datasets')
    sim.add_argument('config_file')
    sim.add_argument('--out', required=True, help='New output directory')
    sim.add_argument('--seed', type=int, required=True, help='Noise seed (nonnegative integer)')
    sim.add_argument('--truth', help='JSON object of generating parameters overriding the starting point')
    args = parser.parse_args(argv)
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
        elif args.command == 'simulate':
            from run import load_config
            config = load_config(args.config_file)
            if not str(config['init'].get('Method', '')).startswith('Sideband'):
                raise ValueError('simulate requires a Sideband init.Method')
            truth = json.loads(Path(args.truth).read_text(encoding='utf-8')) if args.truth else None
            if truth is not None and not isinstance(truth, dict):
                raise ValueError('--truth must hold a JSON object')
            paths = simulate(config, Path(args.config_file).resolve().parent, args.out,
                             seed=args.seed, truth=truth)
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
    return 0


if __name__ == '__main__':
    sys.exit(main())
