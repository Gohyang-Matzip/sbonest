"""Durable checkpoint identity, interruption, and immutable-record checks."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]  # repository root: flat modules live there
sys.path.insert(0, str(ROOT))

import json
import os
import tempfile
from unittest.mock import patch

from numpy.testing import assert_allclose

from sb_checkpoint import Checkpoint


def check_checkpoint():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-checkpoint-')) / 'run'
    identity = {'config': 'a', 'inputs': ['b']}
    try:
        with Checkpoint(folder, identity) as journal:
            journal.save('baseline', {'chi2': 3.5, 'parameters': [1., 2.]})
            with Checkpoint(folder, identity, resume=True):
                raise AssertionError('Concurrent resume obtained the lock')
    except BlockingIOError:
        pass
    assert json.loads((folder / 'baseline.json').read_text())['value']['chi2'] == 3.5
    with Checkpoint(folder, identity, resume=True) as journal:
        assert journal.read('baseline')['parameters'] == [1., 2.]
        journal.save('baseline', journal.read('baseline'))
        try:
            journal.save('baseline', {'chi2': 99})
        except FileExistsError:
            pass
        else:
            raise AssertionError('Completed checkpoint was overwritten')
    for changed in ({'config': 'changed', 'inputs': ['b']}, identity):
        try:
            with Checkpoint(folder, changed, resume=changed != identity):
                pass
        except (ValueError, FileExistsError):
            pass
        else:
            raise AssertionError('Mismatched/new run reused old checkpoint')
    (folder / 'unfinished.pending').write_text('{partial')
    with Checkpoint(folder, identity, resume=True) as journal:
        assert journal.read('baseline')['chi2'] == 3.5
        assert journal.read('absent') is None
        for bad in ('../escape', '/tmp/escape', 'a/b'):
            try:
                journal.save(bad, {})
            except ValueError:
                pass
            else:
                raise AssertionError('Unsafe checkpoint key accepted')
    link = folder.parent / 'link'
    link.symlink_to(folder, target_is_directory=True)
    try:
        with Checkpoint(link, identity, resume=True):
            pass
    except ValueError:
        pass
    else:
        raise AssertionError('Checkpoint directory symlink accepted')
    record = json.loads((folder / 'baseline.json').read_text())
    record['value']['chi2'] = 123456
    (folder / 'baseline.json').write_text(json.dumps(record))
    with Checkpoint(folder, identity, resume=True) as journal:
        try:
            journal.read('baseline')
        except ValueError as exc:
            assert 'checksum' in str(exc)
        else:
            raise AssertionError('Corrupted completed record accepted')


def check_run_resume():
    from sbfit import SidebandModel, run_config
    from test_sb_output import small_config

    folder = Path(tempfile.mkdtemp(prefix='sbonest-resume-'))
    config = small_config(folder)
    config['init']['profile'] = {'v1n_scale': [1.06, 1.10]}
    save = Checkpoint.save

    def interrupt(journal, key, value):
        save(journal, key, value)
        if key == 'profile-v1n_scale-0':
            raise KeyboardInterrupt

    try:
        with patch.object(Checkpoint, 'save', interrupt):
            run_config(config, no_pdf=True)
    except KeyboardInterrupt:
        pass
    else:
        raise AssertionError('No durable per-point completion hook')
    checkpoint = folder / 'fit_checkpoint'
    baseline = json.loads((checkpoint / 'baseline.json').read_text())['value']
    assert not (folder / 'fit_result.json').exists()
    with patch.object(SidebandModel, 'fit', side_effect=AssertionError('Completed baseline fitted again')):
        run_config(config, no_pdf=True, resume=True)
    result = json.loads((folder / 'fit_result.json').read_text())
    assert result['success'] and len(result['profiles']['v1n_scale']) == 2
    assert all(p['success'] for p in result['profiles']['v1n_scale'])
    assert_allclose(result['covariance'], baseline['snapshot']['covariance'])
    originals = {p.name: p.read_bytes() for p in folder.iterdir() if p.is_file()}
    with patch.object(SidebandModel, 'fit', side_effect=AssertionError('Finished run refitted')):
        run_config(config, no_pdf=True, resume=True)
    assert originals == {p.name: p.read_bytes() for p in folder.iterdir() if p.is_file()}
    changed = {**config, 'init': {**config['init'], 'max_nfev': 123}}
    try:
        run_config(changed, no_pdf=True, resume=True)
    except ValueError as exc:
        assert 'identity' in str(exc).lower()
    else:
        raise AssertionError('Changed config resumed a different run')


def check_export_resume_and_failure():
    from sbfit import SidebandModel, run_config
    from test_sb_output import small_config

    folder = Path(tempfile.mkdtemp(prefix='sbonest-export-resume-'))
    config = small_config(folder)
    link = os.link

    def interrupted(source, destination, *args, **kwargs):
        link(source, destination, *args, **kwargs)
        if Path(destination) == folder / 'fit_result.json':
            raise KeyboardInterrupt

    try:
        with patch('sb_checkpoint.os.link', interrupted):
            run_config(config, no_pdf=True)
    except KeyboardInterrupt:
        pass
    else:
        raise AssertionError('Export interruption not exercised')
    before = (folder / 'fit_result.json').read_bytes()
    assert not (folder / 'fit_predictions.csv').exists()
    with patch.object(SidebandModel, 'fit', side_effect=AssertionError('Interrupted export refitted')):
        run_config(config, no_pdf=True, resume=True)
    assert (folder / 'fit_result.json').read_bytes() == before
    assert (folder / 'fit_predictions.csv').exists()
    (folder / 'fit_predictions.csv').write_text('user edited output\n')
    try:
        run_config(config, no_pdf=True, resume=True)
    except FileExistsError:
        pass
    else:
        raise AssertionError('Resume overwrote modified output')
    assert (folder / 'fit_predictions.csv').read_text() == 'user edited output\n'
    config['Project Name'] = str(folder / 'failed')
    config['init']['max_nfev'] = 1
    try:
        run_config(config, no_pdf=True)
    except RuntimeError:
        pass
    else:
        raise AssertionError('Unconverged fit succeeded')
    failure = json.loads((folder / 'failed_result.json').read_text())
    assert not failure['success'] and 'optimizer' in failure and 'provenance' in failure
    with patch.object(SidebandModel, 'fit', side_effect=AssertionError('Completed failed fit retried')):
        try:
            run_config(config, no_pdf=True, resume=True)
        except RuntimeError as exc:
            assert 'failed fit' in str(exc)
        else:
            raise AssertionError('Failed checkpoint reported success')


def check_relative_prefix():
    """A cwd-relative Project Name (as in demo_sideband.py) must stage and publish."""
    from sbfit import SidebandModel, run_config
    from test_sb_output import small_config

    folder = Path(tempfile.mkdtemp(prefix='sbonest-relative-'))
    config = small_config(folder)
    config['Project Name'] = 'relative'
    previous = os.getcwd()
    os.chdir(folder)
    try:
        run_config(config, no_pdf=True)
        result = json.loads((folder / 'relative_result.json').read_text())
        assert result['success']
        assert Path(result['checkpoint']).resolve() == (folder / 'relative_checkpoint').resolve()
        assert (folder / 'relative_predictions.csv').exists()
        with patch.object(SidebandModel, 'fit', side_effect=AssertionError('Relative-prefix run refitted')):
            run_config(config, no_pdf=True, resume=True)
    finally:
        os.chdir(previous)


if __name__ == '__main__':
    check_checkpoint()
    check_run_resume()
    check_export_resume_and_failure()
    check_relative_prefix()
    print('PASS: checkpoint identity, locking, immutable records and partial-write recovery')
