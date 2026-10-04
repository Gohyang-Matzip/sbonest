"""Sideband web runner: upload, preflight, background fit, resume, report and downloads."""
# ruff: noqa: E402 -- Limit numerical libraries before importing them.
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]  # repository root: flat modules live there
sys.path.insert(0, str(ROOT))

import os

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(name, '1')
os.environ.setdefault('MPLBACKEND', 'Agg')

import io
import json
import tempfile
import time

JOBS = Path(tempfile.mkdtemp(prefix='sbonest-web-jobs-'))
os.environ['SBONEST_JOBS_DIR'] = str(JOBS)

import sb_server
from sb_compare import _residue_config
from test_sb_compare import write_synthetic


def upload(client, folder, config, **form):
    data = {'config': (io.BytesIO(json.dumps(config).encode()), 'fit.json')}
    files = [(open(folder / name, 'rb'), name) for name in config['datasets']]
    data['data'] = files
    data.update(form)
    response = client.post('/jobs', data=data, content_type='multipart/form-data')
    for handle, _ in files:
        handle.close()
    return response


def wait_idle(client, job_id, timeout=240):
    process = sb_server._processes[job_id]
    process.wait(timeout=timeout)
    status = client.get(f'/jobs/{job_id}').get_json()
    assert not status['running']
    return status


def check_workflow():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-web-data-'))
    path = write_synthetic(folder, {'A1': (15., 285.), 'G2': (15., 285.), 'S3': (15., 285.)})
    config = json.loads(path.read_text())
    config = _residue_config(config, ['A1', 'G2', 'S3'], {'A1'})
    config['init']['vary'] = ['kab', 'kba', 'v1n_scale', 'A1.R2b']
    client = sb_server.app.test_client()
    assert b'SBONEST Sideband runner' in client.get('/').data
    created = upload(client, folder, config, workers='2', no_pdf='on')
    assert created.status_code == 200, created.get_json()
    job = created.get_json()
    job_id = job['job_id']
    assert job['check']['valid'] and job['check']['n_parameters'] == 4
    assert 'identifiability' in job['check']
    saved = json.loads((JOBS / job_id / 'config.json').read_text())
    assert saved['Project Name'] == 'fit' and saved['datasets'] == config['datasets']
    assert (JOBS / job_id / 'check.json').exists()
    status = client.get(f'/jobs/{job_id}').get_json()
    assert not status['running'] and status['checkpoint'] == {'exists': False} and status['meta']['workers'] == 2
    assert client.post(f'/jobs/{job_id}/resume').status_code == 409
    assert client.post(f'/jobs/{job_id}/report').status_code == 409
    started = client.post(f'/jobs/{job_id}/fit')
    assert started.status_code == 200 and started.get_json()['started']
    assert client.post(f'/jobs/{job_id}/fit').status_code == 409
    status = wait_idle(client, job_id)
    assert status['returncode'] == 0, status['log_tail'][-1500:]
    assert status['result']['success'] and status['checkpoint']['complete'] and status['checkpoint']['baseline']
    assert 'fit_result.json' in status['outputs'] and 'fit_predictions.csv' in status['outputs']
    resumed = client.post(f'/jobs/{job_id}/resume')
    assert resumed.status_code == 200
    status = wait_idle(client, job_id)
    assert status['returncode'] == 0 and 'Resumed completed outputs' in status['log_tail']
    report = client.post(f'/jobs/{job_id}/report')
    assert report.status_code == 200, report.get_json()
    files = report.get_json()['files']
    assert files[0] == 'report_01_summary.json'
    download = client.get(f'/jobs/{job_id}/files/{files[2]}')
    assert download.status_code == 200 and download.data.startswith(b'%PDF-')
    assert client.get(f'/jobs/{job_id}/files/../config.json').status_code == 404
    assert client.get('/jobs/../etc').status_code in (400, 404)
    assert client.get('/jobs/missing-job').status_code == 404
    assert job_id in [row['job_id'] for row in client.get('/jobs').get_json()['jobs']]
    text = client.get(f'/jobs/{job_id}/files/fit_result.txt').data.decode()
    assert 'Results using Sideband' in text
    # PNG preview is rendered from the predictions CSV and cached.
    png = client.get(f'/jobs/{job_id}/preview.png?t=1')
    assert png.status_code == 200 and png.data.startswith(b'\x89PNG')
    assert (JOBS / job_id / 'preview.png').exists()
    listing = client.get('/jobs').get_json()
    entry = next(row for row in listing['jobs'] if row['job_id'] == job_id)
    assert entry['finished'] and not entry['running'] and entry['age_days'] >= 0
    # Archiving moves the job directory without deleting anything.
    archived = client.post(f'/jobs/{job_id}/archive')
    assert archived.status_code == 200
    assert not (JOBS / job_id).exists() and (JOBS / 'archive' / job_id / 'fit_result.json').exists()
    assert client.get(f'/jobs/{job_id}').status_code == 404
    assert job_id in client.get('/jobs').get_json()['archived']


def check_analyses_token_and_expiry():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-web-opts-'))
    path = write_synthetic(folder, {'A1': (15., 285.), 'G2': (15., 285.), 'S3': (15., 285.)})
    config = _residue_config(json.loads(path.read_text()), ['A1', 'G2', 'S3'], {'A1'})
    config['init']['vary'] = ['kab', 'kba', 'v1n_scale', 'A1.R2b']
    client = sb_server.app.test_client()
    # Analysis options become init keys; invalid combinations are rejected before any job exists.
    bad = upload(client, folder, config, multistart_random='2')
    assert bad.status_code == 400 and 'multistart_seed' in bad.get_json()['error']
    bad = upload(client, folder, config, profile_kex='300 -5')
    assert bad.status_code == 400
    created = upload(client, folder, config, multistart_random='1', multistart_seed='3', profile_kex='290 300',
                     profile_interval='kex', bootstrap_replicates='2', bootstrap_seed='9', no_pdf='on')
    assert created.status_code == 200, created.get_json()
    job_id = created.get_json()['job_id']
    assert created.get_json()['analyses'] == ['bootstrap', 'multistart', 'profile', 'profile_interval']
    saved = json.loads((JOBS / job_id / 'config.json').read_text())['init']
    assert saved['multistart'] == {'random_starts': 1, 'seed': 3} and saved['profile'] == {'kex': [290., 300.]}
    assert saved['profile_interval'] == {'parameters': ['kex']} and saved['bootstrap'] == {'replicates': 2, 'seed': 9}
    assert client.post(f'/jobs/{job_id}/fit').status_code == 200
    status = wait_idle(client, job_id)
    assert status['returncode'] == 0, status['log_tail'][-2000:]
    assert set(status['result']['analyses']) == {'multistart', 'profiles', 'profile_intervals', 'bootstrap'}
    assert status['result']['reduced_chi2'] is not None
    # Token protection applies to every endpoint except the page.
    sb_server.app.config['SBONEST_TOKEN'] = 'secret-token'
    try:
        assert client.get('/').status_code == 200
        assert client.get('/jobs').status_code == 401 and client.get(f'/jobs/{job_id}').status_code == 401
        assert client.get('/jobs', headers={'X-SBONEST-Token': 'wrong'}).status_code == 401
        assert client.get('/jobs', headers={'X-SBONEST-Token': 'secret-token'}).status_code == 200
        assert client.get(f'/jobs/{job_id}?token=secret-token').status_code == 200
        png = client.get(f'/jobs/{job_id}/preview.png?token=secret-token&t=1')
        assert png.status_code == 200 and png.data.startswith(b'\x89PNG')
        for action in ('fit', 'resume', 'report', 'archive'):
            rejected = client.post(f'/jobs/{job_id}/{action}', headers={'X-SBONEST-Token': 'wrong'})
            assert rejected.status_code == 401 and rejected.get_json()['error']
        assert (JOBS / job_id).is_dir()
        assert b'access token' in client.get('/').data
    finally:
        sb_server.app.config['SBONEST_TOKEN'] = None
    # Expiry archives only idle jobs older than the limit.
    assert sb_server.archive_expired(None) == []
    old = time.time() - 3 * 86400
    os.utime(JOBS / job_id, (old, old))
    sb_server.app.config['SBONEST_MAX_AGE_DAYS'] = 2
    try:
        moved = client.post('/jobs/archive-expired').get_json()
        assert job_id in moved['archived'] and (JOBS / 'archive' / job_id).exists()
    finally:
        sb_server.app.config['SBONEST_MAX_AGE_DAYS'] = None
    assert sb_server.analysis_settings({'profile_interval': ['pB', 'nonsense']}) == {'profile_interval': {'parameters': ['pB']}}
    assert sb_server.analysis_settings({}) == {}


def check_rejections():
    folder = Path(tempfile.mkdtemp(prefix='sbonest-web-bad-'))
    path = write_synthetic(folder, {'A1': (15., 285.), 'G2': (15., 285.), 'S3': (15., 285.)})
    config = json.loads(path.read_text())
    client = sb_server.app.test_client()
    missing = client.post('/jobs', data={'config': (io.BytesIO(json.dumps(config).encode()), 'fit.json')},
                          content_type='multipart/form-data')
    assert missing.status_code == 400
    incomplete = client.post('/jobs', data={'config': (io.BytesIO(json.dumps(config).encode()), 'fit.json'),
                                            'data': [(open(folder / config['datasets'][0], 'rb'), config['datasets'][0])]},
                             content_type='multipart/form-data')
    assert incomplete.status_code == 400 and 'not uploaded' in incomplete.get_json()['error']
    onest = dict(config, init={'Method': 'Baldwin'})
    onest.pop('sideband')
    wrong = upload(client, folder, onest)
    assert wrong.status_code == 400 and 'Sideband' in wrong.get_json()['error']
    assert upload(client, folder, config, workers='0').status_code == 400
    assert not any(p.is_dir() and (p / 'config.json').exists() and json.loads((p / 'config.json').read_text()).get('init', {}).get('Method') == 'Baldwin'
                   for p in JOBS.iterdir())


def check_running_archive():
    """Rejected archives preserve files while a real child is awaiting input."""
    import subprocess
    import sys

    folder = JOBS / 'running-archive'
    folder.mkdir()
    marker = folder / 'retained.txt'
    marker.write_text('retained')
    process = subprocess.Popen([sys.executable, '-c', 'import sys; sys.stdin.buffer.read()'],
                               stdin=subprocess.PIPE)
    sb_server._processes[folder.name] = process
    try:
        client = sb_server.app.test_client()
        for _ in range(2):
            response = client.post(f'/jobs/{folder.name}/archive')
            assert response.status_code == 409 and response.get_json()['error']
            assert client.get(f'/jobs/{folder.name}').get_json()['running']
            assert marker.read_text() == 'retained'
    finally:
        process.communicate(timeout=10)
        sb_server._processes.pop(folder.name)
    assert process.returncode == 0
    response = client.post(f'/jobs/{folder.name}/archive')
    assert response.status_code == 200
    assert (JOBS / 'archive' / folder.name / marker.name).read_text() == 'retained'


if __name__ == '__main__':
    check_running_archive()
    check_rejections()
    check_workflow()
    check_analyses_token_and_expiry()
    print('PASS: web runner upload, preflight, background fit, resume, report, preview, archive, analyses, token and expiry')
