"""Sideband web runner: upload, preflight, background fit, resume, report and downloads."""
# ruff: noqa: E402 -- Limit numerical libraries before importing them.
import os

for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ.setdefault(name, '1')
os.environ.setdefault('MPLBACKEND', 'Agg')

import io
import json
from pathlib import Path
import tempfile
import time

ROOT = Path(__file__).resolve().parent
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
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        status = client.get(f'/jobs/{job_id}').get_json()
        if not status['running']:
            return status
        time.sleep(0.5)
    raise AssertionError('Fit did not finish in time')


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
    assert job_id in client.get('/jobs').get_json()['jobs']
    text = client.get(f'/jobs/{job_id}/files/fit_result.txt').data.decode()
    assert 'Results using Sideband' in text


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


if __name__ == '__main__':
    check_rejections()
    check_workflow()
    print('PASS: web runner upload, preflight, background fit, resume, report and safe downloads')
