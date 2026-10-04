"""Sideband web interface: upload, preflight, background fitting and reports.

Jobs live under SB_JOBS/<id>/ next to this file (or SBONEST_JOBS_DIR). A job
holds the uploaded configuration (rewritten to the uploaded data names and a
fixed output prefix), its data files, the fit log and every output of run.py,
so each job is a complete, reproducible run with its own checkpoint. Fits run
as background subprocesses of run.py; status is read from the process and the
checkpoint records, and --resume continues an interrupted job. Only files
inside a job directory can be downloaded.
"""
import json
import os
from pathlib import Path
import secrets
import subprocess
import sys
import threading
import time

from flask import Flask, abort, jsonify, render_template_string, request, send_from_directory
from werkzeug.utils import secure_filename

BASE_DIR = Path(__file__).resolve().parent
JOBS_DIR = Path(os.environ.get("SBONEST_JOBS_DIR", BASE_DIR / "SB_JOBS"))
MAX_WORKERS = max(1, os.cpu_count() or 1)

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 64 * 1024 * 1024
_processes = {}
_lock = threading.Lock()


def _job_dir(job_id):
    if not job_id or not all(c.isalnum() or c in "-_" for c in job_id):
        abort(400, "Invalid job id")
    folder = (JOBS_DIR / job_id).resolve()
    if folder.parent != JOBS_DIR.resolve() or not folder.is_dir():
        abort(404, "Unknown job")
    return folder


def _python_command(folder, script, *args):
    env = dict(os.environ, PYTHONPATH=str(BASE_DIR), MPLBACKEND="Agg")
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        env.setdefault(name, "1")
    return [sys.executable, str(BASE_DIR / script), *args], env


def _read_json(path):
    try:
        return json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def _checkpoint_status(folder):
    checkpoint = folder / "fit_checkpoint"
    if not checkpoint.is_dir():
        return {"exists": False}
    names = sorted(p.name[:-5] for p in checkpoint.glob("*.json"))
    counts = {}
    for name in names:
        head = name.rsplit("-", 1)[0] if name[-1].isdigit() else name
        counts[head] = counts.get(head, 0) + 1
    return {"exists": True, "records": counts, "complete": "complete" in names,
            "failed": "failure" in names, "baseline": "baseline" in names}


def _status(job_id, folder):
    with _lock:
        process = _processes.get(job_id)
    running = process is not None and process.poll() is None
    returncode = None if process is None or running else process.returncode
    log = folder / "fit.log"
    tail = ""
    if log.exists():
        text = log.read_text(encoding="utf-8", errors="replace")
        tail = text[-4000:]
    result = _read_json(folder / "fit_result.json")
    outputs = sorted(p.name for p in folder.iterdir()
                     if p.is_file() and p.name.startswith(("fit", "report", "check")) and p.suffix in (".json", ".txt", ".csv", ".pdf", ".log"))
    summary = None
    if result:
        summary = {key: result.get(key) for key in ("success", "chi2", "dof", "n_points", "n_parameters",
                                                     "kex", "pB", "derived_se", "warnings", "model", "message")}
    return {"job_id": job_id, "running": running, "returncode": returncode,
            "checkpoint": _checkpoint_status(folder), "result": summary, "outputs": outputs,
            "log_tail": tail, "meta": _read_json(folder / "job.json")}


PAGE = """<!doctype html><html lang="en"><head><meta charset="utf-8"><title>SBONEST Sideband runner</title>
<style>body{font-family:sans-serif;margin:20px;max-width:980px}fieldset{margin-bottom:14px}
pre{background:#f4f4f4;padding:8px;overflow:auto;max-height:320px}label{display:block;margin:4px 0}
.ok{color:#1a7f37}.bad{color:#b42318}</style></head><body>
<h1>SBONEST Sideband runner</h1>
<p>Upload a Sideband configuration and its data files. The configuration's <code>datasets</code>
entries are matched to the uploaded file names; the output prefix becomes <code>fit</code>.
Each job keeps its own checkpoint, so an interrupted fit can be resumed.</p>
<form id="upload" enctype="multipart/form-data">
<fieldset><legend>1. Inputs</legend>
<label>Configuration JSON <input type="file" name="config" accept=".json" required></label>
<label>Data files (all datasets) <input type="file" name="data" multiple required></label>
<label>Waveform JSON (optional) <input type="file" name="waveform" accept=".json"></label>
</fieldset>
<fieldset><legend>2. Options</legend>
<label><input type="checkbox" name="no_pdf" checked> Skip PDF figures</label>
<label>Workers <input type="number" name="workers" min="1" max="{{ max_workers }}" value="1"></label>
</fieldset>
<button type="submit">Create job and run preflight check</button>
</form>
<div id="job" hidden><h2>Job <span id="jobid"></span></h2>
<pre id="check"></pre>
<button id="fit">Start fit</button> <button id="resume">Resume</button> <button id="report">Make report</button>
<span id="state"></span><pre id="log"></pre><div id="files"></div></div>
<script>
let job=null;
const q=s=>document.querySelector(s);
async function post(url,body){const r=await fetch(url,{method:'POST',body});return r.json();}
q('#upload').addEventListener('submit',async e=>{e.preventDefault();
 const d=new FormData(e.target);const j=await post('/jobs',d);
 if(!j.job_id){q('#check').textContent=JSON.stringify(j,null,2);return;}
 job=j.job_id;q('#jobid').textContent=job;q('#job').hidden=false;
 q('#check').textContent=JSON.stringify(j.check,null,2);poll();});
q('#fit').addEventListener('click',async()=>{await post('/jobs/'+job+'/fit');poll();});
q('#resume').addEventListener('click',async()=>{await post('/jobs/'+job+'/resume');poll();});
q('#report').addEventListener('click',async()=>{const j=await post('/jobs/'+job+'/report');q('#log').textContent=JSON.stringify(j,null,2);poll();});
async function poll(){if(!job)return;const s=await (await fetch('/jobs/'+job)).json();
 q('#state').textContent=s.running?'running':(s.result?('finished: '+(s.result.success?'success':'failed')):'idle');
 q('#state').className=s.result&&s.result.success?'ok':(s.result?'bad':'');
 q('#log').textContent=(s.result?JSON.stringify(s.result,null,2)+'\n':'')+s.log_tail;
 q('#files').innerHTML=s.outputs.map(n=>'<a href="/jobs/'+job+'/files/'+n+'">'+n+'</a>').join(' | ');
 if(s.running)setTimeout(poll,3000);}
</script></body></html>"""


@app.route("/", methods=["GET"])
def index():
    return render_template_string(PAGE, max_workers=MAX_WORKERS)


@app.route("/jobs", methods=["POST"])
def create_job():
    config_file = request.files.get("config")
    data_files = [f for f in request.files.getlist("data") if f and f.filename]
    if config_file is None or not config_file.filename or not data_files:
        return jsonify({"error": "Upload a configuration JSON and at least one data file"}), 400
    try:
        config = json.load(config_file.stream)
    except ValueError:
        return jsonify({"error": "Configuration is not valid JSON"}), 400
    if not isinstance(config, dict) or not isinstance(config.get("datasets"), list):
        return jsonify({"error": "Configuration needs a datasets list"}), 400
    if not str(config.get("init", {}).get("Method", "")).startswith("Sideband"):
        return jsonify({"error": "This runner accepts Sideband configurations only"}), 400
    uploaded = {}
    for item in data_files:
        name = secure_filename(item.filename)
        if not name or name in uploaded:
            return jsonify({"error": f"Invalid or duplicate data file name: {item.filename}"}), 400
        uploaded[name] = item
    expected = [secure_filename(Path(str(p)).name) for p in config["datasets"]]
    missing = [name for name in expected if name not in uploaded]
    if missing:
        return jsonify({"error": "Configuration datasets not uploaded: " + ", ".join(missing)}), 400
    waveform = request.files.get("waveform")
    workers = request.form.get("workers", "1")
    try:
        workers = int(workers)
        if not 1 <= workers <= MAX_WORKERS:
            raise ValueError
    except ValueError:
        return jsonify({"error": f"workers must be an integer between 1 and {MAX_WORKERS}"}), 400
    no_pdf = request.form.get("no_pdf") not in (None, "", "false", "0", "off")
    job_id = time.strftime("%Y%m%d-%H%M%S-") + secrets.token_hex(4)
    folder = JOBS_DIR / job_id
    folder.mkdir(parents=True, exist_ok=False)
    for name, item in uploaded.items():
        item.save(folder / name)
    config["datasets"] = expected
    config["Project Name"] = "fit"
    if waveform is not None and waveform.filename:
        wname = secure_filename(waveform.filename)
        waveform.save(folder / wname)
        sideband = config.get("sideband", {})
        for entry in [sideband.get("decoupling", {}), *sideband.get("datasets", [])]:
            if isinstance(entry, dict) and "waveform_json" in entry:
                entry["waveform_json"] = wname
    (folder / "config.json").write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    (folder / "job.json").write_text(json.dumps({"job_id": job_id, "workers": workers, "no_pdf": no_pdf,
                                                 "created": time.strftime("%Y-%m-%dT%H:%M:%S")}, indent=2) + "\n")
    args = ["config.json", "--check", "--identifiability", "--workers", str(workers)] + (["--no-pdf"] if no_pdf else [])
    command, env = _python_command(folder, "run.py", *args)
    run = subprocess.run(command, cwd=folder, env=env, capture_output=True, text=True, timeout=600)
    (folder / "check.log").write_text(run.stdout + run.stderr, encoding="utf-8")
    try:
        check = json.loads(run.stdout)
    except ValueError:
        check = {"valid": False, "errors": [run.stderr.strip() or "preflight produced no JSON"]}
    (folder / "check.json").write_text(json.dumps(check, indent=2) + "\n", encoding="utf-8")
    return jsonify({"job_id": job_id, "check": check})


def _start(job_id, folder, resume):
    meta = _read_json(folder / "job.json") or {}
    with _lock:
        process = _processes.get(job_id)
        if process is not None and process.poll() is None:
            return jsonify({"error": "A fit is already running for this job"}), 409
        args = ["config.json", "--workers", str(meta.get("workers", 1))]
        if meta.get("no_pdf", True):
            args.append("--no-pdf")
        if resume:
            args.append("--resume")
        command, env = _python_command(folder, "run.py", *args)
        log = (folder / "fit.log").open("a", encoding="utf-8")
        log.write(f"\n$ {' '.join(command[1:])}\n")
        log.flush()
        _processes[job_id] = subprocess.Popen(command, cwd=folder, env=env, stdout=log, stderr=subprocess.STDOUT)
    return jsonify({"job_id": job_id, "started": True, "resume": resume})


@app.route("/jobs/<job_id>/fit", methods=["POST"])
def start_fit(job_id):
    folder = _job_dir(job_id)
    if (folder / "fit_checkpoint").exists():
        return jsonify({"error": "This job already has a checkpoint; use resume or create a new job"}), 409
    return _start(job_id, folder, resume=False)


@app.route("/jobs/<job_id>/resume", methods=["POST"])
def resume_fit(job_id):
    folder = _job_dir(job_id)
    if not (folder / "fit_checkpoint").exists():
        return jsonify({"error": "Nothing to resume; start a fit first"}), 409
    return _start(job_id, folder, resume=True)


@app.route("/jobs/<job_id>/report", methods=["POST"])
def make_report(job_id):
    folder = _job_dir(job_id)
    if not (folder / "fit_result.json").exists():
        return jsonify({"error": "No fit result to report yet"}), 409
    index = 1
    while (folder / f"report_{index:02d}_summary.json").exists():
        index += 1
    prefix = f"report_{index:02d}"
    command, env = _python_command(folder, "sb_workflow.py", "report", "fit_result.json", "--out", prefix)
    run = subprocess.run(command, cwd=folder, env=env, capture_output=True, text=True, timeout=600)
    if run.returncode != 0:
        return jsonify({"error": run.stderr.strip() or run.stdout.strip()}), 500
    return jsonify({"job_id": job_id, "prefix": prefix,
                    "files": [f"{prefix}_summary.json", f"{prefix}_summary.txt", f"{prefix}.pdf"]})


@app.route("/jobs/<job_id>", methods=["GET"])
def job_status(job_id):
    folder = _job_dir(job_id)
    return jsonify(_status(job_id, folder))


@app.route("/jobs", methods=["GET"])
def list_jobs():
    JOBS_DIR.mkdir(parents=True, exist_ok=True)
    jobs = sorted(p.name for p in JOBS_DIR.iterdir() if p.is_dir())
    return jsonify({"jobs": jobs})


@app.route("/jobs/<job_id>/files/<path:filename>", methods=["GET"])
def download(job_id, filename):
    folder = _job_dir(job_id)
    safe = secure_filename(Path(filename).name)
    if not safe or safe != filename or not (folder / safe).is_file():
        abort(404)
    return send_from_directory(folder, safe, as_attachment=True)


def main(argv=None):
    import argparse

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=5050)
    args = parser.parse_args(argv)
    JOBS_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Sideband runner at http://{args.host}:{args.port}/  (jobs in {JOBS_DIR})")
    app.run(host=args.host, port=args.port, debug=False)


if __name__ == "__main__":
    main()
