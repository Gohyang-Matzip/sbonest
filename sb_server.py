"""Sideband web interface: upload, preflight, background fitting, reports and previews.

Jobs live under SB_JOBS/<id>/ next to this file (or SBONEST_JOBS_DIR). A job
holds the uploaded configuration (rewritten to the uploaded data names and a
fixed output prefix), its data files, the fit log and every output of run.py,
so each job is a complete, reproducible run with its own checkpoint. Fits run
as background subprocesses of run.py; status is read from the process and the
checkpoint records, and --resume continues an interrupted job. Only files
inside a job directory can be downloaded.

Access control: when SBONEST_TOKEN (or --token) is set, every request except
the page itself must carry the token in the X-SBONEST-Token header or a
``token`` query/form field. Jobs older than --max-age-days are moved to
SB_JOBS/archive/ when the server starts or when /jobs/archive-expired is called;
nothing is deleted.
"""
import json
import os
from pathlib import Path
import secrets
import shutil
import subprocess
import sys
import threading
import time

from flask import Flask, abort, jsonify, render_template_string, request, send_from_directory
from werkzeug.utils import secure_filename

BASE_DIR = Path(__file__).resolve().parent
JOBS_DIR = Path(os.environ.get("SBONEST_JOBS_DIR", BASE_DIR / "SB_JOBS"))
ARCHIVE_NAME = "archive"
MAX_WORKERS = max(1, os.cpu_count() or 1)
ANALYSIS_FIELDS = ("multistart_random", "multistart_seed", "profile_kex", "profile_interval",
                   "bootstrap_replicates", "bootstrap_seed")

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 64 * 1024 * 1024
app.config["SBONEST_TOKEN"] = os.environ.get("SBONEST_TOKEN") or None
app.config["SBONEST_MAX_AGE_DAYS"] = None
_processes = {}
_lock = threading.Lock()


@app.before_request
def _check_token():
    token = app.config.get("SBONEST_TOKEN")
    if not token or request.endpoint in (None, "index", "static"):
        return None
    supplied = request.headers.get("X-SBONEST-Token") or request.args.get("token") or request.form.get("token")
    if not supplied or not secrets.compare_digest(str(supplied), token):
        return jsonify({"error": "Missing or invalid token"}), 401
    return None


def _job_dir(job_id):
    """Resolve a job id to its directory; reject traversal, archived and unknown ids."""
    if not job_id or not all(c.isalnum() or c in "-_" for c in job_id) or job_id == ARCHIVE_NAME:
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
                     if p.is_file() and p.name.startswith(("fit", "report", "check", "preview"))
                     and p.suffix in (".json", ".txt", ".csv", ".pdf", ".log", ".png"))
    summary = None
    if result:
        summary = {key: result.get(key) for key in ("success", "chi2", "dof", "n_points", "n_parameters",
                                                     "kex", "pB", "derived_se", "warnings", "model", "message")}
        diagnostics = result.get("residual_diagnostics") or {}
        summary["reduced_chi2"] = diagnostics.get("reduced_chi2")
        summary["analyses"] = [key for key in ("multistart", "profiles", "profile_intervals", "bootstrap") if key in result]
    meta = _read_json(folder / "job.json") or {}
    return {"job_id": job_id, "running": running, "returncode": returncode,
            "checkpoint": _checkpoint_status(folder), "result": summary, "outputs": outputs,
            "log_tail": tail, "meta": meta, "age_days": _age_days(folder),
            "preview": (folder / "preview.png").exists() or (folder / "fit_predictions.csv").exists()}


def _age_days(folder):
    try:
        return (time.time() - folder.stat().st_mtime) / 86400.0
    except OSError:
        return None


def archive_expired(max_age_days):
    """Move jobs older than max_age_days into the archive directory; return their ids."""
    if max_age_days is None:
        return []
    archive = JOBS_DIR / ARCHIVE_NAME
    moved = []
    for folder in sorted(p for p in JOBS_DIR.iterdir() if p.is_dir() and p.name != ARCHIVE_NAME) if JOBS_DIR.exists() else []:
        with _lock:
            process = _processes.get(folder.name)
        if process is not None and process.poll() is None:
            continue
        age = _age_days(folder)
        if age is not None and age > max_age_days:
            archive.mkdir(parents=True, exist_ok=True)
            target = archive / folder.name
            if target.exists():
                target = archive / f"{folder.name}-{secrets.token_hex(3)}"
            shutil.move(str(folder), str(target))
            moved.append(folder.name)
    return moved


def analysis_settings(form):
    """Translate the analysis form fields into init keys; invalid values raise ValueError."""
    init = {}

    def integer(name, minimum):
        value = form.get(name)
        if value in (None, ""):
            return None
        try:
            value = int(value)
        except ValueError as exc:
            raise ValueError(f"{name} must be an integer") from exc
        if value < minimum:
            raise ValueError(f"{name} must be at least {minimum}")
        return value

    random_starts = integer("multistart_random", 0)
    if random_starts:
        seed = integer("multistart_seed", 0)
        if seed is None:
            raise ValueError("multistart_seed is required with random restarts")
        init["multistart"] = {"random_starts": random_starts, "seed": seed}
    grid = (form.get("profile_kex") or "").strip()
    if grid:
        try:
            values = [float(v) for v in grid.replace(",", " ").split()]
        except ValueError as exc:
            raise ValueError("profile_kex must list numbers") from exc
        if not values or any(v <= 0 for v in values):
            raise ValueError("profile_kex values must be positive")
        init["profile"] = {"kex": values}
    interval = form.getlist("profile_interval") if hasattr(form, "getlist") else form.get("profile_interval", [])
    interval = [v for v in interval if v in ("kex", "pB", "v1n_scale")]
    if interval:
        init["profile_interval"] = {"parameters": interval}
    replicates = integer("bootstrap_replicates", 1)
    if replicates:
        seed = integer("bootstrap_seed", 0)
        if seed is None:
            raise ValueError("bootstrap_seed is required with bootstrap replicates")
        init["bootstrap"] = {"replicates": replicates, "seed": seed}
    return init


PAGE = r"""<!doctype html><html lang="en"><head><meta charset="utf-8"><title>SBONEST Sideband runner</title>
<style>body{font-family:sans-serif;margin:20px;max-width:1040px}fieldset{margin-bottom:14px}
pre{background:#f4f4f4;padding:8px;overflow:auto;max-height:320px}label{display:block;margin:4px 0}
.ok{color:#1a7f37}.bad{color:#b42318}img{max-width:100%;border:1px solid #ddd}</style></head><body>
<h1>SBONEST Sideband runner</h1>
<p>Upload a Sideband configuration and its data files. The configuration's <code>datasets</code>
entries are matched to the uploaded file names; the output prefix becomes <code>fit</code>.
Each job keeps its own checkpoint, so an interrupted fit can be resumed.{% if token_required %}
This server requires an access token.{% endif %}</p>
{% if token_required %}<label>Access token <input id="token" type="password"></label>{% endif %}
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
<fieldset><legend>3. Optional analyses (added to <code>init</code>)</legend>
<label>Random restarts <input type="number" name="multistart_random" min="0" placeholder="0"> seed <input type="number" name="multistart_seed" min="0"></label>
<label>kex profile grid (s⁻¹, space separated) <input type="text" name="profile_kex" placeholder="e.g. 250 300 350"></label>
<label>Profile intervals <label><input type="checkbox" name="profile_interval" value="kex"> kex</label>
<label><input type="checkbox" name="profile_interval" value="pB"> pB</label></label>
<label>Bootstrap replicates <input type="number" name="bootstrap_replicates" min="1" placeholder="0"> seed <input type="number" name="bootstrap_seed" min="0"></label>
</fieldset>
<button type="submit">Create job and run preflight check</button>
</form>
<pre id="action-error" class="bad" role="alert" hidden></pre>
<div id="job" hidden><h2>Job <span id="jobid"></span></h2>
<pre id="check"></pre>
<button id="fit">Start fit</button> <button id="resume">Resume</button> <button id="report">Make report</button>
<button id="archive">Archive job</button> <span id="state"></span><pre id="log"></pre><div id="files"></div>
<div id="previewbox" hidden><h3>Preview</h3><img id="preview" alt="fit preview"></div></div>
<h2>Jobs</h2><pre id="jobs"></pre>
<script>
let job=null;
const q=s=>document.querySelector(s);
function headers(){const t=q('#token');return t&&t.value?{'X-SBONEST-Token':t.value}:{};}
function actionError(message){q('#action-error').textContent=message;q('#action-error').hidden=!message;}
async function requestJSON(url,options={}){
 let r;
 try{
  r=await fetch(url,{...options,headers:headers()});
  const text=await r.text();let data;
  try{data=JSON.parse(text);}catch{ return {ok:false,status:r.status,error:'HTTP '+r.status+': invalid JSON response'+(text?' - '+text:'')}; }
  if(!r.ok)return {ok:false,status:r.status,data,error:'HTTP '+r.status+': '+(data&&data.error||r.statusText)};
  if(!data||typeof data!=='object')return {ok:false,status:r.status,error:'HTTP '+r.status+': invalid JSON response'};
  return {ok:true,status:r.status,data};
 }catch(e){return {ok:false,status:r?r.status:null,error:(r?'HTTP '+r.status+': ':'Network error: ')+e.message};}
}
async function post(url,body){return requestJSON(url,{method:'POST',body});}
q('#upload').addEventListener('submit',async e=>{e.preventDefault();
 const d=new FormData(e.target);const r=await post('/jobs',d);
 if(!r.ok){actionError(r.error);return;}const j=r.data;actionError('');
 job=j.job_id;q('#jobid').textContent=job;q('#job').hidden=false;
 q('#previewbox').hidden=true;q('#preview').removeAttribute('src');
 q('#check').textContent=JSON.stringify(j.check,null,2);poll();listJobs();});
for(const action of ['fit','resume','report','archive'])q('#'+action).addEventListener('click',async()=>{
 const selected=job;if(!selected)return;const r=await post('/jobs/'+selected+'/'+action);
 if(job!==selected)return;
 if(!r.ok){actionError(r.error);return;}actionError('');
 if(action==='archive'){job=null;q('#job').hidden=true;listJobs();}else{poll();}
});
async function listJobs(){const r=await requestJSON('/jobs');if(!r.ok){actionError(r.error);return;}q('#jobs').textContent=JSON.stringify(r.data,null,2);}
async function poll(){const selected=job;if(!selected)return;const r=await requestJSON('/jobs/'+selected);
 if(job!==selected)return;
 if(!r.ok){actionError(r.error);return;}const s=r.data;
 q('#state').textContent=s.running?'running':(s.result?('finished: '+(s.result.success?'success':'failed')):'idle');
 q('#state').className=s.result&&s.result.success?'ok':(s.result?'bad':'');
 q('#log').textContent=(s.result?JSON.stringify(s.result,null,2)+'\n':'')+s.log_tail;
 const t=q('#token');const suffix=t&&t.value?('?token='+encodeURIComponent(t.value)):'';
 q('#files').innerHTML=s.outputs.map(n=>'<a href="/jobs/'+selected+'/files/'+n+suffix+'">'+n+'</a>').join(' | ');
 q('#previewbox').hidden=!(s.preview&&s.result);
 if(s.preview&&s.result){const url=new URL('/jobs/'+selected+'/preview.png',location.href);
  if(t&&t.value)url.searchParams.set('token',t.value);url.searchParams.set('t',Date.now());q('#preview').src=url.href;}
 if(s.running)setTimeout(()=>{if(job===selected)poll();},3000);}
listJobs();
</script></body></html>"""


@app.route("/", methods=["GET"])
def index():
    """Single-page runner interface."""
    return render_template_string(PAGE, max_workers=MAX_WORKERS, token_required=bool(app.config.get("SBONEST_TOKEN")))


@app.route("/jobs", methods=["POST"])
def create_job():
    """Store the upload, rewrite the configuration and run the preflight check."""
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
    try:
        analyses = analysis_settings(request.form)
    except ValueError as exc:
        return jsonify({"error": str(exc)}), 400
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
    config["init"].update(analyses)
    if waveform is not None and waveform.filename:
        wname = secure_filename(waveform.filename)
        waveform.save(folder / wname)
        sideband = config.get("sideband", {})
        for entry in [sideband.get("decoupling", {}), *sideband.get("datasets", [])]:
            if isinstance(entry, dict) and "waveform_json" in entry:
                entry["waveform_json"] = wname
    (folder / "config.json").write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    (folder / "job.json").write_text(json.dumps({"job_id": job_id, "workers": workers, "no_pdf": no_pdf,
                                                 "analyses": sorted(analyses),
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
    return jsonify({"job_id": job_id, "check": check, "analyses": sorted(analyses)})


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
    """Start the fit as a background run.py process."""
    folder = _job_dir(job_id)
    if (folder / "fit_checkpoint").exists():
        return jsonify({"error": "This job already has a checkpoint; use resume or create a new job"}), 409
    return _start(job_id, folder, resume=False)


@app.route("/jobs/<job_id>/resume", methods=["POST"])
def resume_fit(job_id):
    """Resume an interrupted or finished job through its checkpoint."""
    folder = _job_dir(job_id)
    if not (folder / "fit_checkpoint").exists():
        return jsonify({"error": "Nothing to resume; start a fit first"}), 409
    return _start(job_id, folder, resume=True)


@app.route("/jobs/<job_id>/report", methods=["POST"])
def make_report(job_id):
    """Regenerate the saved-result report with a fresh prefix."""
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


@app.route("/jobs/<job_id>/preview.png", methods=["GET"])
def preview(job_id):
    """PNG preview of the fitted profiles, rendered from the predictions CSV and cached."""
    folder = _job_dir(job_id)
    csv_path = folder / "fit_predictions.csv"
    png = folder / "preview.png"
    if not csv_path.exists():
        abort(404)
    if not png.exists() or png.stat().st_mtime < csv_path.stat().st_mtime:
        import csv

        from sb_report import PREDICTION_COLUMNS, preview_png

        with csv_path.open(newline="", encoding="utf-8") as stream:
            rows = []
            for raw in csv.DictReader(stream):
                row = {key: float(raw[key]) for key in PREDICTION_COLUMNS[2:]}
                row.update(residue=raw["residue"], dataset_index=int(raw["dataset_index"]))
                rows.append(row)
        pending = folder / f".preview-{secrets.token_hex(3)}.png"
        preview_png(pending, rows)
        os.replace(pending, png)
    return send_from_directory(folder, "preview.png", max_age=0)


@app.route("/jobs/<job_id>/archive", methods=["POST"])
def archive_job(job_id):
    """Move a finished job into SB_JOBS/archive/ (nothing is deleted)."""
    folder = _job_dir(job_id)
    with _lock:
        process = _processes.get(job_id)
        if process is not None and process.poll() is None:
            return jsonify({"error": "The job is still running"}), 409
        _processes.pop(job_id, None)
    archive = JOBS_DIR / ARCHIVE_NAME
    archive.mkdir(parents=True, exist_ok=True)
    target = archive / job_id
    if target.exists():
        target = archive / f"{job_id}-{secrets.token_hex(3)}"
    shutil.move(str(folder), str(target))
    return jsonify({"job_id": job_id, "archived_to": str(target)})


@app.route("/jobs/archive-expired", methods=["POST"])
def archive_expired_jobs():
    """Archive every idle job older than the configured maximum age."""
    moved = archive_expired(app.config.get("SBONEST_MAX_AGE_DAYS"))
    return jsonify({"archived": moved, "max_age_days": app.config.get("SBONEST_MAX_AGE_DAYS")})


@app.route("/jobs/<job_id>", methods=["GET"])
def job_status(job_id):
    """Process state, checkpoint records, result summary and log tail of a job."""
    folder = _job_dir(job_id)
    return jsonify(_status(job_id, folder))


@app.route("/jobs", methods=["GET"])
def list_jobs():
    """Names, ages and states of all job directories."""
    JOBS_DIR.mkdir(parents=True, exist_ok=True)
    jobs = []
    for folder in sorted(p for p in JOBS_DIR.iterdir() if p.is_dir() and p.name != ARCHIVE_NAME):
        with _lock:
            process = _processes.get(folder.name)
        jobs.append({"job_id": folder.name, "running": process is not None and process.poll() is None,
                     "finished": (folder / "fit_result.json").exists(),
                     "age_days": round(_age_days(folder) or 0.0, 3)})
    archived = sorted(p.name for p in (JOBS_DIR / ARCHIVE_NAME).iterdir()) if (JOBS_DIR / ARCHIVE_NAME).is_dir() else []
    return jsonify({"jobs": jobs, "archived": archived, "max_age_days": app.config.get("SBONEST_MAX_AGE_DAYS")})


@app.route("/jobs/<job_id>/files/<path:filename>", methods=["GET"])
def download(job_id, filename):
    """Serve one file from inside a job directory."""
    folder = _job_dir(job_id)
    safe = secure_filename(Path(filename).name)
    if not safe or safe != filename or not (folder / safe).is_file():
        abort(404)
    return send_from_directory(folder, safe, as_attachment=True)


def main(argv=None):
    """Start the Flask runner."""
    import argparse

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=5050)
    parser.add_argument("--token", help="Access token (default: SBONEST_TOKEN environment variable)")
    parser.add_argument("--max-age-days", type=float, help="Archive idle jobs older than this many days")
    args = parser.parse_args(argv)
    if args.token:
        app.config["SBONEST_TOKEN"] = args.token
    if args.max_age_days is not None:
        if args.max_age_days <= 0:
            parser.error("--max-age-days must be positive")
        app.config["SBONEST_MAX_AGE_DAYS"] = args.max_age_days
    JOBS_DIR.mkdir(parents=True, exist_ok=True)
    moved = archive_expired(app.config.get("SBONEST_MAX_AGE_DAYS"))
    print(f"Sideband runner at http://{args.host}:{args.port}/  (jobs in {JOBS_DIR}; "
          f"token {'required' if app.config.get('SBONEST_TOKEN') else 'not set'}; archived {len(moved)} expired jobs)")
    app.run(host=args.host, port=args.port, debug=False)


if __name__ == "__main__":
    main()
