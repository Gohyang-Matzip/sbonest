#!/usr/bin/env python
import os
import sys
import uuid
import json
import subprocess
from flask import (
    Flask,
    request,
    jsonify,
    send_from_directory,
    render_template_string,
    url_for,
)
from werkzeug.utils import secure_filename
from estmodel import VALID_METHODS

# --- Configuration ---
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
JOBS_DIR = os.path.join(BASE_DIR, "JOB_FILES")
ALLOWED_EXTENSIONS = {"txt"}

app = Flask(__name__)
app.config["JOBS_DIR"] = JOBS_DIR
app.config["MAX_CONTENT_LENGTH"] = 32 * 1024 * 1024


@app.before_request
def validate_job_id():
    job_id = (request.view_args or {}).get("job_id")
    if job_id is not None and (job_id in (".", "..") or "\\" in job_id):
        return jsonify({"success": False, "error": "Invalid job ID."}), 400


# --- Helper Functions ---
def allowed_file(filename):
    return "." in filename and filename.rsplit(".", 1)[1].lower() in ALLOWED_EXTENSIONS


def positive_int(value):
    if isinstance(value, bool) or not isinstance(value, (int, str)):
        raise ValueError("Expected a positive integer.")
    result = int(value)
    if result <= 0:
        raise ValueError("Expected a positive integer.")
    return result


def run_script(command_args, working_dir, script_name_for_log="script"):
    log_prefix = f"[{script_name_for_log}@{working_dir}]"
    print(f"{log_prefix} Executing: {' '.join(command_args)}")
    try:
        process = subprocess.run(
            [sys.executable] + command_args,
            cwd=working_dir,
            capture_output=True,
            text=True,
            check=False,
            env=dict(os.environ, PYTHONIOENCODING="utf-8"),
        )
        stdout_output = process.stdout.strip() if process.stdout else ""
        stderr_output = process.stderr.strip() if process.stderr else ""
        print(f"{log_prefix} STDOUT:\n{stdout_output}")
        if stderr_output:
            print(f"{log_prefix} STDERR:\n{stderr_output}")
        print(f"{log_prefix} Return Code: {process.returncode}")
        return {
            "success": process.returncode == 0,
            "stdout": stdout_output,
            "stderr": stderr_output,
            "returncode": process.returncode,
        }
    except FileNotFoundError:
        error_msg = f"{log_prefix} Error: Script '{command_args[0]}' not found."
        print(error_msg)
        return {"success": False, "stdout": "", "stderr": error_msg, "returncode": -1}
    except Exception as e:
        error_msg = f"{log_prefix} Exception: {str(e)}"
        print(error_msg)
        return {"success": False, "stdout": "", "stderr": error_msg, "returncode": -1}


def get_project_name_from_config(config_path):
    try:
        with open(config_path, "r") as f:
            config_data = json.load(f)
        return config_data.get("Project Name")
    except Exception as e:
        print(f"Error reading project name from {config_path}: {e}")
        return None


def run_fit_script(job_id, job_dir, script_name, output_suffixes, extra_args=()):
    """Run a fit or MC script and return links to its existing output files."""
    project_name = get_project_name_from_config(os.path.join(job_dir, "config.json"))
    if not project_name:
        return jsonify(
            {"success": False, "error": "Could not read Project Name from config.json."}
        ), 500

    cmd_args = [os.path.join(BASE_DIR, script_name), "config.json", *extra_args]
    result = run_script(cmd_args, job_dir, script_name)
    output_files = []
    if result["success"]:
        for suffix in output_suffixes:
            filename = project_name + suffix
            if os.path.exists(os.path.join(job_dir, filename)):
                output_files.append(
                    {
                        "name": filename,
                        "url": url_for(
                            "download_file", job_id=job_id, filename=filename
                        ),
                    }
                )
    return jsonify(
        {
            "success": result["success"],
            "job_id": job_id,
            "message": f"{script_name} execution {'completed' if result['success'] else 'failed'}.",
            "script_stdout": result["stdout"],
            "script_stderr": result["stderr"],
            "output_files": output_files,
        }
    )


# --- Flask Routes ---
@app.route("/", methods=["GET"])
def index():
    # HTML/JavaScript 부분에 Project Name 입력 필드 추가
    return render_template_string("""
    <!doctype html>
    <html lang="en">
    <head>
        <meta charset="utf-8">
        <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
        <title>ONEST Model Runner</title>
        <style>
            body { font-family: sans-serif; margin: 20px; background-color: #f4f4f4; color: #333; line-height: 1.6; }
            .container { background-color: #fff; padding: 20px; border-radius: 8px; box-shadow: 0 0 10px rgba(0,0,0,0.1); max-width: 800px; margin: auto;}
            h1, h2, h3 { color: #333; border-bottom: 1px solid #eee; padding-bottom: 10px;}
            label { display: block; margin-top: 10px; font-weight: bold; }
            input[type=file], input[type=number], input[type=text], select {
                width: calc(100% - 24px); padding: 10px; margin-top: 5px; border-radius: 4px; border: 1px solid #ddd; box-sizing: border-box;
            }
            button {
                background-color: #007bff; color: white; padding: 10px 15px; border: none;
                border-radius: 4px; cursor: pointer; font-size: 16px; margin-top: 20px;
            }
            button:hover { background-color: #0056b3; }
            .output { margin-top: 20px; padding: 15px; background-color: #e9ecef; border: 1px solid #dee2e6; border-radius: 4px; white-space: pre-wrap; word-wrap: break-word; max-height: 400px; overflow-y: auto; }
            .file-list a { display: block; margin: 5px 0; color: #007bff; text-decoration: none; }
            .file-list a:hover { text-decoration: underline; }
            .error { color: #dc3545; font-weight: bold; }
            .stdout { color: #28a745; }
            .stderr { color: #ffc107; }
            .residue-item { margin-bottom: 8px; padding: 5px; border: 1px solid #eee; border-radius: 3px; }
            .residue-item label { display: inline-block; margin-right: 10px; font-weight: normal;}
            .residue-item input[type=checkbox] { vertical-align: middle; }
            hr { margin: 30px 0; border: 0; border-top: 1px solid #ccc;}
        </style>
    </head>
    <body>
        <div class="container">
            <h1>ONEST Web Interface</h1>

            <section id="initialPrepareSection">
                <h2>1. Initial Preparation (prepare.py)</h2>
                <form id="prepareForm" method="post" enctype="multipart/form-data">
                    <label for="data_files">Select Data Files
                        (e.g.,
                            <a href="{{ url_for('download_example', filename='syn10.txt') }}" download>syn10.txt</a>,
                            <a href="{{ url_for('download_example', filename='syn100.txt') }}" download>syn100.txt</a>
                        ):
                    </label>
                    <input type="file" name="data_files[]" id="data_files_input" multiple required>
                    <button type="button" onclick="submitInitialPrepare()">Run Initial Prepare</button>
                </form>
                <div id="initialPrepareResult" class="output" style="display:none;"></div>
            </section>
            <hr>

            <section id="configureSection" style="display:none;">
                <h2>2. Configure Parameters & Finalize</h2>
                <div id="configureFormContent"></div>
                <div id="finalizeResult" class="output" style="display:none;"></div>
            </section>
            <hr>

            <section id="executionSection" style="display:none;">
                <h2>3. Execute Analysis</h2>
                <label for="job_id_run">Job ID (from previous steps):</label>
                <input type="text" id="job_id_run" name="job_id_run" readonly>

                <div>
                    <h3>Run Standard Fit (run.py)</h3>
                    <button type="button" onclick="submitRunFit()">Run Standard Fit</button>
                </div>
                <div id="runFitResult" class="output" style="display:none;"></div>

                <hr>
                <div>
                    <h3>Run Monte Carlo Simulation (mcrun.py)</h3>
                    <label for="num_runs">Number of MC Runs:</label>
                    <input type="number" id="num_runs" name="num_runs" value="100" required>
                    <label for="num_processes">Number of Processes (optional, 0 or blank for default):</label>
                    <input type="number" id="num_processes" name="num_processes" placeholder="Defaults to CPU cores">
                    <button type="button" onclick="submitMcRun()">Run Monte Carlo</button>
                </div>
                <div id="mcRunResult" class="output" style="display:none;"></div>
            </section>
        </div>

        <script>
            async function submitInitialPrepare() {
                const form = document.getElementById('prepareForm');
                const formData = new FormData(form);
                const resultDiv = document.getElementById('initialPrepareResult');
                resultDiv.style.display = 'block';
                resultDiv.innerHTML = 'Processing initial preparation...';
                document.getElementById('configureSection').style.display = 'none';
                document.getElementById('executionSection').style.display = 'none';

                try {
                    const response = await fetch('/prepare_initial', { method: 'POST', body: formData });
                    const result = await response.json();

                    let html = '<h3>Initial Prepare Results:</h3>';
                    html += '<p><strong>Status:</strong> ' + (result.success ? '<span class="stdout">Success</span>' : '<span class="error">Failed</span>') + '</p>';
                    if (result.job_id) {
                        html += '<p><strong>Job ID:</strong> ' + result.job_id + '</p>';
                        if (result.initial_config_json) {
                             populateConfigureForm(result.initial_config_json, result.job_id);
                             document.getElementById('configureSection').style.display = 'block';
                             html += '<p>Scroll down to configure parameters and finalize.</p>';
                        }
                    }
                    if (result.script_stdout) html += '<h4>Script STDOUT:</h4><pre class="stdout">' + result.script_stdout + '</pre>';
                    if (result.script_stderr) html += '<h4>Script STDERR:</h4><pre class="stderr">' + result.script_stderr + '</pre>';
                    if (result.message) html += '<p><strong>Message:</strong> ' + result.message + '</p>';
                    if (result.error) html += '<p class="error"><strong>Error:</strong> ' + result.error + '</p>';
                    resultDiv.innerHTML = html;

                } catch (e) {
                    resultDiv.innerHTML = '<p class="error">An error occurred: ' + e.message + '</p>';
                }
            }

            function populateConfigureForm(initialConfig, jobId) {
                const configureFormDiv = document.getElementById('configureFormContent');
                let formHtml = '<h3>Configure Parameters</h3>';
                formHtml += `<input type="hidden" id="config_job_id" value="${jobId}">`;

                // Project Name Input
                const initialProjectName = (initialConfig && initialConfig['Project Name']) ? initialConfig['Project Name'] : 'default';
                formHtml += '<label for="project_name_input">Project Name:</label>';
                formHtml += `<input type="text" id="project_name_input" value="${initialProjectName}"><br/><br/>`;

                // Method Selection
                const methods = ['Baldwin', 'Matrix', 'NoEx'];
                let currentMethod = (initialConfig.init && initialConfig.init.Method) ? initialConfig.init.Method : 'Baldwin';
                formHtml += '<label for="selected_method_input">Select Method:</label>';
                formHtml += '<select id="selected_method_input">';
                methods.forEach(method => {
                    formHtml += `<option value="${method}" ${method === currentMethod ? 'selected' : ''}>${method}</option>`;
                });
                formHtml += '</select><br/><br/>';

                formHtml += '<h4>Residue Flags (check to turn ON):</h4>';
                if (initialConfig.residues && initialConfig.residues.length > 0) {
                    initialConfig.residues.forEach(residue => {
                        const residueName = residue.name;
                        const isChecked = residue.flag === 'on';
                        formHtml += `<div class="residue-item">
                                        <input type="checkbox" id="res_flag_${residueName}" name="${residueName}" data-residue-name="${residueName}" ${isChecked ? 'checked' : ''}>
                                        <label for="res_flag_${residueName}">${residueName}</label>
                                     </div>`;
                    });
                } else {
                    formHtml += '<p>No residues found in initial configuration.</p>';
                }
                formHtml += '<button type="button" onclick="submitFinalizeConfig()">Save Final Configuration</button>';
                configureFormDiv.innerHTML = formHtml;
            }

            async function submitFinalizeConfig() {
                const jobId = document.getElementById('config_job_id').value;
                const projectName = document.getElementById('project_name_input').value; // Get Project Name
                const selectedMethod = document.getElementById('selected_method_input').value;
                const resultDiv = document.getElementById('finalizeResult');
                resultDiv.style.display = 'block';
                resultDiv.innerHTML = 'Finalizing configuration...';

                if (!projectName.trim()) {
                    resultDiv.innerHTML = '<p class="error">Project Name cannot be empty.</p>';
                    return;
                }

                const residueFlags = [];
                const checkboxes = document.querySelectorAll('#configureFormContent input[type="checkbox"][data-residue-name]');
                checkboxes.forEach(cb => {
                    residueFlags.push({
                        name: cb.dataset.residueName,
                        flag: cb.checked ? 'on' : 'off'
                    });
                });

                const payload = {
                    project_name: projectName, // Add to payload
                    method: selectedMethod,
                    residue_flags: residueFlags
                };

                try {
                    const response = await fetch('/finalize_config/' + jobId, {
                        method: 'POST',
                        headers: { 'Content-Type': 'application/json' },
                        body: JSON.stringify(payload)
                    });
                    const result = await response.json();
                    let html = '<h3>Finalize Configuration Results:</h3>';
                    html += '<p><strong>Status:</strong> ' + (result.success ? '<span class="stdout">Success</span>' : '<span class="error">Failed</span>') + '</p>';
                    if (result.job_id) {
                        html += '<p><strong>Job ID:</strong> ' + result.job_id + '</p>';
                        document.getElementById('job_id_run').value = result.job_id;
                        if(result.success) document.getElementById('executionSection').style.display = 'block';
                    }
                    if (result.final_config_brief) {
                        html += '<p><strong>Final config.json (Brief - Project, Method, Residues):</strong></p>';
                        html += '<pre>' + JSON.stringify(result.final_config_brief, null, 2) + '</pre>';
                         if(result.download_url) html += '<p><a href="' + result.download_url + '" target="_blank">Download full final config.json</a></p>';
                    }
                    if (result.message) html += '<p><strong>Message:</strong> ' + result.message + '</p>';
                    if (result.error) html += '<p class="error"><strong>Error:</strong> ' + result.error + '</p>';
                    resultDiv.innerHTML = html;
                } catch (e) {
                    resultDiv.innerHTML = '<p class="error">An error occurred: ' + e.message + '</p>';
                }
            }

            // submitRunFit, submitMcRun, displayExecutionResult 함수는 이전과 동일하게 유지
            async function submitRunFit() {
                const jobId = document.getElementById('job_id_run').value;
                const resultDiv = document.getElementById('runFitResult');
                resultDiv.style.display = 'block';
                resultDiv.innerHTML = 'Processing Standard Fit...';
                if (!jobId) { resultDiv.innerHTML = '<p class="error">Please complete previous steps to get a Job ID.</p>'; return; }
                try {
                    const response = await fetch('/run_fit/' + jobId, { method: 'POST' });
                    const result = await response.json();
                    displayExecutionResult(result, resultDiv, "Standard Fit");
                } catch (e) { resultDiv.innerHTML = '<p class="error">An error occurred: ' + e.message + '</p>';}
            }

            async function submitMcRun() {
                const jobId = document.getElementById('job_id_run').value;
                const numRuns = document.getElementById('num_runs').value;
                const numProcesses = document.getElementById('num_processes').value;
                const resultDiv = document.getElementById('mcRunResult');
                resultDiv.style.display = 'block';
                resultDiv.innerHTML = 'Processing Monte Carlo Simulation...';
                if (!jobId) { resultDiv.innerHTML = '<p class="error">Please complete previous steps to get a Job ID.</p>'; return; }

                const payload = { num_runs: parseInt(numRuns) };
                if (numProcesses && parseInt(numProcesses) > 0) { payload.num_processes = parseInt(numProcesses); }

                try {
                    const response = await fetch('/run_mc/' + jobId, {
                        method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(payload)
                    });
                    const result = await response.json();
                    displayExecutionResult(result, resultDiv, "Monte Carlo");
                } catch (e) { resultDiv.innerHTML = '<p class="error">An error occurred: ' + e.message + '</p>';}
            }

            function displayExecutionResult(result, resultDiv, processName) {
                let html = '<h3>' + processName + ' Results:</h3>';
                html += '<p><strong>Status:</strong> ' + (result.success ? '<span class="stdout">Success</span>' : '<span class="error">Failed</span>') + '</p>';
                if (result.message) html += '<p><strong>Message:</strong> ' + result.message + '</p>';
                if (result.script_stdout) html += '<h4>Script STDOUT:</h4><pre class="stdout">' + result.script_stdout + '</pre>';
                if (result.script_stderr && result.script_stderr.trim() !== "") html += '<h4>Script STDERR:</h4><pre class="stderr">' + result.script_stderr + '</pre>';
                if (result.error) html += '<p class="error"><strong>Error:</strong> ' + result.error + '</p>';
                if (result.output_files && result.output_files.length > 0) {
                    html += '<h4>Output Files:</h4><div class="file-list">';
                    result.output_files.forEach(file => {
                        html += '<a href="' + file.url + '" target="_blank">' + file.name + '</a>';
                    });
                    html += '</div>';
                }
                resultDiv.innerHTML = html;
            }
        </script>
    </body>
    </html>
    """)


@app.route("/prepare_initial", methods=["POST"])
def prepare_initial_route():
    if "data_files[]" not in request.files:
        return jsonify({"success": False, "error": "No data files part"}), 400

    files = request.files.getlist("data_files[]")
    if not files or files[0].filename == "":
        return jsonify({"success": False, "error": "No selected files"}), 400

    saved_filenames = [secure_filename(f.filename) for f in files]
    if any(
        not allowed_file(f.filename) or not name
        for f, name in zip(files, saved_filenames)
    ):
        return jsonify({"success": False, "error": "Upload valid .txt files."}), 400
    if len(set(saved_filenames)) != len(saved_filenames):
        return jsonify(
            {"success": False, "error": "Uploaded filenames must be unique."}
        ), 400

    job_id = str(uuid.uuid4())
    job_dir = os.path.join(app.config["JOBS_DIR"], job_id)
    os.makedirs(job_dir, exist_ok=True)

    for file_storage, filename in zip(files, saved_filenames):
        file_storage.save(os.path.join(job_dir, filename))

    script_path = os.path.join(BASE_DIR, "prepare.py")
    cmd_args = [script_path] + saved_filenames
    script_result = run_script(cmd_args, job_dir, "prepare.py")

    if not script_result["success"] or not script_result["stdout"]:
        return jsonify(
            {
                "success": False,
                "error": "prepare.py execution failed or produced no JSON output.",
                "script_stdout": script_result["stdout"],
                "script_stderr": script_result["stderr"],
            }
        ), 500

    initial_config_json_str = script_result["stdout"]
    initial_config_path = os.path.join(job_dir, "initial_config.json")
    try:
        initial_config_data = json.loads(initial_config_json_str)
        with open(initial_config_path, "w") as f:
            json.dump(initial_config_data, f, indent=4)
    except json.JSONDecodeError:
        return jsonify(
            {
                "success": False,
                "error": "prepare.py output was not valid JSON.",
                "script_stdout": initial_config_json_str,
                "script_stderr": script_result["stderr"],
            }
        ), 500

    return jsonify(
        {
            "success": True,
            "job_id": job_id,
            "message": "Initial configuration generated. Please proceed to configure parameters.",
            "initial_config_json": initial_config_data,
            "script_stdout": "",
            "script_stderr": script_result["stderr"],
        }
    )


@app.route("/finalize_config/<job_id>", methods=["POST"])
def finalize_config_route(job_id):
    job_dir = os.path.join(app.config["JOBS_DIR"], job_id)
    initial_config_path = os.path.join(job_dir, "initial_config.json")

    if not os.path.exists(initial_config_path):
        return jsonify(
            {"success": False, "error": "Initial config not found for this job ID."}
        ), 404

    try:
        data = request.get_json(silent=True)
        if not isinstance(data, dict):
            return jsonify({"success": False, "error": "Expected a JSON object."}), 400
        selected_project_name = data.get("project_name")
        selected_method = data.get("method")
        residue_flags_from_client = data.get("residue_flags")

        if (
            not isinstance(selected_project_name, str)
            or not selected_project_name.strip()
        ):
            return jsonify(
                {"success": False, "error": "Project Name cannot be empty."}
            ), 400
        if selected_method not in VALID_METHODS:
            return jsonify(
                {"success": False, "error": "Unknown calculation method."}
            ), 400
        if not isinstance(residue_flags_from_client, list) or not all(
            isinstance(r, dict)
            and isinstance(r.get("name"), str)
            and r.get("flag") in ("on", "off")
            for r in residue_flags_from_client
        ):
            return jsonify(
                {
                    "success": False,
                    "error": "Residues must have a name and an on/off flag.",
                }
            ), 400

        safe_project_name = "".join(
            c if c.isalnum() or c in ("_", "-") else "_" for c in selected_project_name
        )
        with open(initial_config_path, "r") as f:
            config_data = json.load(f)

        config_data["Project Name"] = safe_project_name
        config_data.setdefault("init", {})["Method"] = selected_method

        client_flags_map = {
            res_info["name"]: res_info["flag"] for res_info in residue_flags_from_client
        }
        known_names = {r["name"] for r in config_data["residues"]}
        if (
            len(client_flags_map) != len(residue_flags_from_client)
            or not client_flags_map.keys() <= known_names
        ):
            return jsonify(
                {"success": False, "error": "Duplicate or unknown residue name."}
            ), 400

        for residue in config_data["residues"]:
            if residue["name"] in client_flags_map:
                residue["flag"] = client_flags_map[residue["name"]]

        if not any(r["flag"] == "on" for r in config_data["residues"]):
            return jsonify(
                {"success": False, "error": "Select at least one residue."}
            ), 400
        final_config_path = os.path.join(job_dir, "config.json")
        with open(final_config_path, "w") as f:
            json.dump(config_data, f, indent=4)

        final_config_brief = {
            "Project Name": config_data.get("Project Name"),
            "init": {"Method": config_data.get("init", {}).get("Method")},
            "residues": config_data.get("residues"),
        }

        return jsonify(
            {
                "success": True,
                "job_id": job_id,
                "message": "Configuration finalized and 'config.json' saved.",
                "final_config_brief": final_config_brief,
                "download_url": url_for(
                    "download_file", job_id=job_id, filename="config.json"
                ),
            }
        )

    except Exception as e:
        return jsonify(
            {"success": False, "error": f"Error finalizing configuration: {str(e)}"}
        ), 500


@app.route("/run_fit/<job_id>", methods=["POST"])
def run_fit_route(job_id):
    job_dir = os.path.join(app.config["JOBS_DIR"], job_id)
    config_path = os.path.join(job_dir, "config.json")

    if not os.path.isdir(job_dir) or not os.path.exists(config_path):
        return jsonify(
            {"success": False, "error": "Job ID or final config.json not found."}
        ), 404

    return run_fit_script(
        job_id, job_dir, "run.py", ("_data.pdf", "_result.txt", ".pdf")
    )


@app.route("/run_mc/<job_id>", methods=["POST"])
def run_mc_route(job_id):
    job_dir = os.path.join(app.config["JOBS_DIR"], job_id)
    config_path = os.path.join(job_dir, "config.json")

    if not os.path.isdir(job_dir) or not os.path.exists(config_path):
        return jsonify(
            {"success": False, "error": "Job ID or final config.json not found."}
        ), 404

    data = request.get_json(silent=True)
    if not isinstance(data, dict) or "num_runs" not in data:
        return jsonify(
            {"success": False, "error": "Missing 'num_runs' in request."}
        ), 400

    try:
        num_runs = positive_int(data["num_runs"])
    except ValueError:
        return jsonify(
            {"success": False, "error": "'num_runs' must be a positive integer."}
        ), 400

    num_processes = data.get("num_processes")
    try:
        if num_processes is not None and str(num_processes).strip() != "":
            num_processes = positive_int(num_processes)
        else:
            num_processes = None
    except ValueError:
        return jsonify(
            {"success": False, "error": "'num_processes' must be a positive integer."}
        ), 400

    extra_args = [str(num_runs)]
    if num_processes is not None:
        extra_args.append(str(num_processes))

    return run_fit_script(
        job_id,
        job_dir,
        "mcrun.py",
        ("_data.pdf", ".pdf", "_mc.txt", "_mcmean.pdf"),
        extra_args,
    )


@app.route("/download/<job_id>/<path:filename>")
def download_file(job_id, filename):
    job_dir = os.path.join(app.config["JOBS_DIR"], job_id)
    return send_from_directory(job_dir, filename, as_attachment=True)


@app.route("/example/<path:filename>")
def download_example(filename):
    example_dir = os.path.join(BASE_DIR, "example")
    safe_filename = secure_filename(filename)
    return send_from_directory(example_dir, safe_filename, as_attachment=True)


def create_app():  # Gunicorn·uWSGI용 팩토리
    return app  # 이미 만든 app 객체 그대로 반환


if __name__ == "__main__":
    os.makedirs(JOBS_DIR, exist_ok=True)
    print(f"Server running. Base directory for scripts: {BASE_DIR}")
    print(f"Job files will be stored under: {JOBS_DIR}")
    print("Access the web interface at http://127.0.0.1:5001")
    app.run(port=5001)
