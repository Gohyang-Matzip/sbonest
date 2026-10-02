"""Refit the existing full 1.2 GHz profiles through the CLI, without H-rate inputs.

Run: .venv/bin/python refit_automatic_proton.py --out results/auto_H_refit/fits
"""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import copy
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parent
PEAKS = ("A1", "G2", "S3")
PAIRS = {"two_RF": [0, 2], "three_RF": [0, 1, 2]}


def run_job(path):
    began = time.monotonic()
    env = dict(
        os.environ,
        OMP_NUM_THREADS="1",
        OPENBLAS_NUM_THREADS="1",
        VECLIB_MAXIMUM_THREADS="1",
        MPLBACKEND="Agg",
    )
    with path.with_suffix(".log").open("w") as log:
        proc = subprocess.run(
            [sys.executable, str(ROOT / "run.py"), str(path), "--no-pdf"],
            cwd=path.parent,
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
    if proc.returncode:
        raise RuntimeError(f"CLI failed; inspect {path.with_suffix('.log')}")
    result = json.loads(path.with_name(path.stem + "_result.json").read_text())
    return path.stem, result, time.monotonic() - began


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    source_path = ROOT / "results/peakwise_H_fit/summary.json"
    source = json.loads(source_path.read_text())
    # Archived metadata records the original machine; use this checkout's data.
    source["input_files"] = {
        kind: [str(source_path.parent / Path(p).name) for p in paths]
        for kind, paths in source["input_files"].items()
    }
    old_meta = json.loads((ROOT / "results/two_RF_check/metadata.json").read_text())
    truth = old_meta["truth"] | {"kex": 300.0, "pB": 0.05}
    rng = np.random.default_rng(20261005)
    starts = [
        {
            f"{peak}.{rate}": float(np.exp(rng.uniform(np.log(lo), np.log(hi))))
            for peak in PEAKS
            for rate, lo, hi in (("R1H", 0.2, 20.0), ("R2H", 5.0, 100.0))
        }
        for _ in range(3)
    ]
    paths = []
    for label, indices in PAIRS.items():
        for kind, suffix, start in [("noisy", "", {}), ("exact", "_exact", {})] + [
            ("noisy", f"_start_{i + 1:02}", s) for i, s in enumerate(starts)
        ]:
            name = label + suffix
            cfg = copy.deepcopy(source["config"])
            cfg["Project Name"] = str(out / name)
            cfg["datasets"] = [source["input_files"][kind][i] for i in indices]
            for key in ("R1H", "R2H"):
                cfg["sideband"]["decoupling"].pop(key, None)
            cfg["init"]["initial"].update(start)
            path = out / (name + ".json")
            path.write_text(json.dumps(cfg, indent=2) + "\n")
            paths.append(path)
    files = [
        Path(__file__),
        ROOT / "sbfit.py",
        ROOT / "sideband.py",
        ROOT / "fit.py",
        ROOT / "run.py",
        ROOT / "estmodel.py",
        source_path,
    ] + [Path(p) for paths_ in source["input_files"].values() for p in paths_]
    meta = dict(
        truth=truth,
        random_start_seed=20261005,
        random_H_starts=starts,
        H_mode="Peakwise fitted nuisance parameters; shared across RF and states A/B",
        default_H_start={"R1H": 2.0, "R2H": 25.0},
        source_sha256={
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in files
        },
    )
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    results = {}
    with ThreadPoolExecutor(max_workers=2) as pool:
        for future in as_completed([pool.submit(run_job, path) for path in paths]):
            name, r, elapsed = future.result()
            results[name] = r
            print(
                name,
                "chi2",
                round(r["chi2"], 6),
                "kex",
                round(r["kex"], 5),
                "RF",
                round(r["parameters"]["v1n_scale"]["value"], 7),
                "bounds",
                r["at_bounds"],
                "seconds",
                round(elapsed, 1),
                flush=True,
            )
    for label, indices in PAIRS.items():
        for name, r in results.items():
            if not name.startswith(label):
                continue
            assert r["n_points"] == 147 * 3 * len(indices) and r["n_parameters"] == 24
            assert r["proton_relaxation_mode"] == "fit" and r["jacobian_rank"] == 24
            if name.endswith("_exact"):
                assert r["chi2"] < 1e-8, (name, r["chi2"])
                for n, v in r["parameters"].items():
                    assert abs(v["value"] - truth[n]) < 0.003, (name, n, v)
            if "_start_" in name:
                assert abs(r["chi2"] - results[label]["chi2"]) < 0.02, name
    (out / "summary.json").write_text(
        json.dumps(results, indent=2, allow_nan=False) + "\n"
    )
    summarize(out, results, truth)
    print(
        "PASS: 10 CLI refits, full profiles, exact recovery and three random starts per RF design",
        flush=True,
    )


def summarize(out, results, truth):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sbfit import SidebandModel

    rows = []
    fig, axes = plt.subplots(
        3, 2, figsize=(11, 8), sharex=True, sharey=True, layout="constrained"
    )
    colors = {25.0: "#0072B2", 50.0: "#D55E00", 100.0: "#009E73"}
    for col, label in enumerate(PAIRS):
        result = results[label]
        model = SidebandModel(result["config"])
        p = np.array([result["parameters"][n]["value"] for n in model.parameter_names])
        P = model.seParam(p)
        for i, residue in enumerate(model.dataset.res):
            ax = axes[i, col]
            for es in residue.estSpecs:
                x = np.asarray(es.offset)
                observed = np.asarray(es.int)
                predicted = model.calc(P, i, x, es)
                fine = np.linspace(105, 135, 3001)
                ax.plot(
                    fine,
                    model.calc(P, i, fine, es),
                    color=colors[es.v1],
                    lw=1.1,
                    label=f"{es.v1:g} Hz",
                )
                ax.plot(x, observed, ".", color=colors[es.v1], ms=2.5, alpha=0.65)
                for a, b, c, d in zip(x, observed, predicted, es.intstd):
                    rows.append(
                        dict(
                            design=label,
                            peak=residue.label,
                            nominal_RF_Hz=es.v1,
                            offset_ppm=a,
                            observed=b,
                            fitted=c,
                            sigma=d,
                            residual_sigma=(c - b) / d,
                        )
                    )
            ax.set_title(
                f"{residue.label} · {'25/100' if col == 0 else '25/50/100'} Hz",
                fontsize=10,
            )
            ax.set_xlim(135, 105)
            if col == 0:
                ax.set_ylabel("Normalized intensity")
            if i == 2:
                ax.set_xlabel("¹⁵N offset (ppm)")
            if i == 0:
                ax.legend(frameon=False, fontsize=8, loc="lower left")
    fig.suptitle(
        "1.2 GHz full CEST profiles including 90°x–240°y–90°x sidebands\n"
        "Dots: existing synthetic data; lines: refit with automatic peakwise proton relaxation",
        fontsize=12,
    )
    fig.savefig(out / "full_profile_fit.png", dpi=220)
    fig.savefig(out / "full_profile_fit.pdf")
    plt.close(fig)
    with (out / "predictions.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    for label in PAIRS:
        chi2 = sum(r["residual_sigma"] ** 2 for r in rows if r["design"] == label)
        np.testing.assert_allclose(chi2, results[label]["chi2"], atol=1e-6)
    lines = [
        "# R1H/R2H 입력 없이 전체 sideband profile 재fitting",
        "",
        "1.2 GHz·105–135 ppm·약 25 Hz 간격의 기존 3-peak **합성 데이터**를 그대로 사용했다.",
        "90°x–240°y–90°x decoupling sideband를 포함한 모든 점을 fitting했다. 새로운 실측 결과는 아니다.",
        "",
        "## 프로그램 동작",
        "",
        "사용자는 R1H/R2H를 입력하지 않아도 된다. 프로그램이 (2,25) s⁻¹에서 시작해",
        "각 peak의 두 값을 질소·교환·RF 파라미터와 함께 최적화한다. 임의의 값을 고정하는 방식이 아니다.",
        "H rates는 서로 다른 peak 사이에 공유하지 않으며, 같은 peak의 RF 파일과 A/B 상태 사이에서는 공유한다.",
        "24개 자유 파라미터의 covariance에 H-rate 불확도를 포함한다. 경계에서의 대칭 국소 오차는 정식 confidence interval이 아니다.",
        "",
        "## 기본 자동 fitting 결과",
        "",
        "| RF (Hz) | 관측점 | 자유 파라미터 | kex (s⁻¹) | pB | v1n_scale | Reduced χ² |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for label, indices in PAIRS.items():
        r = results[label]
        rf = "/".join(str([25, 50, 100][i]) for i in indices)
        lines.append(
            f"| {rf} | {r['n_points']} | {r['n_parameters']} | {r['kex']:.5f} | {r['pB']:.7f} | {r['parameters']['v1n_scale']['value']:.7f} | {r['chi2'] / r['dof']:.5f} |"
        )
    lines += [
        "",
        "생성값은 kex=300 s⁻¹, pB=0.05, v1n_scale=1.08이다.",
        "",
        "| 파라미터 | 생성값 | 25/100 Hz 추정값 ± SE | 25/50/100 Hz 추정값 ± SE |",
        "|---|---:|---:|---:|",
    ]
    for n in results["two_RF"]["parameters"]:
        if n in ("kab", "kba") or n.endswith((".R1H", ".R2H")):
            continue
        cells = []
        for label in PAIRS:
            v = results[label]["parameters"][n]
            cells.append(f"{v['value']:.8g} ± {v['stderr']:.4g}")
        lines.append(f"| {n} | {truth[n]:g} | " + " | ".join(cells) + " |")
    lines += [
        "",
        "![전체 profile 재fitting](full_profile_fit.png)",
        "",
        "## 초기값 재검사",
        "",
        "R1H는 0.2–20, R2H는 5–100 s⁻¹에서 peak별로 독립 log-uniform 난수를 뽑아",
        "세 번 재시작했다(seed 20261005). 두 RF 조합에 같은 초기값 세트를 사용했다.",
        "",
        "| 조합 | 기본 fit과의 최대 절댓값 Δχ² | kex 범위 (s⁻¹) | v1n_scale 범위 |",
        "|---|---:|---:|---:|",
    ]
    for label in PAIRS:
        rr = [results[label]] + [results[f"{label}_start_{i:02}"] for i in range(1, 4)]
        delta = max(abs(r["chi2"] - rr[0]["chi2"]) for r in rr)
        kval = [r["kex"] for r in rr]
        rf = [r["parameters"]["v1n_scale"]["value"] for r in rr]
        lines.append(
            f"| {label} | {delta:.6g} | {min(kval):.6f}–{max(kval):.6f} | {min(rf):.8f}–{max(rf):.8f} |"
        )
    lines += [
        "",
        "무잡음 대조군 두 개에서 전체 24개 생성값 복원도 확인했다.",
        "",
        "## 내부 H-rate 진단",
        "",
        "입력이 필요 없다는 것이 H rates가 중요하지 않거나 정확히 결정된다는 뜻은 아니다.",
        "개별 H rates가 약하게 결정되거나 하한에 도달하면 진단 파일에 기록한다.",
        "",
        "| 조합 | Peak | R1H ± SE (s⁻¹) | R2H ± SE (s⁻¹) |",
        "|---|---|---:|---:|",
    ]
    for label in PAIRS:
        r = results[label]
        for peak in PEAKS:
            values = [r["parameters"][f"{peak}.{key}"] for key in ("R1H", "R2H")]
            lines.append(
                f"| {label} | {peak} | "
                + " | ".join(f"{v['value']:.5g} ± {v['stderr']:.4g}" for v in values)
                + " |"
            )
    for label in PAIRS:
        lines += ["", f"{label} 진단: " + "; ".join(results[label]["warnings"]), ""]
    lines += [
        "## 실행 및 산출물",
        "",
        "기본 설정에는 R1H/R2H 입력이 없다. 동일 설정을 재실행하려면 Project Name을 새 출력 경로로 바꾼다.",
        "",
        "```bash",
        "cd /Users/donghanlee/work/projects/sbonest",
        ".venv/bin/python refit_automatic_proton.py --out results/auto_H_repeat_01",
        "```",
        "",
        "[두 RF 설정](two_RF.json) · [두 RF 결과](two_RF_result.json) · [세 RF 설정](three_RF.json) · [세 RF 결과](three_RF_result.json) ·",
        "[전체 결과](summary.json) · [측정점별 예측·잔차](predictions.csv) · [그림 PDF](full_profile_fit.pdf) · [입력·코드 hash](metadata.json)",
        "",
    ]
    (out / "REFIT_REPORT.md").write_text("\n".join(lines))


if __name__ == "__main__":
    main()
