"""Create descriptive sensitivity tables/plot from the completed fixed-H refits."""
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
meta = json.loads((HERE/"metadata.json").read_text())
results = json.loads((HERE/"results.json").read_text())
by_id = {r["id"]:r for r in results}
baseline = by_id["baseline_noisy"]
truth = meta["truth"] | {"kex":300., "pB":.05}
arms = ["near_both","broad_R1H","broad_R2H","broad_both"]
labels = ["Both ±20%", "R1H ×0.5–2", "R2H ×0.5–2", "Both ×0.5–2"]
groups = {
    "kex": ["kex"], "dw": [f"{p}.dw_ppm" for p in ("A1","G2","S3")],
    "v1n": ["v1n_scale"], "R1N": [f"{p}.R1" for p in ("A1","G2","S3")],
    "R2Na": [f"{p}.R2a" for p in ("A1","G2","S3")],
    "R2Nb": [f"{p}.R2b" for p in ("A1","G2","S3")],
}
targets = ["kex","pB","v1n_scale"] + [n for keys in groups.values() for n in keys if n not in ("kex","v1n_scale")]


def value(r,n):
    return r[n] if n in ("kex","pB") else r["parameters"][n]


stats = {}
for arm in arms:
    rows = [r for r in results if r["arm"]==arm]
    assert len(rows)==meta["draws_per_arm"] and all(r["success"] for r in rows)
    stats[arm] = {}
    for n in targets:
        values = np.array([value(r,n) for r in rows])
        changes = values-value(baseline,n)
        relative = 100*changes/abs(truth[n])
        stats[arm][n] = dict(median_abs_shift_pct=float(np.median(abs(relative))),
            max_abs_shift_pct=float(max(abs(relative))),
            mean_signed_shift_pct=float(np.mean(relative)),
            rms_shift_pct=float(np.sqrt(np.mean(relative**2))),
            max_abs_shift=float(max(abs(changes))),
            max_abs_shift_in_baseline_SE=float(max(abs(changes))/baseline["stderr"][n]),
            min_value=float(min(values)),max_value=float(max(values)),
            max_abs_error_from_truth_pct=float(max(abs(values-truth[n]))/abs(truth[n])*100))
    stats[arm]["fit_quality"] = dict(
        min_reduced_chi2=min(r["chi2"]/r["dof"] for r in rows),
        median_reduced_chi2=float(np.median([r["chi2"]/r["dof"] for r in rows])),
        max_reduced_chi2=max(r["chi2"]/r["dof"] for r in rows),
        boundary_fits=[r["id"] for r in rows if r["at_bounds"]])
(HERE/"statistics.json").write_text(json.dumps(stats,indent=2,allow_nan=False)+"\n")

matrix = np.array([[max(stats[a][n]["max_abs_shift_pct"] for n in keys) for a in arms]
                   for keys in groups.values()])
assert (matrix>0).all()
fig,ax=plt.subplots(figsize=(8.2,4.8),layout="constrained")
low=10**np.floor(np.log10(matrix.min()))
high=10**np.ceil(np.log10(matrix.max()))
im=ax.imshow(matrix,cmap="YlGnBu",norm=LogNorm(vmin=low,vmax=high),aspect="auto")
ax.set_xticks(range(4),labels=labels)
ax.set_yticks(range(6),labels=list(groups))
ax.set_title("Sensitivity to random, fixed proton relaxation rates",pad=15)
for i in range(6):
    for k in range(4):
        x=matrix[i,k]
        text=f"{x:.3f}%" if x<.01 else f"{x:.2f}%"
        ax.text(k,i,text,ha="center",va="center",color="white" if im.norm(x)>.65 else "black",fontsize=10)
ax.set_xlabel(f"{meta['draws_per_arm']} draws per condition; maximum across peaks for local parameters",labelpad=12)
fig.colorbar(im,ax=ax,label="Maximum |shift| (% of true value)\nLogarithmic color scale",pad=.03)
fig.savefig(HERE/"parameter_sensitivity.png",dpi=220)
fig.savefig(HERE/"parameter_sensitivity.pdf")
plt.close(fig)

lines=["# 무작위 peak별 R1H·R2H 고정값의 영향 검사", "",
"검사일: 2026-10-02. 1.2 GHz·30 ppm의 3-peak 합성 데이터에 대한 검사다.","",
"## 결과 요약", "",
f"무작위 H 고정값 {4*meta['draws_per_arm']}회와 대조·재확인을 포함한 총 {len(results)}회 fitting을 완료했다.",
"요청한 파라미터 모두에 수치적인 영향이 있었으며, 이번 범위에서는 R2H 가정의 영향이 더 컸다.",
f"R1H만 0.5–2배로 바꾼 경우 기록한 파라미터의 변화는 기준 국소 표준오차의 최대 {max(stats['broad_R1H'][n]['max_abs_shift_in_baseline_SE'] for n in targets):.2f}배였다.",
f"두 H rate를 ±20% 범위로 바꾼 경우에도 최대 {max(stats['near_both'][n]['max_abs_shift_in_baseline_SE'] for n in targets):.2f}배였다.",
f"그러나 R2H만 0.5–2배로 바꾸면 R2Nb는 최대 {matrix[5,2]:.2f}%, v1n은 최대 {matrix[2,2]:.2f}% 변했다.",
f"둘 다 0.5–2배로 바꾼 경우 v1n 변화는 기준 표준오차의 최대 {stats['broad_both']['v1n_scale']['max_abs_shift_in_baseline_SE']:.2f}배였으므로, 상대 변화가 1% 미만이라는 이유만으로 무시할 수 없다.",
"이 조건에서는 R1H 고정값의 영향이 작았다. R2H는 peak별 fitting 또는 고정값 민감도 검사를 유지하는 편이 타당하다.","",
"## 검사 정의", "",
"`R1H`, `R2H`에 무작위 **고정값**을 넣고 `kab`, `kba`, `v1n_scale` 및",
"peak별 `peak_ppm`, `dw_ppm`, `R1`, `R2a`, `R2b`를 다시 fitting했다.",
"무작위 초기값을 넣은 뒤 H rates도 자유롭게 fitting한 검사가 아니다.",
"같은 peak의 A/B 상태와 여러 RF 파일은 같은 H rates를 사용하며, 서로 다른",
"peak의 H rates는 독립적으로 뽑았다. 각 fit의 자유 파라미터는 18개다.","",
"측정 profile과 잡음 realization은 매번 동일하다. 따라서 올바른 H rates를",
"고정한 기준 fit과의 차이는 H rates 가정의 영향을 나타낸다. 생성값으로부터의",
"전체 오차와 기준 fit으로부터의 변화는 `statistics.json`에서 따로 기록했다.","",
"| 조건 | R1H | R2H | 반복 수 |","|---|---|---|---:|",
f"| Both ±20% | 참값 × U(0.8,1.2) | 참값 × U(0.8,1.2) | {meta['draws_per_arm']} |",
f"| R1H ×0.5–2 | 참값 × log-U(0.5,2) | 참값 고정 | {meta['draws_per_arm']} |",
f"| R2H ×0.5–2 | 참값 고정 | 참값 × log-U(0.5,2) | {meta['draws_per_arm']} |",
f"| Both ×0.5–2 | 참값 × log-U(0.5,2) | 참값 × log-U(0.5,2) | {meta['draws_per_arm']} |","",
"넓은 범위의 세 조건은 동일한 난수쌍을 사용한다. 난수 seed는 20261004다.",
"이는 지정한 민감도 시험 범위이며, 모든 단백질의 일반적인 실험 범위라는 뜻은 아니다.","",
"H rate 생성값 (s⁻¹): A1=(2,25), G2=(1.2,18), S3=(3,35).",
"각 profile은 105–135 ppm의 147점, 간격 24.985 Hz, RF 25/50/100 Hz,",
"T=0.4 s, intensity σ=0.001이다. 전체 관측점은 1323개다.",
"RR decoupling은 90°x–240°y–90°x 반복이며 p90=70 μs다.","",
"## 요청한 파라미터별 최대 변화", "",
"아래 값은 `100 × |무작위 H 고정 fit − 기준 fit| / |생성값|`이다.",
"각 조건에서 관찰한 최댓값이며, 잔기별 파라미터는 3개 peak 전체의 최댓값을 표시한다.",
"미지의 모든 H rates에 대한 보장 범위나 confidence interval이 아니다.","",
"| 파라미터 | Both ±20% | R1H ×0.5–2 | R2H ×0.5–2 | Both ×0.5–2 |",
"|---|---:|---:|---:|---:|"]
for label,vals in zip(groups,matrix):
    lines.append("| "+label+" | "+" | ".join(f"{v:.4f}%" for v in vals)+" |")
lines += ["", "`v1n`은 공통 RF scale을 fitting했다. 실제 25/50/100 Hz RF 값의 상대 변화도",
          "동일하다. `dw`는 부호를 유지한 ppm 파라미터이며 상대 변화의 분모는 |dw_true|다.","",
          "![무작위 H rate 민감도](results/random_H_rates/parameter_sensitivity.png)","",
          "## 넓은 범위에서 둘 다 무작위로 고정한 결과", "",
          "Median과 max는 기준 fit으로부터의 절댓값 변화다. 마지막 열은 동일한 변화를",
          "기준 fit의 국소 표준오차로 나눈 값이며, H-rate 불확도를 포함한 총 오차가 아니다.","",
          "| 파라미터 | 생성값 | 기준 fit | Median 변화 (%) | 최대 변화 (%) | 최대 변화 / 기준 SE |",
          "|---|---:|---:|---:|---:|---:|"]
for n in targets:
    s=stats['broad_both'][n]
    lines.append(f"| {n} | {truth[n]:.6g} | {value(baseline,n):.7g} | {s['median_abs_shift_pct']:.4f} | {s['max_abs_shift_pct']:.4f} | {s['max_abs_shift_in_baseline_SE']:.2f} |")
lines += ["", "## Fit quality", "",
          f"기준 fit: χ²={baseline['chi2']:.6f}, dof={baseline['dof']}, reduced χ²={baseline['chi2']/baseline['dof']:.6f}.","",
          "| 조건 | Reduced χ² min | median | max | 경계 도달 fit 수 |",
          "|---|---:|---:|---:|---:|"]
for arm,label in zip(arms,labels):
    s=stats[arm]['fit_quality']
    lines.append(f"| {label} | {s['min_reduced_chi2']:.4f} | {s['median_reduced_chi2']:.4f} | {s['max_reduced_chi2']:.4f} | {len(s['boundary_fits'])} |")
lines += ["", "## 잡음 없는 데이터 및 다른 초기값으로 재검사", "",
          "파라미터 변화가 큰 두 fit과 χ²가 가장 큰 fit을 선정했다(중복 제외).",
          "같은 H 고정값으로 무잡음 데이터를 fitting하고, 잡음 데이터에서는 다른 시작점으로",
          "재fitting했다. 각 결과는 `exact_check_*`와 `restart_check_*` JSON에 보존했다.","",
          "| 선택한 fit | 잡음 fit χ² | 다른 시작점 χ² | 무잡음 fit χ² |",
          "|---|---:|---:|---:|"]
for r in results:
    if r['arm']=='exact_check':
        sid=r['source_fit']
        lines.append(f"| {sid} | {by_id[sid]['chi2']:.5f} | {by_id['restart_check_'+sid]['chi2']:.5f} | {r['chi2']:.5f} |")
lines += ["", "## 범위와 재현", "",
          "- H rates를 참값으로 고정한 무잡음 fit은 모든 질소/RF 파라미터를 복원했다.",
          "- 무작위 입력 fitting의 수렴 실패와 경계 도달을 기록했다. 실패를 조용히 제외하지 않았다.",
          "- 기준 covariance는 중앙 차분과 dense Jacobian으로 계산했으며 reduced χ²로 재조정하지 않았다.",
          "- 한 시료 조건, 3개 합성 peak, 고정된 잡음 realization에 대한 민감도다. 실제 시료나 다른 H shift/자기장/T에 일반화할 수 없다.",
          "- 각 peak의 A/B 상태에 같은 H relaxation을 가정했다. 독립적인 상태별 relaxation은 검사하지 않았다.",
          "- 기존 production fitting CLI는 변경하지 않았다.","",
          "```bash", "cd /Users/donghanlee/work/projects/sbonest",
          ".venv/bin/python check_random_proton_rates.py --out results/random_H_repeat_01 --draws 16 --workers 3", "```","",
          "새 출력 폴더를 사용한다. 결과 파일은 덮어쓰지 않는다.","",
          "[검사 코드](check_random_proton_rates.py) · [전체 수치 CSV](results/random_H_rates/results.csv) ·",
          "[전체 fitting JSON](results/random_H_rates/results.json) · [통계 JSON](results/random_H_rates/statistics.json) ·",
          "[조건·source hash](results/random_H_rates/metadata.json) · [그림 PDF](results/random_H_rates/parameter_sensitivity.pdf) ·",
          "[보고서/그림 생성 코드](results/random_H_rates/summarize.py)", ""]
(ROOT/"RANDOM_PROTON_RELAXATION_CHECK.md").write_text("\n".join(lines))
print("Maximum percent shifts:")
for group,row in zip(groups,matrix): print(group, dict(zip(arms,row)))
print("Saved statistics, report, PNG, PDF")
