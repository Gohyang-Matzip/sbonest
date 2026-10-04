# SBONEST: step-by-step dummy 데이터 가이드

[English](DUMMY_GUIDE.md) · [전체 매뉴얼](SBONEST_MANUAL.ko.md) · [README](README.md)

처음 사용한다면 이 가이드부터 따라 한다. 저장소의 **합성 데이터**를 복사하고,
3개 peak를 fitting한 뒤 결과와 잔차를 확인한다. 선택 단계에서는 다른 초기값과
likelihood scan, 소규모 bootstrap도 실행한다. 실험 데이터나 인접한 OC 저장소는
필요 없다. 2026-10-04 기준 Python 3.12 실행 절차다.

예제는 하나의 1.2 GHz proton 자기장, 두 nitrogen RF amplitude(25/100 Hz)를
사용한다. RF별로 A1·G2·S3 각각에 105–135 ppm의 offset 147개가 있다.
전체 관측점 882개와 peak별 R1H/R2H를 포함한 파라미터 24개를 fitting한다.
“두 RF”는 같은 자기장에서의 두 amplitude다. 합성 데이터 재현은 프로그램
사용법과 실행을 확인하며, 실험 시료에 대한 검증을 대신하지 않는다.

## 1. 터미널을 열고 설치하기

Git과 Python 3.12가 설치된 macOS 또는 Linux에서 실행한다. Windows에서는
WSL의 Linux 환경을 사용할 수 있다. Python 입력창이 아닌 터미널에 명령을
입력하고, 코드 앞에 `$`를 추가하지 않는다. 이미 저장소가 있으면 clone은
생략하고 `run.py`가 있는 저장소 최상위 폴더로 이동한다.

```bash
git clone https://github.com/Gohyang-Matzip/sbonest.git
cd sbonest
```

`.venv`에 작동하는 Python 3.12 환경이 없다면 아래처럼 만든다.
기존의 다른 환경을 그 자리에서 덮어쓰지는 않는다.
`.venv`의 Python 버전이 다르면 새로 받은 저장소에서 이 가이드를 진행한다.

```bash
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements-sideband.txt -c constraints-sideband.txt
.venv/bin/python --version
```

마지막 줄에 Python 3.12.x가 나와야 한다. `optimalcontrol-nmr`도 함께 설치된다.
이후 단계는 **같은 터미널을 유지하고 저장소 최상위 폴더에서** 진행한다.
모든 명령이 `.venv/bin/python`을 사용하므로 별도의 `activate`는 필요 없다.

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg
```

## 2. 새 연습 폴더 만들기

도우미 명령이 **새** 폴더에 입력 파일 두 개를 복사하고 `fit.json`을 만든다.
원본은 보존한다. 같은 터미널에서 `SBONEST_DEMO_DIR` 값을 유지한다.
`dummy_01`이 이미 있으면 `dummy_02`처럼 사용하지 않은 이름으로 바꾼다.

```bash
export SBONEST_DEMO_DIR="$PWD/session_artifacts/dummy_01"
.venv/bin/python sb_workflow.py init-demo --out "$SBONEST_DEMO_DIR"
```

새 폴더에는 `fit.json`, `data/noisy_25.txt`, `data/noisy_100.txt`가 생긴다.
데이터 경로는 JSON 폴더 기준이다. `Project Name`은 `.json`으로 끝나는 파일명이
아니라 출력 접두사이며, 도우미가 절대 경로로 지정하므로 결과가 한곳에 모인다.
이 도우미는 별도 waveform 파일이 없는 제공 예제를 사용한다.

Fitting 전에 입력과 출력 예정 경로를 검사한다.

```bash
.venv/bin/python run.py "$SBONEST_DEMO_DIR/fit.json" --check
```

JSON 요약에서 `valid: true`, 관측점 882개, 자유 파라미터 24개를 확인한다.
실제 입력 경로, 자기장·RF, 자유·고정 파라미터, bounds와 출력 충돌도 표시한다.
파일을 쓰거나 optimizer를 실행하지 않지만 초기값 격자는 모델을 평가할 수 있다.
검사 통과가 fitting 수렴을 보장하지는 않는다.
검사에 `--identifiability`를 추가하면 초기값에서 Jacobian을 한 번 계산해 약하게
결정되는 파라미터와 강한 상관을 나열한다. 이 예제에서는 이미 세 `R1H` 속도와
`R1H/R2H` 쌍이 지목된다.

## 3. Fitting 실행하기

```bash
.venv/bin/python run.py "$SBONEST_DEMO_DIR/fit.json"
```

명령이 끝나고 터미널 입력 상태로 돌아올 때까지 기다린다. Iteration과 함께
proton 파라미터의 약한 식별성 또는 경계 경고가 나올 수 있다. 이 잡음 포함
예제에서는 수렴에 성공해도 이런 경고가 가능하므로 내용을 읽는다.

성공하면 같은 폴더에 다음 출력과 checkpoint 폴더가 추가된다.

| 파일 | 확인할 내용 |
|---|---|
| `fit_result.json` | 파라미터, 불확실성, 진단, 실행 이력 |
| `fit_result.txt` | 텍스트 보고서. Intensity 표는 반올림된 값 |
| `fit_predictions.csv` | 전체 정밀도의 예측값과 잔차 |
| `fit.pdf` | Fitting profile과 표준화 잔차 패널 |
| `fit_data.pdf` | 데이터만 표시한 profile |
| `fit_checkpoint/` | 실행 식별 정보와 완료된 fitting·분석 기록 |

PDF가 필요 없으면 검사와 첫 fitting 모두에 `--no-pdf`를 붙인다. 수치 파일
세 개와 checkpoint는 생성된다. 새 실행은 2단계의 폴더명을 `dummy_02`처럼
바꿔 시작하고 이전 폴더는 보존한다.
코어가 여러 개인 컴퓨터라면 fitting 명령에 `--workers 4`(또는 사용할 수 있는
코어 수)를 붙인다. 결과는 동일하고 fitting·restart·scan·bootstrap replicate가
더 빨리 끝난다. 1단계의 thread 환경변수는 1로 유지한다.

중단한 실행은 설정을 바꾸지 않고 다음처럼 재개한다.

```bash
.venv/bin/python run.py "$SBONEST_DEMO_DIR/fit.json" --resume
```

첫 실행에 `--no-pdf`를 사용했다면 재개 명령에도 붙인다. 설정, 입력·waveform,
실행 소스, Python·package, 플랫폼, thread 설정, PDF 모드가 같아야 재개된다.
완료된 기본 fitting, restart, profile 점, bootstrap 반복은 재사용하고 미완료
작업을 이어서 수행한다. 관계없는 기존 출력은 보호한다. `fit_checkpoint/`의
기록은 checksum으로 검증하므로 편집하지 않고 온전히 보관한다. 중단된 fitting은
이어서 실행하지만 이미 실패로 끝난 fitting은 실패 결과를 복원하고 재시도 없이
실패 종료한다. 실패한 fitting을 다시 시도하려면 설정을 수정하고 새 `Project Name`을
사용한다. 실행환경을 바꿀 때도 새 출력 접두사가 필요하다.

## 4. 결과가 맞는지 확인하기

```bash
.venv/bin/python - <<'PY'
import json
import os
from pathlib import Path

folder = Path(os.environ["SBONEST_DEMO_DIR"])
r = json.loads((folder / "fit_result.json").read_text())
assert r["success"], r.get("message")
assert (r["n_points"], r["n_parameters"], r["dof"]) == (882, 24, 858)
print("kex (s^-1):", r["kex"], "SE:", r["derived_se"]["kex"])
pb_se = r["derived_se"]["pB"]
print("pB (%):", 100 * r["pB"], "SE (percentage points):",
      None if pb_se is None else 100 * pb_se)
print("RF scale:", r["parameters"]["v1n_scale"]["value"])
print("Actual RF (Hz):", r["v1n_hz"])
print("Reduced chi2:", r["chi2"] / r["dof"])
print("Rank:", r["jacobian_rank"], "of", r["n_parameters"])
print("At bounds:", r["at_bounds"])
print("Warnings:", r["warnings"])
print("Provenance fields:", sorted(r["provenance"]))
PY
```

아래 값과 대략 일치하는지 확인한다. 플랫폼에 따라 미세한 차이는 가능하다.
근거는 [보존된 두 RF fitting](results/auto_H_refit/fits/two_RF_result.json)과
[합성 데이터 생성 파라미터](results/peakwise_H_fit/summary.json)다.

| 항목 | 대표 fitting 값 | 합성 데이터 생성값 |
|---|---:|---:|
| kex (s⁻¹) | 299.68 | 300 |
| pB (%) | 4.983 | 5.0 |
| RF scale | 1.08227 | 1.08 |
| 실제 RF (Hz), 입력 파일 순서 | 27.057 / 108.227 | 27 / 108 |
| Reduced χ² | 0.93569 | — |
| Jacobian rank | 24 / 24 | — |

JSON의 `pB`와 SE는 **분율**이다. 0.05는 5%이며, 두 값에 각각 100을 곱하면
percent와 percentage point 단위의 SE가 된다. `derived_se`에는 두 교환 속도의
covariance가 포함되어 있으므로 두 속도의 오차를 독립적으로 더하지 않는다.
국소 SE는 입력한 intensity의 절대 오차를 사용하고 reduced χ²로 재조정하지
않는다. `null`은 산출 불가이며, 고정 파라미터의 SE=0은 가정을 나타낸다.
Rank가 충분해도 모든 H rate가 정밀하다는 뜻은 아니다. 이 오차는 모델 내
불확실성이며, 잘못된 고정 입력이나 모델 불일치는 포함하지 않는다.

## 5. 그림을 확인하고 저장 결과 보고서 만들기

2단계에서 출력한 연습 폴더를 파일 관리자에서 연다. `fit.pdf`와 `fit_data.pdf`를
더블클릭한다. 각각 peak당 한 페이지로, 이 예제에서는 세 페이지다.
`--no-pdf`를 사용했어도 아래 보고서 명령으로 PDF를 만들 수 있다.
Fit PDF는 위에 전체 profile, 아래에 표준화 잔차를 표시한다. RF 범례는
fitting된 amplitude이며, 데이터 전용 PDF의 RF 범례는 명목 amplitude다.

두 dip, sideband, 먼 offset의 baseline 부근에 잔차가 반복되는 패턴을 보이는지
확인한다. Reduced χ²가 작다는 것만으로 모델의 타당성이 입증되지는 않는다.
계산선은 측정 offset의 예측값을 연결한 것이며 더 촘촘한 simulation이 아니다.

저장한 JSON과 해당 predictions CSV로 보고서를 만든다.

```bash
.venv/bin/python sb_workflow.py report "$SBONEST_DEMO_DIR/fit_result.json" \
  --out "$SBONEST_DEMO_DIR/report_01"
```

명령은 CSV 잔차 계산과 JSON의 관측점 수·χ²·실제 RF가 일치하는지 확인한다.
`report_01_summary.json`, `report_01_summary.txt`, `report_01.pdf`를 만들며
잔기·dataset별 잔차 통계를 포함한다. 재fitting하거나 원래 fitting 입력을
요구하지 않는다. 결과 JSON과 predictions CSV는 함께 보관한다. CSV 이름을
바꿨다면 `--predictions`로 경로를 지정한다. 기존 보고서는 보호하므로 다시
생성할 때는 `report_02`를 사용한다.

`residual_sigma = (관측값 - 예측값) / sigma`다. CSV에는 잔기, dataset index,
자기장, 시간, 명목·실제 RF, offset도 저장된다. 연습 폴더 전체를 보관한다.
입력 사본·설정·JSON 실행 이력도 결과와 함께 있어야 한다.
`session_artifacts/`는 Git에서 제외되므로 필요하면 따로 백업한다.

## 6. 선택: 다른 초기값, kex scan, bootstrap 실행하기

24개 파라미터를 자유롭게 두고 추가 fitting을 수행한다. 기본 초기값과 다른
초기값으로 **두 번** 시도하고, kex 격자 세 점의 제약 fitting과 parametric
bootstrap 5회를 실행한다. 원래 데이터는 유지한다. 아래 블록은 연습 폴더당
한 번 사용하며 기존 `analysis.json`을 보호하고 새 출력 접두사를 쓴다.

```bash
.venv/bin/python - <<'PY'
import json
import os
from pathlib import Path

folder = Path(os.environ["SBONEST_DEMO_DIR"])
config = json.loads((folder / "fit.json").read_text())
config["Project Name"] = str(folder / "analysis")
config["init"]["multistart"] = {
    "starts": [{"kab": 20.0, "kba": 380.0, "v1n_scale": 1.0}]
}
config["init"]["profile"] = {"kex": [295.0, 300.0, 305.0]}
config["init"]["bootstrap"] = {"replicates": 5, "seed": 20261004, "confidence": 0.95}
with (folder / "analysis.json").open("x") as stream:
    json.dump(config, stream, indent=2)
    stream.write("\n")
PY
.venv/bin/python run.py "$SBONEST_DEMO_DIR/analysis.json" --check --no-pdf
.venv/bin/python run.py "$SBONEST_DEMO_DIR/analysis.json" --no-pdf
```

모든 초기값·격자점·반복 결과를 담은 보고서를 만든다.

```bash
.venv/bin/python sb_workflow.py report "$SBONEST_DEMO_DIR/analysis_result.json" \
  --out "$SBONEST_DEMO_DIR/analysis_report_01"
```

`analysis_report_01_summary.txt` 또는 PDF를 연다. `selected=true`는 수렴한
초기값 중 χ²가 가장 작은 해다. 주 fitting 성공이 모든 profile 점이나 bootstrap
반복의 성공을 뜻하지 않는다. 실패한 profile 점은 null χ²를 보존하고, 음수
Δχ²는 기준보다 좋은 해를 찾았다는 뜻이다. 세 점 격자는 신뢰구간이나 전역
최적해의 증명이 아니다.

Bootstrap은 선택된 해와 입력한 절대 sigma로 seed를 지정한 Gaussian 데이터를
만들어 각각 다시 fitting한다. `bootstrap.samples`는 성공·실패를 모두 보존하고
`bootstrap.intervals`는 성공한 반복값의 percentile을 담는다. 5회 반복은 사용법
시연일 뿐이며, 성공한 반복이 100회 미만이면 꼬리 추정 경고가 나온다. 구간은
선택 모델·고정 입력에 조건부이고 모델 불일치를 포함하지 않으며 명목 신뢰수준의
coverage를 보장하지 않는다. 추가 설정과 합성 국소 SE coverage는 매뉴얼
8.1–8.2절을 참고한다.

중단되면 변경하지 않은 `analysis.json`으로 `--resume --no-pdf`를 사용한다.
초기값이 모두 실패하면 명령은 실패 종료하고 `analysis_result.json`과 checkpoint에
시도 기록을 남긴다. `--resume`는 그 실패를 복원하고 시도를 반복하지 않은 채
실패 종료한다. 기록을 확인하고 설정을 수정한 뒤 새 접두사로 실행한다. 실패
결과에는 보고서에 필요한 predictions CSV가 없다.

## 7. 선택: 단일 fitting의 실행 비용 측정하기

```bash
.venv/bin/python benchmark.py "$SBONEST_DEMO_DIR/fit.json" profile \
  --profile-output "$SBONEST_DEMO_DIR/fit.prof"
.venv/bin/python - <<'PY'
import os
import pstats
from pathlib import Path

path = Path(os.environ["SBONEST_DEMO_DIR"]) / "fit.prof"
pstats.Stats(str(path)).sort_stats("cumulative").print_stats(10)
PY
```

벤치마크는 fitting 한 번을 실행해 profiling 정보를 저장하고 기존 fitting
출력은 유지한다. 선택적 multistart·profile likelihood·bootstrap은 실행하지 않는다.
기존 `.prof` 파일은 보호되므로 다시 측정할 때는 새 이름을 사용한다.
Profiler 자체의 비용과 PC 환경이 시간에 영향을 준다. 일반 실행 시간과
동일하게 보거나 다른 환경의 시간과 단순 비교하지 않는다.

## 8. 문제가 생겼을 때

| 증상 | 다음 조치 |
|---|---|
| `python3.12: command not found` | Python 3.12 설치 후 1단계 반복 |
| `.venv/bin/python` 또는 `run.py`를 찾지 못함 | 저장소 최상위 폴더로 이동하고 1단계 확인 |
| `No module named optimalcontrol` | 동일한 `.venv/bin/python`으로 1단계 설치 명령 실행 |
| 터미널을 새로 열었더니 `SBONEST_DEMO_DIR`가 비어 있음 | 저장소 최상위 폴더로 이동해 thread export를 반복하고 `export SBONEST_DEMO_DIR="/2단계에서/출력된/절대경로"` 지정 후 기존 결과 확인. 새 실행은 2단계 사용 |
| `Output already exists` 또는 `FileExistsError` | 기존 실행은 보존하고 2단계로 새 폴더 생성 |
| JSON을 옮긴 뒤 입력 파일을 찾지 못함 | 입력 사본의 `data/` 폴더를 JSON과 함께 두거나 dataset 경로 수정 |
| Checkpoint identity mismatch | 원래 설정·실행환경을 복원해 재개하거나 기존 실행을 보존하고 새 접두사 사용 |
| 최대 평가 횟수 초과 | 초기값·단위·경고 확인. 제한을 늘리기 전에 매뉴얼 10절 참고 |
| H-rate·rank·경계 경고 | 4단계의 오차 한계와 잔차 확인. 수렴만으로 판단하지 않음 |

## 9. 실험 데이터로 넘어가기

별도 설정·출력 폴더를 사용하고 매뉴얼 2–6절과 9절을 따른다. 합성 intensity,
양수인 절대 오차, 잔기명, ¹H shift, 두 핵의 자기장 주파수, saturation time,
pulse 조건을 측정값으로 바꾼다. 데이터 첫 줄은 **¹⁵N MHz**, offset은 절대
¹⁵N ppm, 시간은 초, RF는 Hz다. 제외할 잔기는 명시적으로 `off`로 지정한다.
R1H/R2H를 모두 생략하면 fitting하고, 둘 다 입력하면 고정한다. 자동 H-rate
fitting은 현재 하나의 proton 자기장을 지원한다. 입력 자료를 보존하고 모델
적합성을 검토한 뒤 실험의 교환 파라미터를 해석한다.

실행 절차의 구현: [sb_workflow.py](sb_workflow.py), [sbfit.py](sbfit.py),
[sb_checkpoint.py](sb_checkpoint.py), [sb_report.py](sb_report.py),
[sb_bootstrap.py](sb_bootstrap.py).
