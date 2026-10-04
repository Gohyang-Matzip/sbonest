# SBONEST fitting 매뉴얼

[English](SBONEST_MANUAL.md) · [README](README.md) · [기술 설명](SIDEBAND.md)

이 매뉴얼은 현재 코드의 CLI `Sideband` 모델을 설명하며, 2026-10-03에
실행을 확인했다. OC의 NH spin operator를 사용해 ¹H decoupling sideband를
포함한 두 상태 ¹⁵N CEST profile을 fitting한다. 주 예제는
90°x–240°y–90°x 반복 decoupling, 1.2 GHz 장비, 전체 30 ppm offset 범위다.
제공된 예제 데이터는 모두 **합성 데이터이며 실제 측정 데이터가 아니다.**

복사되어 있는 [ONEST 매뉴얼](MANUAL.ko.md)은 다른 모델용이다.
기존 웹 인터페이스, `prepare.py`, `mcrun.py`는 현재 `Sideband` fitting의
지원 실행 경로가 아니다.

**R1H/R2H는 입력하지 않아도 된다.** Decoupling 설정에서 두 값을 생략하면
프로그램이 peak별로 독립 fitting한다. H-rate 입력이 없는
[이식 가능한 두 RF 설정](example/sideband_auto_H/two_RF.json)과
[전체 profile 재fitting 결과](results/auto_H_refit/fits/REFIT_REPORT.md)를 참고한다.

## 1. 실행환경과 첫 fitting

Python 3.12를 CI 기준으로 사용한다. 새로 받은 저장소에서는 다음처럼 설치한다.

```bash
git clone https://github.com/Gohyang-Matzip/sbonest.git
cd sbonest
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-sideband.txt
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg
.venv/bin/python run.py example/sideband_auto_H/two_RF.json
```

OC는 `optimalcontrol-nmr` 패키지로 설치되므로 인접한 OC 폴더가 필요 없다.
이미 작동하는 환경은 그대로 사용한다. 수치 결과만 필요하면 마지막 명령에
`--no-pdf`를 추가한다.

이 예제는 25/100 Hz nitrogen RF에서 3개 peak를 fitting한다. 각 peak·RF마다
105–135 ppm에 균일한 offset 147개가 있으며 실제 간격은 약 24.985 Hz다.
Sideband를 포함한 전체 관측점 882개에 대해 공통 교환 속도 2개, RF scale
1개, peak별 파라미터 7개(질소 5개와 R1H/R2H)를 합한 24개를 fitting한다.
**H-rate 입력은 필요 없다.** 두 RF란 같은 1.2 GHz proton 자기장에서의
두 saturation amplitude를 뜻하며, 서로 다른 두 자기장을 뜻하지 않는다.

| 항목 | 보존된 fitting 결과 | 합성 데이터 생성값 |
|---|---:|---:|
| `kex` (s⁻¹) | 299.68025 | 300 |
| `pB` (분율) | 0.0498302 | 0.05 |
| `v1n_scale` | 1.0822687 | 1.08 |
| Reduced χ² | 0.93569 | — |
| Jacobian rank | 24 / 24 | — |

실행환경에 따라 미세한 수치 차이는 가능하다. 이 잡음 포함 예제에서는
H-rate의 약한 식별성과 경계 경고가 나올 수 있다. Jacobian rank가 충분해도
개별 H rates가 정밀하게 결정된다는 뜻은 아니다. 예제 재현은 프로그램의
실행을 확인하는 것이며 실험 시료에 대한 검증을 대신하지 않는다.

출력 접두사는 `results/auto_H_two_RF`이며 파일 종류는 8절을 참고한다.
재실행하려면 JSON 사본의 `Project Name`을 새 접두사로 변경한다. 기존
출력은 덮어쓰지 않는다. 25/50/100 Hz는
`example/sideband_auto_H/three_RF.json`을 사용한다(1323개 관측점, 24개 파라미터).

이전의 single-peak·H-rate 고정 예제는
[full.json](manuscript/sideband_30ppm/results/1200/full.json)에 보존되어 있으며
관측점 441개, fitting parameter 8개다. `results/`의 보존된 JSON에는 원래
PC의 절대 경로가 있을 수 있으므로 새 실행에는 위의 이식 가능한 설정을 사용한다.

## 2. 입력 데이터와 단위

자기장, saturation time, 명목 ¹⁵N RF의 조합마다 텍스트 파일 하나를 사용한다.
한 파일에는 여러 잔기를 넣을 수 있다.
[full_25.txt](manuscript/sideband_30ppm/results/1200/full_25.txt)의 앞부분은 다음과 같다.

```text
121.5949416
0.4
25 0
# offset intensity error
# A1 R2a: 13 R2b: 18 dw: 2.7
105 0.54861876445228 0.001
105.20547945205 0.5469790572768 0.001
105.41095890411 0.54623604725812 0.001
```

위 내용은 형식 설명을 위한 발췌이며, fitting용 전체 데이터가 아니다.

| 위치 | 의미와 단위 |
|---|---|
| 첫째 줄 | 양수인 **¹⁵N Larmor 주파수, MHz**. ¹H 장비 주파수를 넣지 않는다 |
| 둘째 줄 | Saturation duration `T`, 초 |
| 셋째 줄 첫 값 | 명목 ¹⁵N RF amplitude `v1`, Hz |
| 셋째 줄 둘째 값 | RF 불균일성 폭 `v1err`, Hz. 0이면 RF 평균화를 하지 않는다 |
| 넷째 줄 | 반드시 필요한 컬럼 제목 줄. Parser가 이 줄을 건너뛴다 |
| 잔기 헤더 | 잔기명과 초기 `R2a`, `R2b` (s⁻¹), `dw` (ppm) |
| 데이터 첫 열 | **절대 ¹⁵N chemical-shift offset, ppm** |
| 데이터 둘째 열 | 정규화한 intensity `I/I0` |
| 데이터 셋째 열 | 같은 intensity 단위의 양수인 절대 표준오차 |

따라서 `0.001`은 `I/I0`의 오차가 0.001이라는 뜻이며, 0.001%라는 뜻이 아니다.
Field, `T`, `v1`, intensity error는 양수여야 하고 `v1err`는 0 이상이어야 한다.
모든 수치 입력은 유한해야 한다.

다른 잔기는 `# G2 R2a: 12 R2b: 15 dw: -2` 같은 헤더와 데이터 행을 이어서
추가한다. 같은 잔기를 한 파일에 중복해서 넣지 않는다. `A1`, `G23` 등의
잔기명은 JSON과 정확히 일치해야 한다. 한 파일에서는 헤더 형식을 통일한다.
`# A1` 같은 짧은 헤더도 가능하지만 초기값이 `R2a=10`, `R2b=100`, `dw=0.1`이
되므로, full header에 적절한 값을 쓰는 편이 좋다.
잔기 데이터 중간에 별도의 주석 줄을 넣지 않는다. `#`로 시작하는 줄은
해당 데이터 블록의 끝으로 처리된다. 그 뒤의 데이터 행에는 새 잔기 헤더가
필요하며, 헤더가 없으면 점을 조용히 누락하는 대신 오류를 낸다. 빈 데이터
파일, 관측값이 없는 잔기 헤더, 필수 넷째 줄의 컬럼 제목 대신 잔기 헤더나
데이터 행을 넣은 파일도 거부한다. 완결된 잔기 블록 사이의 주석과 빈 줄은
계속 사용할 수 있다.

모델의 초기 상태는 `(pA Nz, pB Nz)`이고 검출값은 `Nz(A,T)/pA`다.
평형 자화로 회복시키는 source term은 없다. 이 모델에서 saturation 영향이
없는 먼 offset의 신호는 `exp(-R1*T)`로 접근한다. 실험 reference와 정규화가
이 convention에 맞는지 확인하고, 먼 offset의 baseline을 일괄적으로 1로
재정규화하지 않는다.

## 3. 전체 설정 파일

아래 내용을 `full_25.txt`, `full_50.txt`, `full_100.txt`가 있는 폴더에
`config.json`으로 저장하거나, 데이터 경로를 사용할 파일로 바꾼다.
아래 수치는 합성 예제의 조건이므로 실제 실험에서는 측정 조건으로 바꿔야 한다.

```json
{
  "Project Name": "results/my_sideband_fit",
  "datasets": ["full_25.txt", "full_50.txt", "full_100.txt"],
  "residues": [{"name": "A1", "flag": "on"}],
  "init": {
    "Method": "Sideband",
    "kex": {"min": 200, "max": 400, "nsteps": 3},
    "pB": {"min": 0.03, "max": 0.07, "nsteps": 3},
    "initial": {"A1.peak_ppm": 120.05, "A1.R1": 1.2},
    "max_nfev": 500
  },
  "sideband": {
    "decoupling": {
      "h_larmor_mhz": 1200.0,
      "h_carrier_ppm": 8.5,
      "p90_s": 0.000070,
      "cycle": "RR",
      "b1_scale": 1.0,
      "J_hz": 92.0
    },
    "residues": {
      "A1": {"h_ppm_a": 6.2, "h_ppm_b": 6.5}
    },
    "v1n": {"mode": "scale", "initial": 1.0, "bounds": [0.8, 1.2]}
  }
}
```

프로젝트 폴더에서 실행한다.

```bash
.venv/bin/python run.py /absolute/path/to/config.json
```

`sbfit.py`도 같은 방식으로 실행할 수 있다. 수치 결과만 필요하면 명령 뒤에
`--no-pdf`를 붙인다. JSON에는 주석이나 마지막 항목 뒤의 쉼표를 넣지 않는다.

**경로 규칙:** 입력 데이터와 waveform의 상대 경로는 설정 파일이 있는 폴더를
기준으로 한다. 상대 경로로 지정한 `Project Name`은 **명령을 실행한 작업
디렉터리**를 기준으로 한다. 이 값은 확장자가 붙은 파일명이 아니라 출력
접두사다. 같은 이름의 결과가 있으면 중단하므로 새 접두사 또는 새 작업
폴더를 사용해 이전 결과를 보존한다.

**잔기 선택:** 읽어들인 잔기는 기본적으로 모두 활성화되어 있다.
최상위 `residues` 목록에서 생략한다고 비활성화되지 않는다. 제외할 잔기는
명시적으로 `"flag": "off"`로 설정하고, 모든 활성 잔기에 ¹H shift를 제공한다.

## 4. 실제 decoupling 순서 지정

이 매뉴얼에서는 `R = 90°x–240°y–90°x`이며 `cycle="RR"`은 같은 R 블록을
반복한다. `RRbar`는 위상 순서가 다른 조건이며 RR과 같은 뜻이 아니다.

| 설정 | 의미 |
|---|---|
| `h_larmor_mhz` | ¹H Larmor 주파수, MHz. 명목 1.2 GHz 예제에서는 1200 |
| `h_carrier_ppm` | ¹H decoupler carrier, ppm |
| `p90_s` | 명목 90° pulse 길이, 초. 70 μs는 `0.000070` |
| `cycle` | 예제는 `RR`. `R`, `RRbar`, `MLEV4`도 지원 |
| `b1_scale` | ¹H RF amplitude 배율, 기본값 1. `v1n_scale`과 다름 |
| `J_hz` | NH scalar coupling, Hz. 기본값 92 |
| `R1H`, `R2H` | 선택적인 기존 고정값, s⁻¹. 둘 다 생략하면 peak별 자동 fitting |
| `h_ppm_a`, `h_ppm_b` | 활성 잔기별 상태 A/B의 ¹H chemical shift, ppm |

세 pulse의 길이는 `p90_s`, `(8/3)*p90_s`, `p90_s`다.
Proton RF amplitude는 `b1_scale/(4*p90_s)` Hz이며, `b1_scale`을 바꿔도
pulse 길이는 유지된다. 70 μs, scale 1이면 amplitude는 약 3571.43 Hz이고
R 블록 길이는 326.67 μs다.

Sideband profile은 두 상태 NH 시스템을 실제 pulse segment에 따라 전파해
계산한다. RF, offset, J coupling, relaxation, exchange가 동시에 작용하며,
`T`의 마지막에 남는 pulse 일부도 포함한다. 별도의 sideband amplitude를
fitting하거나 sideband 위치를 수동으로 지정하지 않는다.

Pulse 조건과 H chemical shift는 고정 입력값이다. `decoupling`과 dataset
override에서 R1H/R2H를 모두 생략하면 `run.py`가 두 값을 peak별 nuisance
parameter로 함께 fitting한다. 초기값은 (2,25) s⁻¹이며, 측정값이나 임의의
고정값이 아니다. 비음수로 제한하며 같은 peak의 RF 파일과 A/B 상태는
같은 값을 공유한다. 서로 다른 peak끼리는 공유하지 않는다.
현재 자동 모드는 단일 proton 자기장만 지원한다.

추가 설정은 필요 없다. 고급 사용자는 `sideband` 안에
`"proton_relaxation": {"mode": "fit"}` 또는 `{"mode": "fixed"}`를 지정할 수 있다.
H rate를 하나라도 명시한 기존 설정은 mode를 따로 주지 않으면 고정 동작을
유지한다. `mode="fit"`과 decoupling의 명시적 H rates를 함께 주면 오류를 낸다.
자동 모드의 시작값·범위는 필요할 때만 `init.initial`/`init.bounds`에서
`A1.R1H`, `A1.R2H`처럼 지정한다. 일반 자동 workflow에서는 `init.vary`를
생략한다. 명시적인 `vary` 목록은 여전히 해당 파라미터만 fitting한다.

결과에는 H 추정값·표준오차·상관관계와 경계/약한 식별성 진단을 남긴다.
교환/RF 결과가 안정적이어도 개별 H rates는 부정확할 수 있다. H rates와의
공분산은 다른 파라미터의 오차 계산에 포함된다. 이전 식별성 검사는
[proton relaxation 보고서](PROTON_RELAXATION_FIT_CHECK.md)에 있다.
Minor state의 ¹H shift를 모르면 가정한 값을 명시하고 가능한 범위에서
바꿔가며 fitting을 반복한다. ¹H shift나 decoupling amplitude가 틀리면
그 오차를 nitrogen RF나 교환 파라미터가 대신 흡수할 수 있다.

파일별 조건이 다르면 `sideband.datasets`에 decoupling override 목록을
추가한다. 파일이 3개일 때의 예는 다음과 같다.

```json
[{"h_carrier_ppm": 8.4}, {}, {"h_carrier_ppm": 8.6}]
```

목록 길이는 최상위 `datasets`의 길이와 같아야 한다. 각 항목은 공통
`sideband.decoupling` 값을 덮어쓰며, 반드시 데이터 파일 순서와 일치해야 한다.

## 5. Nitrogen RF, ν1N의 fitting 또는 고정

`sideband.v1n` 객체 전체를 다음 중 하나로 교체한다.

| Mode | RF fitting parameter | 의미 |
|---|---|---|
| `fixed` | 없음 | 각 파일의 명목 `v1` 사용 |
| `scale` | `v1n_scale` | 모든 파일과 활성 잔기에 공통 배율 하나 적용 |
| `per_dataset` | `v1n[0]`, `v1n[1]`, … | 파일마다 RF amplitude를 Hz로 fitting; 파일 내 잔기는 공유 |

RF 고정:

```json
{"mode": "fixed"}
```

공통 scale fitting:

```json
{"mode": "scale", "initial": 1.0, "bounds": [0.8, 1.2]}
```

25/50/100 Hz 파일 3개의 RF를 각각 fitting:

```json
{
  "mode": "per_dataset",
  "initial": [25, 50, 100],
  "bounds": [[20, 30], [40, 60], [80, 120]]
}
```

`scale`에서 실제 RF는 `명목 v1 × v1n_scale`이다. `per_dataset`의 index는
0부터 시작하며 입력 파일 순서를 따른다. 초기값과 bounds는 양수여야 하고
초기값은 bounds 안에 있어야 한다. 범위를 생략하면 scale은 0.5–1.5,
파일별 RF는 명목값의 0.5–1.5배를 사용한다.

`v1err`는 Gaussian RF 불균일성 sampling 폭이다. RF fitting 중에는
`v1err/v1` 비율을 유지한다. 이 값은 **fitting된 평균 RF의 불확도가 아니다.**
Fitting parameter의 표준오차는 결과의 `stderr`에 따로 기록된다.

합성 예제로 ν1N을 fitting할 수 있음을 확인했다. 실측에서 식별 가능한지는
sampling, SNR, RF level 수, 고정된 proton 입력의 정확도에 달려 있다.
먼저 RF 고정과 공통 scale fitting을 비교한다. 파일마다 다른 보정이 필요하다는
calibration 근거가 있고 진단 결과도 이를 뒷받침할 때 파일별 RF를 사용한다.

## 6. 파라미터 이름, 초기값, 범위, 고정

| 파라미터 이름 | 공유 범위 | 단위와 정의 |
|---|---|---|
| `kab` | 모든 활성 잔기와 파일 | A → B 교환 속도, s⁻¹ |
| `kba` | 모든 활성 잔기와 파일 | B → A 교환 속도, s⁻¹ |
| `v1n_scale` 또는 `v1n[i]` | 위에서 선택한 방식 | 무차원 또는 Hz |
| `A1.peak_ppm` | A1 잔기의 모든 파일 | 상태 A의 ¹⁵N chemical shift, ppm |
| `A1.dw_ppm` | A1 잔기의 모든 파일 | δB − δA, ppm. 양수와 음수 모두 가능 |
| `A1.R1` | A1 잔기의 두 상태 | ¹⁵N longitudinal rate, s⁻¹ |
| `A1.R2a`, `A1.R2b` | A1 잔기의 상태 A/B | ¹⁵N transverse rate, s⁻¹ |
| `A1.R1H`, `A1.R2H` | 자동 모드; A1의 두 상태와 RF 파일 | 내부 proton relaxation rate, s⁻¹ |

결과에는 `kex = kab + kba`, `pB = kab/kex`를 보고한다.
`pB=0.05`는 5%다. `init.initial`, `init.bounds`, `init.vary` 안의
파라미터 이름으로는 `kex`, `pB`를 사용할 수 없고 `kab`, `kba`를 사용한다.

`init.kex`와 `init.pB`는 시작 교환 속도를 고르는 grid를 정의한다.
여기의 min/max는 **최종 fitting 범위를 제한하지 않는다.** 잔기 초기값은
첫 spectrum과 헤더에서 얻으며, `init.initial`은 grid 단계 이후의 초기 벡터를
덮어쓴다. Sideband를 main minimum으로 오인할 수 있다면 적절한
`peak_ppm` 초기값을 명시하는 것이 도움이 된다.

범위를 좁히려면 다음과 같은 객체를 `init.bounds`로 추가한다.

```json
{"kab": [0, 100], "kba": [1, 1000], "A1.dw_ppm": [1, 5]}
```

기본 물리적 범위는 `kab ≥ 0`, `kba ≥ 1e-8`, relaxation rate ≥ 0이며,
shift에는 기본 제한이 없다. 사용자가 지정하는 bounds는 이 범위와
`sideband.v1n`의 RF bounds 안에 있어야 한다. 하한은 상한보다 작아야 하며
모든 초기값은 범위 안에 있어야 한다. 파라미터를 고정할 때는 하한과 상한을
같게 하지 말고 `vary`를 사용한다.

**`init.vary`를 생략하면 모든 파라미터를 fitting한다. `initial`에 값만
쓴다고 고정되지 않는다.** 비어 있지 않은 `vary` 목록을 주면 그 이름만
fitting하고 나머지는 초기값으로 고정한다. 합성 예제에서 RF만 복원하는
시험을 하려면 `init` 객체 전체를 다음으로 교체한다.

```json
{
  "Method": "Sideband",
  "kex": {"min": 300, "max": 300, "nsteps": 1},
  "pB": {"min": 0.05, "max": 0.05, "nsteps": 1},
  "initial": {
    "kab": 15,
    "kba": 285,
    "v1n_scale": 1.0,
    "A1.peak_ppm": 120.0,
    "A1.dw_ppm": 3.0,
    "A1.R1": 1.5,
    "A1.R2a": 12.0,
    "A1.R2b": 15.0
  },
  "vary": ["v1n_scale"],
  "max_nfev": 500
}
```

이는 교환과 relaxation을 이미 안다고 가정한 calibration 예제다.
실험 데이터에서 이 값들을 추정하는 절차를 대신하지 않는다. 관측점 수는
자유 파라미터 수보다 많아야 한다. 여러 물리적으로 타당한 초기값으로
별도 fitting을 수행한다. 국소 최적화이므로 전역 최소점을 보장하지 않는다.

## 7. 전체 30 ppm에 대한 균일 offset 간격

제공된 1.2 GHz 예제는 120 ppm을 중심으로 105–135 ppm을 사용한다.
¹⁵N 주파수 `νN`을 MHz로 표시하면 30 ppm의 폭은 `30 × νN` Hz다.
`νN = 121.5949416 MHz`에서는 약 3647.85 Hz다.

이 연구 예제는 양 끝과 중심점을 유지하면서 25 Hz를 목표 간격으로 사용한다.
중심이 `δ0`이면 다음처럼 계산한다.

```text
Nintervals = 2 × ceil(15 × νN / 25)
δi = δ0 − 15 + 30 × i / Nintervals,  i = 0, …, Nintervals
실제 간격 (Hz) = 30 × νN / Nintervals
```

| 명목 ¹H 주파수 | 예제의 ¹⁵N 주파수 (MHz) | RF profile당 점 수 | 실제 간격 (Hz) |
|---|---:|---:|---:|
| 600 MHz | 60.7974708 | 75 | 24.648 |
| 800 MHz | 81.0632944 | 99 | 24.815 |
| 1200 MHz | 121.5949416 | 147 | 24.985 |

실험에서는 장비의 실제 ¹⁵N 주파수로 변환한다. 25 Hz는 이 연구의 기준
측정 간격이며, 모든 좁은 dip을 충분히 분해한다는 보편적인 기준은 아니다.
필요한 구간은 더 촘촘한 pilot profile로 확인한다. 30 ppm과 대략 일정한 Hz
간격을 유지하면 고자장에서 측정점 수가 늘어나므로, 점당 scan 수가 같아도
총 측정 시간은 같지 않다.

Fitter는 입력 offset을 그대로 사용한다. 보간, 재sampling, 균일 간격 강제,
sideband 자동 masking은 하지 않는다. 전체 profile 분석에서는 유효한 모든
측정점을 입력한다. Masking과 비교하려면 선택한 행을 제거한 별도 입력 파일을
만들고, 정규화와 오차 추정은 유지한 채 별도 config와 출력 접두사를 사용한다.
파라미터 불확도와 residual 구조를 비교한다. 점 수가 다른 데이터의 raw χ²를
직접 비교해서는 안 된다.

`demo_sideband.py --offset-step 25`는 다른 예제로, ±2400 Hz를 사용하므로
이 자기장에서는 약 39.48 ppm이다. **30 ppm** 시험에는 1절의
`example/sideband_auto_H/two_RF.json` 또는 `three_RF.json`을 사용한다.

### 7.1 600 및 800 MHz에서 SBONEST와 ONEST 비교

[600/800 MHz 비교](results/field_comparison_600_800_20261003_02/REPORT.txt)는
두 자기장을 **각각 따로** fitting한다. 같은 자기장 내에서는 두 모델에 동일한
잡음 포함 합성 관측값, 절대 σ = 0.001, 파라미터 범위와 optimizer를 적용한다.
자유 파라미터는 `kab`, `kba`, `v1n_scale`과 peak의 `peak_ppm`, `dw_ppm`,
`R1`, `R2a`, `R2b`를 합한 8개다.

단일 peak는 120 ppm, ΔδN = 3 ppm, T = 0.4 s이고 명목 nitrogen RF는
25/50/100 Hz다. 105–135 ppm 창에서 RF당 점 수는 600 MHz에서 75개(총 225개),
800 MHz에서 99개(총 297개)다. SBONEST의 R1H = 2, R2H = 25 s⁻¹은
논문의 생성값으로 고정한다. 이 절은 **고정-H 비교**이며, 현재 기본 기능인
자동 peak별 H-rate fitting은 1절을 참고한다.

1절의 의존성 설치와 thread 환경변수 설정 후 저장소 루트에서 실행한다.

```bash
.venv/bin/python compare_field_models.py \
  --source results/field_comparison_600_800_20261003_02/inputs \
  --out results/field_comparison_repeat_01
```

`--out`은 항상 새 폴더로 지정한다. Script의 기본 입력 경로는 과거 로컬 논문
데이터를 가리키므로 새 checkout에서는 위 `--source`를 명시해야 한다.
함께 배포하는 입력 사본은 원본 데이터를 보존하며 다른 환경에서도 사용할 수 있다.

| ¹H 자기장 | 모델 | kex ± SE (s⁻¹) | pB ± SE (%) | Reduced χ² |
|---|---|---:|---:|---:|
| 600 MHz | SBONEST | 297.89 ± 1.63 | 5.0287 ± 0.0299 | 1.058 |
| 600 MHz | ONEST Matrix | 297.95 ± 1.64 | 5.0142 ± 0.0357 | 1.055 |
| 800 MHz | SBONEST | 299.85 ± 1.23 | 4.9994 ± 0.0313 | 1.035 |
| 800 MHz | ONEST Matrix | 301.41 ± 1.23 | 4.9702 ± 0.0330 | 1.058 |

참값은 kex = 300 s⁻¹, pB = 5%, RF scale = 1.08이다. ±는 절대 intensity
오차와 파라미터 간 공분산을 사용한 국소 표준오차 1 SE이며, reduced χ²로
재조정하지 않는다. 표의 pB 오차는 퍼센트포인트 단위이고,
`parameters.csv`의 `pB` 값과 SE는 분율 단위다.

ONEST는 음수 예측을 0으로 제한하는 실제 Matrix worker를 사용하며,
RF scale 처리와 최적화 조건을 SBONEST에 맞춘다. 기존 ONEST CLI의 기본 설정을
그대로 실행한 비교는 아니다. 별도 `N_signed_control`은 이 zero clamp를 제거한다.
모델·잡음 조건별 초기값 3개를 보존하며 총 36회 fitting한다. 최저 χ²를 채택하되,
잡음 포함 ONEST에서는 초기값에 따라 조금 다른 최솟값에 도달했다.
모든 초기값의 수렴점 일치나 전역 최적성을 가정해서는 안 된다.

600 MHz의 잡음 포함 결과는 거의 같다. 800 MHz 무잡음 ONEST의 kex는
301.535 s⁻¹(+0.51%), signed 대조는 301.551 s⁻¹(+0.52%)이고,
SBONEST는 두 자기장에서 모두 300 s⁻¹을 복원한다. 120 ppm 기준
|offset| = 1250–1850 Hz 마스크는 이 창의 관측점을 제외하지 않으므로,
동일한 masked fitting을 반복하지 않았다. 한 파라미터 조합과 자기장별 잡음 표본
하나에 대한 결과다. ONEST의 국소 SE는 모델 불일치를 포함하지 않으며,
다른 peak 위치·측정 창·pulse 조건과 실제 측정 데이터는 별도 검증이 필요하다.

`summary.json`에는 모든 초기값·공분산·진단·소스 해시가, `parameters.csv`에는
전체 파라미터가, `predictions_*.npz`에는 그림의 수치 데이터가 있다.
`verification.json`은 저장 결과를 확인한 기록이다.
[600 MHz](results/field_comparison_600_800_20261003_02/comparison_600.png)와
[800 MHz](results/field_comparison_600_800_20261003_02/comparison_800.png) 그림에서
profile과 표준화 residual을 확인할 수 있다. 새 실행은 fitting·표·그림·입력 사본을
생성한다. `REPORT.txt`와 `verification.json`은 배포 결과에 별도로 작성한 기록이며,
비교 script가 자동 생성하는 파일은 아니다.

## 8. 출력 파일과 결과 해석

`Project Name = results/my_sideband_fit`이면 다음 파일이 만들어진다.

| 파일 | 내용 |
|---|---|
| `results/my_sideband_fit_result.json` | 파라미터, 실제 RF (Hz), 진단값, 입력 config |
| `results/my_sideband_fit_result.txt` | 사람이 읽는 결과와 관측/계산 intensity 표 |
| `results/my_sideband_fit.pdf` | 활성 잔기당 한 페이지의 데이터와 계산 profile |
| `results/my_sideband_fit_data.pdf` | 데이터만 그린 profile |

PDF는 입력된 전체 offset 범위를 보여준다. 계산선은 측정 offset에서 계산한
값을 연결한 것이며, 별도의 촘촘한 simulation이 아니다. 측정점 사이의 좁은
구조는 드러나지 않을 수 있다. 범례의 RF는 **명목값**이며 실제 fitting RF는
JSON의 `v1n_hz`에서 확인한다. 텍스트 intensity 표는 소수 셋째 자리까지
출력하므로 정량 재분석에는 원래 입력 데이터를 보존해 사용한다.
`plot_full_profile.py`는 별도 demo의 파일과 `validation.json`을 전제로 하므로,
임의 fitting 결과에 사용하는 범용 plotter가 아니다.

1절의 예제 실행 후, 프로젝트 폴더에서 다음으로 결과를 확인할 수 있다.

```bash
.venv/bin/python - <<'PY'
import json
from pathlib import Path
r = json.loads(Path("results/auto_H_two_RF_result.json").read_text())
print("kex:", r["kex"], "pB:", r["pB"])
print("actual RF (Hz):", r["v1n_hz"])
print("reduced chi2:", r["chi2"] / r["dof"])
print("rank:", r["jacobian_rank"], "of", r["n_parameters"])
print("warnings:", r["warnings"])
PY
```

교환 파라미터를 해석하기 전에 다음을 확인한다.

1. Sideband까지 전체 profile을 살핀다. Main dip, minor dip, 먼 offset
   baseline에서 체계적인 불일치가 있는지 확인한다.
2. `n_points`, 활성 잔기, RF 파일 순서, `v1n_hz`를 확인한다.
3. `warnings`, `at_bounds`, `jacobian_rank`, `scaled_condition`을 읽는다.
   Rank가 `n_parameters`보다 작으면 각 파라미터를 국소적으로 분리해
   식별할 수 없다. RF correlation의 절댓값이 0.95를 넘거나 scaled
   condition number가 10⁶을 넘으면 경고가 나온다.
4. `v1n_correlations`를 살피고 타당한 다른 초기값으로 반복한다.
   Solver가 수렴했다는 메시지만으로 유일한 물리적 해를 입증하지 못한다.

`stderr`는 입력 intensity error를 절대 σ로 사용하는 국소 선형 covariance
추정값이며, **reduced χ²로 재조정하지 않는다.** Rank가 부족하면 자유
파라미터의 오차는 `null`이다. 고정 파라미터의 `vary=false`, `stderr=0`은
가정을 나타낼 뿐 실험적으로 정밀하게 측정했다는 뜻이 아니다.
이 오차에는 고정 proton 입력의 불확도나 모델 불일치가 포함되지 않는다.
JSON은 전체 covariance matrix를 출력하지 않으므로, `kab`와 `kba`를
독립이라고 가정해 `kex`나 `pB` 오차를 계산하지 않는다.
현재 버전에는 profile likelihood나 bootstrap을 수행하는 내장 CLI 명령이 없다.

## 9. 실험 fitting 순서와 모델 한계

1. 두 핵의 장비 주파수, carrier, pulse 위상과 길이, 보정한 ¹H amplitude,
   saturation time, 명목 ¹⁵N RF, reference 정규화, 잔기 shift, intensity
   error를 기록한다.
2. 분리된 잔기와 여러 RF level로 시작한다. 전체 profile을 확인한 뒤
   RF 고정 또는 공통 scale fitting을 선택한다.
3. 출력 접두사를 구분해 RF 고정과 fitting을 비교한다. Residual,
   파라미터 변화, rank, correlation, calibration과의 일치도를 확인한다.
4. 고정한 ¹H shift, amplitude, coupling을 타당한 범위에서 바꿔 fitting한다.
   자동 모드에서는 peak별 H-rate 초기값을 바꿔 재시작하고 내부 H 진단을
   확인한다. 고정-H 모드를 선택한 경우에는 H rates도 바꿔 본다. 이러한
   민감도는 국소 `stderr`와 별도로 보고한다.
5. 동일 `kab`, `kba`를 공유한다는 과학적 근거가 있을 때 잔기를 함께
   fitting한다. 입력 데이터, config, 비교한 결과를 모두 보존한다.

여러 데이터셋을 함께 fitting하면 각 잔기의 nitrogen shift뿐 아니라
**`R1`, `R2a`, `R2b`도 모든 파일에서 공유**한다. 자기장별 독립 nitrogen
relaxation parameter는 없다. 여러 field를 합치면 이 제약이 적용되므로,
영향을 확인하기 위해 field별로 따로 fitting할 필요가 있다.
`sideband.datasets` override로 이 공유 규칙을 바꿀 수는 없다.

모델은 교환하는 두 상태 각각에 N 하나와 H 하나를 포함하고, 현상론적
product-operator relaxation을 사용한다. 상태 교환에 따라 N과 H의 shift가
모두 달라질 수 있다. 추가 proton, 별도의 proton/water exchange 과정,
CSA–DD cross-correlation, 시간에 따른 RF drift는 포함하지 않는다.
각 segment를 정확히 전파한다고 해서 모든 실험 효과를 정확히 나타내는
모델이 되는 것은 아니다.

현재 1.2 GHz 결과는 sideband 정보를 활용할 수 있다는 합성 데이터의
개념 검증이다. 실제 fitting이 1.2 GHz에서 항상 더 좋아진다는 결론은 아니다.
계획한 실험에서 sampling, SNR, 총 측정 시간, decoupling 모델의 적합성을
함께 검증해야 한다.

## 10. 문제 해결

| 증상 | 확인 및 조치 |
|---|---|
| `No module named optimalcontrol` | `.venv/bin/python` 사용 여부와 `requirements-sideband.txt` 또는 로컬 OC 설치 확인 |
| `Output already exists` | 새 `Project Name`을 지정하거나 새 폴더에서 실행 |
| 파일을 찾지 못함 | 데이터 경로는 JSON 폴더, 출력 경로는 실행 작업 폴더 기준인지 확인 |
| 비교 실행에서 `results/600/full.json` 또는 `results/800/full.json` 누락 | 7.1절의 배포 입력 `--source` 경로를 명시 |
| 잔기 또는 ¹H shift 누락 | 잔기명 정확히 일치시키고 제외할 잔기는 명시적으로 off |
| `Missing column-header line` | 넷째 줄에 컬럼 제목을 유지하고 그 다음에 첫 잔기 헤더를 배치 |
| `outside a residue block` | 잔기 헤더를 확인하고 데이터 행 중간의 별도 주석 줄을 제거 |
| `No residue data found` 또는 `No data points for residue` | 모든 잔기 블록에 관측값을 넣고 파일마다 최소 한 블록을 포함 |
| `v1n` initial/bounds 오류 | Mode와 배열 길이를 맞추고 `fixed`로 바꿀 때 scale용 설정 제거 |
| 파라미터 이름 오류 | 6절 이름 사용. `vary`나 bounds에는 `kex`/`pB` 대신 `kab`/`kba` 사용 |
| 초기값이 bounds 밖 | 헤더의 `dw`/R2와 명시적 초기값이 모든 제한 안에 있는지 확인 |
| 최대 평가 횟수 초과 | 단위, 초기값, 식별성을 먼저 확인한 뒤 필요할 때 `init.max_nfev` 조정 |
| Rank 경고 또는 `stderr: null` | Correlation과 sampling을 살펴 불필요한 자유 파라미터를 줄이거나 정보가 있는 데이터 추가 |
| 실행이 오래 걸림 | 수치 라이브러리 thread를 1로 유지하고 `--no-pdf` 사용. 0이 아닌 `v1err`는 RF 평균화 비용 추가 |

프로젝트 폴더에서 수치 회귀 검증을 다시 실행할 수 있다.

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
  MPLBACKEND=Agg .venv/bin/python test_sideband.py
```

검증 범위는 [VALIDATION.md](VALIDATION.md)를 참고한다.

## 11. 선택 사항: OC waveform 입력

측정하거나 설계한 waveform을 사용할 때는 OC의
`optimalcontrol.io.export_json` 출력 형식을 읽을 수 있다.
예를 들어 decoupling 객체 전체를 다음으로 교체한다.

```json
{
  "h_larmor_mhz": 1200.0,
  "h_carrier_ppm": 8.5,
  "waveform_json": "my_oc_period.json",
  "rf_hz": 5000.0,
  "b1_scale": 1.0,
  "J_hz": 92.0
}
```

파일은 의도한 supercycle을 포함하는 **실제로 반복할 전체 주기**를 담아야
한다. Channel은 정확히 `x`/`y` 두 개, 시간축은 0에서 시작하는 균일 간격,
`metadata.pulse_dt`는 양수여야 한다. `units="a.u."`이면 amplitude를 변환할
`rf_hz`가 필요하고, `units="Hz"`이면 `rf_hz`를 생략한다.
파일별로 상속되는 설정까지 포함해 `waveform_json`과 `p90_s`/`cycle`을
동시에 지정하지 않는다. 서로 다른 waveform segment가 많으면 계산 시간이
늘어난다. `cestdec.Scheme` JSON은 별도 형식이므로 OC waveform으로 직접
입력할 수 없다.

구현 근거: [sbfit.py](sbfit.py), [sideband.py](sideband.py),
[run.py](run.py), [est_data.py](est_data.py).
