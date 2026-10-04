# SBONEST fitting 매뉴얼

[English](SBONEST_MANUAL.md) · [Step-by-step dummy 가이드](DUMMY_GUIDE.ko.md) ·
[README](README.md) · [기술 설명](SIDEBAND.md)

이 매뉴얼은 2026-10-04 기준 현재 코드의 CLI `Sideband` 모델을 설명한다. OC의 NH spin operator를 사용해 ¹H decoupling sideband를
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

처음 사용한다면 입력 복사, 새 출력 경로 생성, 결과 확인까지 포함한
[dummy 데이터 가이드](DUMMY_GUIDE.ko.md)를 따른다. 아래의 간단한 명령은
예제의 고정 출력 접두사를 사용하므로 해당 설정의 최초 실행용이다.

```bash
git clone https://github.com/Gohyang-Matzip/sbonest.git
cd sbonest
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements-sideband.txt -c constraints-sideband.txt
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg
.venv/bin/python run.py example/sideband_auto_H/two_RF.json --check
.venv/bin/python run.py example/sideband_auto_H/two_RF.json
```

OC는 `optimalcontrol-nmr` 패키지로 설치되므로 인접한 OC 폴더가 필요 없다.
이미 작동하는 환경은 그대로 사용한다. 수치 결과만 필요하면 검사와 fitting
명령 모두에 `--no-pdf`를 추가한다.

`--check`는 실제 데이터·waveform 경로, 관측점 수, 자기장·RF, 자유·고정
파라미터, 초기값, bounds와 출력 충돌을 JSON으로 표시한다. 잘못된 설정이나
충돌이 있으면 실패 종료한다. 파일을 쓰거나 optimizer를 실행하지 않지만
초기값 격자는 모델을 평가할 수 있다. 선택적 분석 설정도 fitting 전에 검사한다.
이는 준비 상태를 확인하며 수렴·식별성을 보장하지 않는다. 새 실행의 출력을
만들기 전에 사용한다.

`--check --identifiability`는 초기값에서 grouped 유한차분 Jacobian을 한 번
계산해(이 예제에서 1–2초) `identifiability` 항목에 열 정규화 특이값, rank,
조건수, 입력 absolute sigma 기준 `(J^T J)^-1`의 자유 파라미터별 기대 표준오차와
상대오차, 감도가 0인 파라미터, 약하게 결정되는 파라미터(상대오차 100% 초과
또는 유한한 오차 없음), 0.95 이상의 쌍별 상관, 유도된 kex/pB 오차를 보고한다.
이는 초기값과 샘플링 설계의 성질이므로 optimizer를 돌리기 전에 R1H/R2H가
분리되지 않는 상황을 경고할 수 있지만, fitting 후의 공분산을 대신하지는
않는다. 병렬 실행은 8.4절을 참고한다.

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
재실행하려면 JSON 사본의 `Project Name`을 새 접두사로 변경한다. JSON을 다른
폴더로 옮겼다면 데이터·waveform 경로도 변경하거나 입력 파일을 함께 복사한다.
입력 경로는 JSON 폴더 기준이고, 상대 출력 경로는 작업 폴더 기준이다.
Dummy 가이드 2단계는 제공 예제의 입력을 복사하고 새 절대 출력 접두사를 만든다.
기존 출력은 덮어쓰지 않는다. 25/50/100 Hz는
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

**두 자기장 공동 fitting.** `compare_joint_fields.py --source INPUTS --out NEW_DIR
[--workers N]`은 두 자기장의 "full" dataset 여섯 개를 12.1절의 자기장 그룹 모델로
함께 fitting해 위의 개별 fitting과 비교한다. 보관된 실행은
[`results/field_comparison_joint_20261004/summary.json`](results/field_comparison_joint_20261004/summary.json)이다.
생성에 쓴 양성자 이완을 고정하면(벤치마크 가정) 공동 fitting은 kex = 298.9 ± 0.9 s⁻¹
(참값 300)를 주며, 600 MHz 단독 297.9 ± 1.6, 800 MHz 단독 299.9 ± 1.2과 비교된다.
자기장 그룹별 양성자 이완을 fitting해도 kex·pB는 그대로지만(298.9 ± 1.0, 0.0501)
양성자 이완 자체는 이 자기장에서 결정되지 않는다. 상한이 없으면 R1H는 하한으로,
R2H는 수만 s⁻¹로 흘러가고(질소 이완까지 자기장별이면 질소 이완도 편향된다)
경계 변형은 `init.bounds` 50·500 s⁻¹(`--proton-bounds`)를 쓴다. AICc는 고정 이완
설명을 선호한다(자기장별 양성자 이완 ΔAICc +7.8, 질소 이완까지 +12.5). 추가
파라미터가 정보를 담지 않을 때 예상되는 결과다. 자기장당 잡음 한 번의 합성
참값이므로 추가 자유도의 비용을 보여 줄 뿐 실험 데이터에서의 이득을 보여 주지
않는다.

## 8. 출력 파일과 결과 해석

`Project Name = results/my_sideband_fit`이면 다음 파일이 만들어진다.

| 파일 | 내용 |
|---|---|
| `results/my_sideband_fit_result.json` | 파라미터, 전체 공분산, 파생 표준오차, 진단값, 입력 config 및 실행 이력 |
| `results/my_sideband_fit_predictions.csv` | 전체 정밀도의 관측값·예측값·σ·표준화 잔차, 잔기 및 dataset 식별자 |
| `results/my_sideband_fit_result.txt` | 사람이 읽는 결과와 관측/계산 intensity 표 |
| `results/my_sideband_fit.pdf` | 활성 잔기당 profile 및 표준화 잔차 패널 |
| `results/my_sideband_fit_data.pdf` | 데이터만 그린 profile |
| `results/my_sideband_fit_checkpoint/` | 실행 식별 정보와 완료된 fitting·분석 기록(8.3절) |

PDF는 입력된 전체 offset 범위를 보여준다. 계산선은 측정 offset에서 계산한
값을 연결한 것이며, 별도의 촘촘한 simulation이 아니다. 측정점 사이의 좁은
구조는 드러나지 않을 수 있다. Fit PDF 범례의 RF는 **fitting된 값**이며 데이터 전용 PDF는 명목값이다.
JSON의 `v1n_hz`에도 실제 RF를 기록한다. 텍스트 intensity 표는 소수 셋째
자리까지 출력하므로 정량 재분석에는 전체 정밀도의 CSV와 원래 입력을 사용한다.
CSV의 `residual_sigma`는 `(관측−예측)/σ`이며 제곱합은 χ²와 일치한다.
`plot_full_profile.py`는 별도 demo의 파일과 `validation.json`을 전제로 하므로,
임의 fitting 결과에 사용하는 범용 plotter가 아니다.

1절의 예제 실행 후, 프로젝트 폴더에서 다음으로 결과를 확인할 수 있다.

```bash
.venv/bin/python - <<'PY'
import json
from pathlib import Path
r = json.loads(Path("results/auto_H_two_RF_result.json").read_text())
assert r["success"], r.get("message")
print("kex (s^-1):", r["kex"], "SE:", r["derived_se"]["kex"])
print("pB (fraction):", r["pB"], "SE:", r["derived_se"]["pB"])
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
JSON의 `covariance`는 행과 열이 `parameter_order` 순서인 전체 행렬이며,
`derived_se.kex`와 `derived_se.pB`는 `kab`·`kba`의 공분산까지 반영한다.
계산할 수 없는 원소·표준오차는 `null`로 기록한다. `schema_version=2`는 기존
scalar 필드를 유지하고 이 정보를 추가한다. `provenance`에는 fitting 전 입력·
waveform·실행 소스의 SHA-256, canonical config 해시, Python/package 버전,
플랫폼·thread 설정·Git 상태와 실행 시간을 저장한다. 해시는 원본 보관을 대신하지 않는다.
`parameter_order`는 고정 파라미터까지 포함하지만 `n_parameters`는 자유
파라미터만 센다. `pB`와 `derived_se.pB`는 분율이므로 각각 100을 곱해
population은 %, SE는 percentage point로 보고한다. 이전 보존 결과에는 새
필드가 없을 수 있다. 원본은 유지하고 필요하면 현재 코드로 새 fitting을 수행한다.
`constraints-sideband.txt`는 검증한 Python 3.12 의존성 조합이다.
선택적 profile likelihood와 parametric bootstrap은 8.1–8.2절에서 설명한다.
두 방법 모두 모델·고정 입력의 불확실성을 제거하지 않는다.

## 8.1. 다중 초기값, profile likelihood 및 성능 검사

아래 선택 항목은 `init` 안에 넣는다. 생략하면 기존의 단일 fitting을 수행한다.
이 블록은 기존 `init`에 합칠 부분이며 독립적인 전체 설정 파일이 아니다.
[Dummy 가이드](DUMMY_GUIDE.ko.md)의 6단계는 JSON을 직접 편집하지 않고도
완전한 설정 파일을 만들고 실행하는 명령을 제공한다.

```json
"multistart": {
  "starts": [{"kab": 20.0, "kba": 380.0, "v1n_scale": 1.0}],
  "random_starts": 2,
  "seed": 20261004
},
"profile": {
  "kex": [295.0, 300.0, 305.0],
  "pB": [0.045, 0.05, 0.055],
  "v1n_scale": [1.06, 1.08, 1.10]
}
```

`multistart`는 설정된 초기 fitting을 항상 포함한다. `starts`는 6절의
파라미터 이름으로 지정한 초기값 override 목록이다. 난수 초기값에는 명시적
seed가 필요하다. 자유 파라미터만 바꿀 수 있으며 고정 파라미터는 유지한다.
유한한 양쪽 bounds 안에서는 균등분포를, 한쪽/무한 bounds에서는
`sb_analysis.py`에 명시된 국소 섭동을 사용한다. 초기값 민감도를 살피는
절차이며 전역 최적해를 입증하지 않는다. 모든 시도의 전체 초기값·결과 벡터,
상태·메시지·χ²를 `multistart`에 기록하고, 수렴한 시도 중 최저 χ²를
`selected=true`로 표시한다. 위 예제는 기본 1회, 명시적 초기값 1회, 난수
초기값 2회로 총 4회 시도한다. 수렴 상태가 서로 다를 수 있어 각 기록을 확인한다.
전부 실패하면 CLI는 실패 종료하고
`success=false`인 `_result.json`에 시도와 실행 이력을 보존한다.
재시도에는 새 출력 접두사를 사용한다. 잘못된 설정은 fitting 시도로 취급하지 않고 거부한다.

`profile`은 kex, pB 또는 v1n_scale(scale 모드만)의 명시적 격자를 받는다.
각 점에서 나머지 자유 파라미터를 다시 최적화하며 물리적·사용자 bounds와
고정 파라미터를 지킨다. kex/pB는 kab/kba를 정확히 변환하여 제약하며 penalty
잔차를 추가하지 않는다. 불가능하거나 수렴하지 않은 점은 target·실패 메시지와
null χ²를 남긴다. `profiles`에는 기준 χ², 각 결과 벡터, nuisance 이름과
원래 Δχ²를 저장한다. 음수 Δχ²는 0으로 바꾸지 않고 경고한다. 이는 기준보다
더 좋은 해를 찾았다는 뜻이며, 기본 보고 결과를 몰래 교체하지 않는다.
각 점은 선택된 기본 해에서 출발하므로 국소 수렴이 제약하의 전역 최적해를
보장하지 않는다. 경계·약한 식별성·격자 범위·모델 적합성을 검토해야 하므로
격자만으로 confidence interval을 자동 산출하지 않는다. 이는 절대 σ를
사용하는 모델 내 likelihood 진단이며 실험 검증이나 고정 입력의 불확도는
포함하지 않는다. Parametric bootstrap은 8.2절을 참고한다.

최상위 `success`는 선택된 fitting의 성공 여부이며 모든 scan 점의 성공을
뜻하지 않는다. `profiles.<name>`의 각 항목에 있는 `success`, `message`와
profile 경고를 읽는다. Multistart 벡터의 순서는 `parameter_order`, profile
벡터의 순서는 `profiles.parameter_names`에 기록된다.

기본 fitting과 제약 profile fitting은 임의의 `init.vary` 부분집합·순서에
묶음 수치미분을 사용하며 bounds에서는 가능한 안쪽 방향으로 미분한다. 잔기별로
독립인 좌표만 묶으며 교환 제약의 결합된 영향은 유지한다. 벤치마크는
Sideband·ONEST config를 모두 받아 단일 fitting 시간을 측정한다(선택적 추가 분석은 실행하지 않는다).

```bash
.venv/bin/python benchmark.py example/sideband_auto_H/two_RF.json
mkdir -p session_artifacts
.venv/bin/python benchmark.py example/sideband_auto_H/two_RF.json profile \
  --profile-output session_artifacts/sideband_01.prof
```

상위 폴더를 먼저 만든다. 기존 profile 파일은 덮어쓰지 않는다. 예전의
`config.json profile` 형식도 유지하며 경로가 비어 있을 때
`benchmark_profile.prof`를 저장한다. Python `pstats`로 읽을 수 있다.
벤치마크는 일반 fitting의 JSON/CSV/PDF를 생성하거나 덮어쓰지 않는다.
Profiling 시간에는 profiler와 보고 출력의 비용이 포함되므로 일반 fitting
시간과 직접 비교하지 않는다.

## 8.2. Parametric bootstrap과 합성 불확실성 검사

새 출력 접두사를 지정한 설정 파일의 `init` 안에 다음 부분을 추가한다.

```json
"bootstrap": {"replicates": 5, "seed": 20261004, "confidence": 0.95}
```

`replicates`는 양의 정수, `seed`는 명시적인 0 이상의 정수, `confidence`는
0과 1 사이 값이어야 한다(기본 0.95). 각 반복에서 선택된 기본 해의 예측값에
입력한 절대 sigma를 사용한 독립 Gaussian 잡음을 더한다. Offset·자기장·RF와
고정 입력은 유지한다. 각 합성 데이터를 선택된 기본 해에서 다시 fitting하며
반복 안에서 restart나 profile scan을 추가로 실행하지 않는다. 이후 원래 관측값과
선택된 fitting 상태를 복원한다.

`bootstrap.samples`는 실패까지 포함해 각 index·상태·메시지·fitting 벡터·
국소 SE·χ²·경계 정보를 보존한다. 벡터 순서는 `parameter_names`다.
`bootstrap.intervals`는 성공한 반복의 lower/median/upper percentile과
`n_success`·`fixed`를 각 모델 파라미터와 파생 kex/pB에 대해 기록한다.
고정 구간은 가정을 나타낸다. `successful`, `replicates`, `warnings`로 실패를
확인하며, 실패를 제외하면 분포에 편향이 생길 수 있다. 성공한 반복이 100회
미만이면 꼬리 추정이 불안정하다는 경고가 나온다. 5회는 실행 절차의 시연이다.
반복 횟수가 많아도 명목 신뢰수준의 coverage를 보장하지 않는다. 구간은 선택
모델·입력한 절대 sigma·고정 입력에 조건부이며 모델 불일치를 포함하지 않는다.

별도의 반복 데이터 연구에는 `.venv/bin/python validate_uncertainty.py --config PATH --truth PATH
--replicates N --seed N --confidence 0.95 --out NEW_DIRECTORY`와 직접 준비한
합성 설정·알려진 생성값을 사용한다. Truth JSON은 **모든** 모델 파라미터 이름을
유한한 값에 대응한 map 또는 `{"truth": {...}}` 형식이어야 하며 추가 이름은
허용하지 않는다. 값은 설정 bounds 안에 있어야 하며 고정 파라미터는 설정된
초기값과 일치해야 한다. `--check`로 파라미터 이름을 확인할 수 있다. Fitting
추정값은 독립적으로 알려진 생성값이 아니다.

연구는 생성값에서 새 Gaussian 관측값을 만들고 매번 생성값을 시작점으로
fitting한다. `coverage.json`은 **국소 normal-SE 구간의 coverage**를 측정하며
bootstrap percentile 구간의 coverage를 측정하지 않는다. 자유 변수별로 유효
횟수·coverage 비율·binomial SE·Wilson 95% 구간·산출 불가/0인 SE 횟수·경계
빈도를 보고하며 실패도 별도로 기록한다. Coverage 분모에서 실패와 양수가
아니거나 산출 불가인 SE를 제외하므로 제외 내역을 확인한다. 고정 변수에는
coverage를 주장하지 않는다. 출력에는 `study.json`, config/truth 사본,
`samples.json`, 개별 `samples/` 기록도 보존한다. 작은 연구는 표본 불확실성이
크며 실험 타당성이나 보편적인 coverage를 입증하지 않는다.

`--profile-interval kex pB`와 `--inner-bootstrap B`는 연구를 확장한다. 각 replicate가
자신의 데이터에서 likelihood-ratio 구간(8.5절)과 B회 parametric bootstrap(seed
`seed·1000003 + index`)도 받고, `coverage.json`에 profile·bootstrap 구간의 coverage
비율·binomial SE·Wilson 구간·평균 폭을 담은 `intervals` 블록이 국소 SE coverage
옆에 추가된다. 열린(한쪽) profile 구간과 실패한 내부 bootstrap은 따로 세고 비율에서
제외한다. replicate마다 fitting 한 번 + 구간 양당 약 여섯 번의 재fitting + B번의
fitting이 들므로 `--workers N`을 쓴다(replicate를 병렬로 돌리며 직렬 결과와
동일하다). 100회 미만은 여전히 시연이다.

## 8.3. Checkpoint, 재개, 저장 결과 보고서

일반 fitting은 기본적으로 `PROJECT_checkpoint/`를 만든다. 변경하지 않는
실행 식별 정보 `manifest.json`, `provenance.json`, 완료된 기본 fitting의
`baseline.json` snapshot을 보존한다. 선택적 분석은 `attempt-0.json`,
`profile-kex-0.json`, `bootstrap-0.json`과 이후 번호의 기록을 추가한다.
실패한 fitting은 가능한 경우 `failure.json`과 실패 결과 JSON에 진단을 남긴다.
출력 준비 폴더(`export-*`), `exports.json`, `complete.json`은 출력 파일의
저장을 추적한다. 각 단계가 완료될 때 파일이 생기므로 폴더 전체를 결과와
함께 보관한다. Checkpoint 기록은 checksum으로 검증하므로 편집하지 않는다.
변경된 기록으로는 재개할 수 없다.

1절의 예제를 중단했다면 다음처럼 재개한다.

```bash
.venv/bin/python run.py example/sideband_auto_H/two_RF.json --resume
```

완료된 기본 fitting·restart·profile·bootstrap 작업은 재사용하고 Ctrl-C 등으로
중단된 미완료 작업을 실행한다. 모든 초기값이 실패한 multistart를 포함해 이미
실패로 끝난 fitting은 실패 JSON을 복원하고 재시도 없이 실패 종료한다. 실패를
확인하고 설정을 수정한 뒤 새 `Project Name`으로 다시 시도한다. 재개하려면
설정, 실제 입력·waveform 경로와 해시, 실행 소스 해시,
Python·package 버전, 플랫폼, thread 설정, `--no-pdf` 모드가 일치해야 한다.
처음 `--no-pdf`를 사용했다면 재개할 때도 붙인다. `--check`와 `--resume`는
함께 사용할 수 없으며 기존 실행을 검사하면 보호된 경로가 충돌로 표시된다.
설정·소스·실행환경을 바꾸려면 새 접두사가 필요하다. 재개는 관계없거나 변경된
출력을 덮어쓰지 않는다. 일반 재실행도 실패한 실행을 포함해 기존 checkpoint가
있으면 거부한다. 설정을 바꿔 새로 시작하기 전에 실패 기록을 보존한다.

성공한 저장 결과와 predictions CSV에서 보고서를 다시 만든다.

```bash
.venv/bin/python sb_workflow.py report results/auto_H_two_RF_result.json \
  --out results/auto_H_two_RF_report_01
```

이 명령은 저장 JSON/CSV만으로 `_summary.json`, `_summary.txt`, `.pdf`를
만들며 최적화나 원래 데이터 파일을 요구하지 않는다. 관측점 수·χ² 등의
일치 여부를 검사하고 잔기·dataset별 잔차 통계, fitting 진단, 저장된 restart·
profile·bootstrap 기록을 담는다. 실패한 profile 점과 음수 Δχ²도 유지한다.
파일명을 바꿨다면 `--predictions PATH`로 CSV를 지정한다. 실패 기록만 있는
JSON에는 보고할 predictions가 없다. 기존 파일은 보호하므로 매번 새 보고서
접두사를 사용한다. 원래 fitting을 `--no-pdf`로 실행했어도 PDF를 만들 수 있다.

## 8.4. `--workers`로 병렬 실행하기

`run.py CONFIG --workers N`(`sbfit.py`와 `--check --identifiability`에서도 사용
가능)은 같은 설정과 데이터를 가진 모델 사본을 하나씩 든 worker 프로세스 N개를
시작한다. 주 프로세스가 optimizer와 checkpoint를 담당하고, worker는 모든
fitting의 grouped Jacobian 열, 명시적·무작위 restart, profile 점, bootstrap
replicate를 계산한다. 결과는 index 순서로 모으므로 checkpoint 기록은 그대로
순서 있는 prefix이며 `--resume`도 직렬 실행과 똑같이 동작한다. 결과가 worker
수에 의존하지 않으므로 worker 수는 checkpoint 식별 정보에 포함하지 않는다.
병렬 실행은 직렬 실행의 모든 수치를 재현하며, `test_sb_parallel.py`가 결과
JSON과 checkpoint 기록에서 이를 검증한다.

`OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS`는 1로
유지한다. 병렬성은 pool이 제공한다. 10코어 노트북에서 제공 예제의 882점
fitting은 worker 1개로 약 31초, 8개로 10초가 걸렸다. pool을 열 때 모든 worker를
미리 시작하고, Jacobian 열은 worker당 벡터 하나씩 계산하며, optimizer 자체의
잔차 계산은 잔기/dataset 블록으로 나눠 worker에 분산한다. restart·profile 점·
bootstrap replicate는 항목 수까지 worker 수에 거의 비례해 빨라진다. worker는
spawn 방식으로 시작하므로
`run_config(..., workers=N)`를 호출하는 Python 스크립트는 진입점을
`if __name__ == "__main__":`로 감싸야 한다. restart·profile·bootstrap에는 경과
시간과 남은 시간 추정이 담긴 진행 표시가 출력된다. Ctrl-C 후 worker는 진행
중인 항목을 마치고 종료하며, 완료된 항목은 이미 checkpoint에 있다.

## 8.5. Profile에서 구하는 likelihood-ratio 구간

`init.profile`(8.1절)은 정해진 격자를 계산한다. `init.profile_interval`은
정확한 nuisance 재fitting profile이 chi-square 임계값과 만나는 두 점을 찾아
likelihood가 2차식이라는 가정 없이 모델 내 신뢰구간을 준다.

```json
"profile_interval": {
  "parameters": ["kex", "pB"],
  "confidence": 0.95,
  "max_evaluations": 40,
  "relative_tolerance": 0.001,
  "max_doublings": 8
}
```

`parameters`에는 `kex`, `pB`, scale RF 모드에서는 `v1n_scale`을 쓸 수 있다.
임계값은 자유도 1의 chi-square 분위수(95%에서 3.84)다. 추정값에서 양쪽으로
`z × 국소 SE`(유한한 국소 오차가 없으면 추정값의 10%)에서 시작해 두 배씩
넓혀 가며 교차점을 가두고, Brent 방법으로 `relative_tolerance × |추정값|`
정밀도까지 찾는다. 계산한 모든 profile 점은 계산 순서대로 결과의
`profile_intervals.<name>.points`와 checkpoint의 `profile_interval-<name>-N`
기록에 남으므로 `--resume`는 재fitting 없이 탐색을 재생한다. 파라미터 bound나
`max_doublings` 안에서 profile이 임계값에 닿지 않으면 그쪽은 숫자 대신
메시지와 함께 열린 구간으로 보고한다. 기준 chi2보다 낮은 점이 나오면 기준
fit이 최적이 아니며 구간을 신뢰할 수 없다고 표시한다. nuisance 재fitting이
실패하면 `message`에 기록하고 나머지 파라미터는 계속 계산한다. 계산 한 번은
제약 재fitting 한 번이며(제공 예제에서 직렬 10–17초; `--workers`로 Jacobian
열을 병렬화한다) 보고서에는 구간, 계산한 모든 점, 임계값을 함께 그린 profile
그림이 들어간다. 이 구간은 모델·고정 입력·입력 absolute sigma에 조건부이며
전역 최적을 보장하지 않는다.

## 8.6. 실험 설계: 측정 전에 기대 오차 계산하기

`python sb_workflow.py design DESIGN_JSON --out NEW_DIRECTORY [--workers N]`은
알려진 참값에서 모델과 국소 Fisher 정보만으로 계획한 측정이 파라미터를 얼마나
잘 결정하는지 평가한다. optimizer는 실행하지 않는다. 설계 파일은 기준 설정
(decoupling, 잔기, RF·proton 설정, bounds, `vary`를 재사용), 참값, 시나리오를
담는다.

```json
{
  "config": "fit.json",
  "truth_result": "results/auto_H_two_RF_result.json",
  "scenarios": [
    {"name": "two_rf_147", "datasets": [
      {"v1n_hz": 25, "T": 0.4, "sigma": 0.01, "offsets_ppm": {"min": 105, "max": 135, "n": 147}},
      {"v1n_hz": 100, "T": 0.4, "sigma": 0.01, "offsets_ppm": {"min": 105, "max": 135, "n": 147}}]},
    {"name": "three_rf_60", "datasets": [
      {"v1n_hz": 25, "T": 0.4, "sigma": 0.01, "offsets_rel_ppm": {"min": -15, "max": 15, "n": 60}},
      {"v1n_hz": 50, "T": 0.4, "sigma": 0.01, "offsets_rel_ppm": {"min": -15, "max": 15, "n": 60}},
      {"v1n_hz": 100, "T": 0.4, "sigma": 0.01, "offsets_rel_ppm": {"min": -15, "max": 15, "n": 60}}]}
  ]
}
```

`truth`는 모든 모델 파라미터 이름에 값을 대응시키고, 대신 `truth_result`로
저장된 결과 JSON의 fitting 값을 참값으로 쓸 수 있다(fitting 추정값은 계획용
가정이지 독립적으로 알려진 값이 아니다). 시나리오의 각 dataset은 질소 RF
세기, saturation 시간, absolute sigma, offset 격자(절대 `offsets_ppm` 또는
각 잔기 `peak_ppm` 기준 상대 `offsets_rel_ppm`)를 주며 `v1err_hz`, `field_mhz`,
`decoupling` 덮어쓰기는 선택이다. 시나리오마다 잡음 없는 완전한 합성 입력과
`design_config.json`을 `scenarios/<name>/`에 쓰고, 참값에서 grouped Jacobian을
계산해 `design.json`, `design.txt`, `design.pdf`에 관측점 수, 측정 시간 대리값
Σ(점 수 × T), rank와 조건수, 자유 파라미터별 기대 표준오차와 상대오차, 유도
kex/pB 오차, 약하게 결정되는 파라미터, 강한 상관을 보고한다. 기대 오차는
sigma에 정확히 비례하는 참값에서의 국소 선형값이며 모델 안에서 설계를
비교할 뿐 시료를 검증하지 않는다. 제공 예제에서 sigma 0.01일 때 RF 두 세기와
offset 147개는 kex 기대 오차 6.7 s⁻¹, 세기당 60개는 saturation 시간 41%로
10.8 s⁻¹, 100 Hz 한 세기만 쓰면 116 s⁻¹이다.

**어디를 측정할지 최적화하기.** 선택 항목 `optimize`는 한 시나리오의 조밀한
후보 격자에서 측정 예산만큼을 고른다.

```json
"optimize": {"scenario": "two_rf_147", "budget": 60, "criterion": "kex", "min_per_dataset": 4}
```

후보 시나리오는 절대 `offsets_ppm` 격자여야 한다. saturation offset 하나가
모든 잔기를 한꺼번에 주는 스펙트럼 행 하나이므로 선택 단위는 행 전체다. 모든
후보 행에서 시작해 제거해도 기준값이 가장 적게 나빠지는 행을 Fisher 정보의
정확한 Woodbury 갱신으로 반복 제거하며, `budget`개 행이 남고 dataset마다 최소
`min_per_dataset`개가 유지될 때 멈춘다. `criterion`은 `kex`(기본)·`pB`·자유
파라미터 이름의 기대 분산 또는 행렬식 기준 `D`다. 결과에는 보통 시나리오와
똑같이 평가한 `<이름>_optimized_<budget>`와 `<이름>_uniform_<budget>`(같은
예산을 균일 배치) 두 시나리오, dataset별 선택 offset, 기준값 경로와 그림이
추가된다. 제공 예제에서 sigma 0.01일 때 후보 294행(RF 두 세기) 중 60행은 kex
기대 오차 7.6 s⁻¹로, 균일 60행의 17.5 s⁻¹, 전체 294행의 6.7 s⁻¹과 비교된다.
선택은 참값과 모델에 국소적인 계획 보조 도구이며, 한 양에 최적인 설계가 다른
양에는 더 나쁠 수 있다.

## 8.7. 공유 교환 대 잔기별 교환 비교

`python sb_workflow.py compare CONFIG --out NEW_DIRECTORY [--workers N] [--pdf]`는
설정한 모델을 두 방식으로 fitting한다. 하나는 설정대로 모든 활성 잔기가
kab/kba(와 RF scale)를 공유하는 전역 모델이고, 다른 하나는 다른 잔기를 모두
끄고 잔기별로 따로 fitting한 개별 모델이다. 잔기별 `initial`, `bounds`,
`vary`, multistart 시작값은 해당 잔기의 것만 유지하고 profile·bootstrap·구간
분석은 반복하지 않는다. 각 하위 fitting은 `global/`과 `individual/<잔기>/`
아래에서 보통의 checkpoint 경로로 실행된다. `comparison.json`,
`comparison.txt`, `comparison.pdf`에는 모델별 chi2, 파라미터 수, 국소 오차가
붙은 kex·pB, 개별 chi2 합, 두 설명의 AICc·BIC(입력 absolute sigma의 Gaussian,
공통 상수 제외), AICc 기준 선호 모델, 공유 모델 대 잔기별 모델의 중첩 F-검정이
담긴다. 하위 fitting이 실패하면 보고하고 비교는 비워 둔다. 이 통계는 입력
sigma와 고정 입력 아래에서 설명을 비교할 뿐이다. 공유 속도가 선호되면 하나의
교환 과정과 부합하지만 증명은 아니며, 개별 속도가 선호되는 것은 모델 불일치나
잘못 보정된 오차 때문일 수도 있다.

**2상태인가 3상태인가.** `sb_workflow.py compare CONFIG --out DIR --models Sideband
Sideband_3st_Linear [--h-ppm-c A1=7.1 ...]`은 같은 데이터를 2상태 설정과 거기서
유도한 3상태 모델(12.2절)로 fitting한다. 유도 설정은 `--h-ppm-c`가 없으면 C 상태에
B 상태의 양성자 이동을 재사용하고, `kbc`/`kcb`와 `dwC_ppm`을 2상태 교환 속도
주변의 명시적 multistart 조합 다섯 개에서 시작하며, 2상태 분석은 제외한다.
`comparison.json`/`.txt`에는 모델별 chi²·AICc·BIC·분포·속도·경계 표시와 2상태
대비 AICc/BIC 차이가 담긴다. 3상태 모델은 파라미터 공간의 경계에서만 2상태로
환원되므로 F-검정은 보고하지 않는다. 복귀 속도(`kcb`)가 하한에 닿으면 C 상태가
흡수 상태가 되므로 경고가 그런 fitting을 퇴화로 표시하고 분포는 의미가 없다.
합성 3상태 데이터에서는 ΔAICc ≈ −26,000으로 3상태가 선호되고 속도가 복원되며,
2상태 데이터에서는 3상태 fitting이 퇴화하고 AICc가 2상태를 선호한다(+6.9).

## 8.8. 잔차·sigma 진단

이 매뉴얼의 모든 불확실성 진술은 입력 absolute sigma와 데이터를 설명하는 모델을
가정한다. fitting이 끝나면 결과 JSON의 `residual_diagnostics`가 표준화 잔차만으로
다음을 계산한다.

- 기대 산포 `sqrt(2/dof)`가 붙은 `reduced_chi2`와 그로부터 유도한
  `sigma_scale_estimate = sqrt(chi2/dof)`. 1에서 산포의 세 배 넘게 벗어나면 sigma가
  너무 작거나(또는 모델이 구조를 놓치거나) 너무 크다고 경고한다.
- 잔기별·dataset별·잔기/dataset 블록별 점 수와 reduced chi-square. 블록마다
  offset 순으로 정렬한 잔차 부호의 Wald–Wolfowitz runs test, `2/sqrt(n)` 기준이
  붙은 lag-1 자기상관, 최대 표준화 잔차. 체계적인 부호 run이나 자기상관은
  sigma가 아니라 모델 부적합을 경고한다.
- Gaussian 기대치와 비교한 3 sigma 초과 이상점 수.
- `rescaled_stderr`와 자유 파라미터별 `stderr_rescaled`: 국소 표준오차에
  `sqrt(chi2/dof)`를 곱한 값. sigma가 균일하게 잘못 추정되었다고 믿을 때 쓰는
  관례적 대안이며 모델이 맞다고 가정한다. 기본값이 아니다.

텍스트 보고서와 `sb_workflow.py report`도 같은 진단을 출력하고(보고서는 predictions
CSV에서 다시 계산한다) 그 경고는 결과 경고에 합쳐진다. 이 검정들은 지표일
뿐이다. 1에서 먼 reduced chi-square는 sigma와 모델 중 무엇이 틀렸는지 말해 주지
않으며, runs test는 짧은 블록에서 검정력이 떨어진다.

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

여러 데이터셋을 함께 fitting하면 각 잔기의 nitrogen shift를 모든 파일에서
공유하고, 기본값에서는 `R1`, `R2a`, `R2b`도 공유한다. 자기장 그룹, 자기장별
양성자·질소 이완, 3상태 모델은 12절에 설명한다. 자기장을 합치기 전에 공유
fitting과 자기장별 fitting을 비교한다.

모델은 교환하는 각 상태(기본 2개, 12.2절의 3상태 모델은 3개)에 N 하나와 H 하나를 포함하고, 현상론적
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
| `Output already exists` | 새 `Project Name` 사용. 변경하지 않은 중단 실행은 `--resume` |
| Checkpoint identity mismatch | 원래 설정·실행환경을 복원해 재개하거나 기존 실행을 보존하고 새 접두사 사용 |
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
| `No multistart attempt converged` | 실패 결과 JSON의 상태·메시지를 확인한 뒤 타당한 초기값·bounds와 새 접두사로 재시도 |
| Profile Δχ²가 음수 | Scan이 기준 해를 개선함. 해당 벡터를 확인하고 새 출력 접두사로 그 해에서 다시 fitting |
| `derived_se`, `covariance`, `provenance` 없음 | 이전 보존 결과인지 확인. 원본을 유지하고 현재 코드로 새 fitting 수행 |
| `Profile already exists` | 기존 측정은 보존하고 새 `--profile-output` 파일명 지정 |
| 실행이 오래 걸림 | 수치 라이브러리 thread를 1로 유지하고 `--no-pdf`와 `--workers N`(8.4절) 사용. 0이 아닌 `v1err`는 RF 평균화 비용 추가 |

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

## 12. 여러 양성자 자기장과 3상태 교환

### 12.1. 자기장 그룹과 자기장별 이완

decoupling 항목의 `h_larmor_mhz`가 같은 dataset들은 하나의 **자기장 그룹**이 된다.
그룹은 자기장 오름차순으로 번호가 붙고 결과(`field_groups`)와 텍스트 보고서에
표시된다. 그룹이 하나면 달라지는 것이 없다. 그룹이 여럿이면:

- 자동 양성자 이완(`proton_relaxation.mode = "fit"`)은 잔기마다 **그룹별로**
  `R1H`/`R2H` 한 쌍을 fitting하며 이름은 `A1.R1H[0]`, `A1.R2H[0]`, `A1.R1H[1]`,
  …이다. 따라서 자동 양성자 이완으로 여러 자기장을 함께 fitting할 수 있고 이전의
  단일 자기장 제한은 없어졌다.
- `sideband` 절의 `"nitrogen_relaxation": {"mode": "per_field"}`는 잔기마다
  그룹별 `R1`, `R2a`, `R2b`(와 `R2c`)를 두며 이름은 `A1.R1[0]`, `A1.R2a[1]`,
  …이다. 기본값 `"shared"`는 이전처럼 모든 자기장에 한 벌을 쓴다. peak 위치,
  화학적 이동 차이, 교환 속도, RF 파라미터는 항상 공유한다.
- `init.initial`, `init.bounds`, `init.vary`, multistart 시작값에서 `A1.R1`처럼
  그룹 표시가 없는 이름은 그 키의 모든 그룹을 가리킨다. 그룹 이름을 직접 쓰면
  그 값이 우선한다.

```json
"sideband": {
  "decoupling": {"h_carrier_ppm": 8.5, "p90_s": 7e-05},
  "datasets": [{"h_larmor_mhz": 600}, {"h_larmor_mhz": 600},
               {"h_larmor_mhz": 800}, {"h_larmor_mhz": 800}],
  "proton_relaxation": {"mode": "fit"},
  "nitrogen_relaxation": {"mode": "per_field"},
  "residues": {"A1": {"h_ppm_a": 6.2, "h_ppm_b": 6.5}},
  "v1n": {"mode": "scale"}
}
```

자기장별 이완은 파라미터를 늘리므로 `--check --identifiability`로 식별성을
확인하고, 공유·자기장별 fitting을 서로 다른 출력 접두사로 비교한 뒤 결론을
내린다. 7.1절은 이 옵션으로 실행한 600/800 MHz 공동 벤치마크를 보고한다. kex·pB는
그대로이고, 양성자 이완은 그 자기장에서 결정되지 않아 bounds가 필요하며, AICc는
고정 이완 설명을 선호한다.

### 12.2. 3상태 교환

`init.Method = "Sideband_3st_Linear"`(A ⇄ B ⇄ C) 또는 `"Sideband_3st_Triangle"`
(A ⇄ C 추가)은 2상태 모델과 같은 decoupling·이완·RF 처리로 48차원 NH
Liouvillian을 전파한다. 속도는 `kab, kba, kbc, kcb`(삼각형은 `kca, kac` 추가)이고
잔기마다 `dwC_ppm`(A 기준 C 상태의 이동)과 `R2c`가 더해지며 `sideband.residues`
항목에 `h_ppm_c`가 필요하다. 분포는 속도 네트워크의 정상 분포에서 구해
`exchange.populations`에 보고한다. `kex`와 `pB`는 2상태 모델만 요약하므로
null이며 `exchange.kex_AB`/`kex_BC`가 쌍별 합을 준다. 3상태에는 교환 속도
격자 탐색이 없으므로 `init.initial`에 시작 속도를 준다. 3상태 모델에서는
`kex`/`pB`의 profile과 profile 구간을 거부하며, 대신 bootstrap percentile을
모든 속도에 대해 보고한다. 3상태 minor-state fitting은 2상태보다 많은 RF 세기나
자기장을 요구하고 시작값에 민감하므로 `multistart`를 쓰고 `--check
--identifiability`를 확인하며, 0에 가까운 분포는 데이터가 지지하지 않는 것으로
본다. C 상태가 비어 있도록 `kbc`, `kcb`를 두면 2상태 결과를 재현하며,
`test_sb_models.py`가 이 사실과 합성 선형 3상태 참값의 복원을 검증한다.

```json
"init": {
  "Method": "Sideband_3st_Linear",
  "initial": {"kab": 12, "kba": 300, "kbc": 80, "kcb": 60,
              "A1.dwC_ppm": -3.5, "A1.R2c": 18}
},
"sideband": {"residues": {"A1": {"h_ppm_a": 6.2, "h_ppm_b": 6.5, "h_ppm_c": 7.1}}, ...}
```

## 13. `sbonest` 명령 설치

저장소는 설치 가능한 패키지이기도 하다. Python 3.12 환경에서:

```bash
python -m pip install -e . -c constraints-sideband.txt
sbonest --help
sbonest version
```

콘솔 명령은 모든 도구를 묶는다. `sbonest check CONFIG [--identifiability]`,
`sbonest fit CONFIG [--no-pdf] [--workers N]`, `sbonest resume CONFIG`,
`sbonest report RESULT_JSON --out PREFIX`, `sbonest init-demo --out DIR`,
`sbonest design DESIGN_JSON --out DIR`, `sbonest compare CONFIG --out DIR`,
`sbonest import-bruker ...`(14절), `sbonest serve`(15절),
`sbonest benchmark CONFIG [profile]`, 그리고 패키지 버전과 provenance에 기록되는
실행 소스 해시를 출력하는 `sbonest version`이다. 모든 명령은 스크립트(`run.py`,
`sb_workflow.py` 등)와 같은 함수를 호출하므로 출력·checkpoint·provenance가
동일하며, 설치 없이 checkout에서 스크립트를 그대로 써도 된다. 모듈은 저장소
루트에 평평하게 유지되어 `provenance`의 소스 해시가 의미를 잃지 않는다.

## 14. Bruker pseudo-2D 데이터 가져오기

`python sb_import.py PDATA OFFSETS --out NEW_FILE --peak LABEL=ppm[:dw_ppm] ...`
(또는 `sbonest import-bruker ...`)은 처리된 Bruker pseudo-2D CEST 실험을 SBONEST
dataset 파일 하나로 변환한다. `PDATA`는 `procs`, `proc2s`, `2rr`가 있는 처리
폴더(`.../pdata/1`)이며, 판독기는 NumPy만으로 submatrix 배치, 두 byte order,
int32·float64 저장, `NC_proc` 스케일을 처리한다. 검출 차원은 양성자 축(`OFFSET`,
`SW_p`, `SF`에서 ppm)이고 각 행은 하나의 saturation offset이다.

필수 입력은 명시적이다. `OFFSETS`는 행마다 하나의 saturation offset(포화 핵의
ppm, 또는 `--offset-unit hz --carrier-ppm X`와 함께 Hz)을 나열한다.
`--peak A1=8.30:2.5`는 잔기의 양성자 위치와 선택적 시작 `dw`를 준다(반복 가능).
`--reference-row K`는 세기를 정규화하고 출력에서 제외되는 기준 스펙트럼이고,
`--noise-region LO HI`는 기준 행에서 표준편차를 구해 absolute error로 쓰는 신호
없는 양성자 구간이며, `--saturation-s`와 `--v1-hz`가 헤더를 채운다.
`--half-width`(ppm)와 `--mode max|sum`은 창 통계를, `--r2a`, `--r2b`는 헤더
시작값을, `--exclude-row`는 추가 제외 행을 정하고, 포화 핵 자기장은 기본적으로
`SF × γ(15N)/γ(1H)`이다(`--nucleus`, `--field-mhz`). 출력은 2절의 텍스트 형식이며
JSON 요약이 출력된다. 위상·baseline 품질, peak 겹침, 기준 행 선택은 사용자의
책임이므로 fitting 전에 변환된 profile을 확인한다.

## 15. Sideband fitting 웹 실행기

`python sb_server.py [--host 127.0.0.1 --port 5050]`(또는 `sbonest serve`)은
터미널보다 브라우저를 선호하는 사용자를 위한 Flask 페이지를 띄운다. Sideband
설정과 데이터 파일(선택적으로 OC waveform)을 올리면 작업이 만들어진다. 서버는
`datasets`를 올린 파일 이름으로, 출력 접두사를 `fit`으로 바꾸고 모든 것을
`SB_JOBS/<작업 id>/`(또는 `SBONEST_JOBS_DIR`)에 저장한 뒤 식별성 진단을 포함한
사전 검사를 실행해 JSON을 보여 준다. 버튼으로 선택한 worker 수의 백그라운드
`run.py` 프로세스로 fitting을 시작하고, checkpoint로 중단된 작업을 재개하며,
보고서를 다시 만든다. 페이지는 작업 상태(프로세스 상태, checkpoint 기록, 결과
요약, 로그 끝부분)를 주기적으로 조회하고 모든 출력의 다운로드 링크를 보여 준다.
작업 폴더 안의 파일만 제공한다. 실행기는 신뢰할 수 있는 로컬 네트워크용이며
인증이 없다.

구현 근거: [sbfit.py](sbfit.py), [sideband.py](sideband.py),
[run.py](run.py), [est_data.py](est_data.py), [sb_analysis.py](sb_analysis.py),
[sb_checkpoint.py](sb_checkpoint.py), [sb_workflow.py](sb_workflow.py),
[sb_report.py](sb_report.py), [sb_bootstrap.py](sb_bootstrap.py),
[sb_design.py](sb_design.py), [sb_compare.py](sb_compare.py), [sb_import.py](sb_import.py),
[sb_server.py](sb_server.py), [sb_cli.py](sb_cli.py), [validate_uncertainty.py](validate_uncertainty.py).
