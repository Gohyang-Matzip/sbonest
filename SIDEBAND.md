# SBONEST: sideband를 포함한 ¹⁵N CEST fitting

`ONEST`의 복사본에 `Sideband` 모델을 추가했다. 기존 ONEST 입력 파일을 사용하며,
main/minor dip과 ¹H decoupling sideband를 모두 동일한 residual에 포함한다.
원본 ONEST와 OC의 소스는 변경하지 않았다. 새 모델은 CLI에서 사용한다.

## 실행

새 checkout에서 실행환경을 설치한 뒤, 이식 가능한 두 RF 예제를 실행한다.
OC는 PyPI의 `optimalcontrol-nmr` 패키지로 설치되며 인접한 OC 폴더는 필요 없다.

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-sideband.txt
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg
.venv/bin/python run.py example/sideband_auto_H/two_RF.json
```

예제는 1.2 GHz에서 25/100 Hz RF, 3개 peak, 105–135 ppm 전체 범위를 사용하며
R1H/R2H를 입력하지 않는다. 자세한 절차는 [영문](SBONEST_MANUAL.md) 또는
[한글 매뉴얼](SBONEST_MANUAL.ko.md)을 참고한다.

결과는 **실행한 작업 디렉터리**의 `<Project Name>_result.txt`, `_result.json`,
`.pdf`, `_data.pdf`에 저장한다. `--no-pdf`로 그림을 생략할 수 있다.
동일 이름의 결과가 있으면 덮어쓰지 않고 중단하므로 `Project Name`을 바꾼다.
`sbfit.py config.json`도 같은 fitting 경로를 사용한다.
입력 데이터·OC waveform의 상대 경로는 config 파일의 디렉터리를 기준으로 한다.

## 입력과 모델

ONEST 파일의 첫 세 줄은 **¹⁵N Larmor 주파수(MHz)**, saturation time(s),
`v1 v1err`(Hz)이고, 네 번째 줄은 컬럼 제목이다. 데이터 offset은 절대 ¹⁵N ppm이다.
각 잔기는 `# A1 R2a: 12 R2b: 15 dw: 3`과 같은 헤더로 시작한다.
`dw = δB − δA`이며, initial 값의 부호도 데이터에 맞춰 지정한다.

`init.Method`는 `Sideband`로 지정하고 다음 section을 추가한다.
아래 수치는 실행 예제의 **합성 실험 조건**이다. 실제 실험값으로 교체해야 한다.

```json
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
```

- 모든 active residue의 `h_ppm_a`, `h_ppm_b`가 필요하다. B-state의 ¹H shift가
  알려져 있지 않다면 A와 같다고 명시적으로 가정할 수 있지만 그 가정의 민감도를 확인해야 한다.
- `h_larmor_mhz`, `h_carrier_ppm`, `p90_s`는 필수 입력이다. `J_hz=92`, `b1_scale=1`은
  생략 시 사용하는 모델값이다. R1H/R2H를 모두 생략하면 초기값 (2,25) s⁻¹에서
  peak별로 독립 fitting한다. 같은 peak의 RF 파일과 A/B 상태는 같은 H rates를 쓴다.
  기존처럼 decoupling에 H rate를 명시하면 고정 모드를 유지한다.
  `sideband.proton_relaxation.mode`로 `fit`/`fixed`를 명시할 수도 있다.
  자동 모드는 단일 proton 자기장만 지원하며, `fit` 모드에서는 decoupling의
  R1H/R2H를 제거한다. 선택적인 시작값·범위는 `init`에 `<peak>.R1H/R2H`로 지정한다.
- `RR`은 같은 R=90x–240y–90x를 반복하는 stock CPD다. `R`, `RRbar`, `MLEV4`도
  지원한다. RRbar는 두 번째 R의 모든 위상을 180° 바꾼 별도 순서다.
- RF amplitude는 `b1_scale/(4*p90_s)`이고 pulse duration은 고정한다.
  따라서 B1 miscalibration을 바꿀 때 pulse length까지 바꾸는 오류를 피한다.
- carrier/decoupling이 데이터마다 다르면 `sideband.datasets`에 데이터 파일과
  같은 순서의 override object 목록을 준다. 예: `[{}, {"h_carrier_ppm": 8.4}]`.
  서로 다른 B₀의 자동 H fitting은 따로 실행한다. 고정 H 모드에서는 field와
  H rates를 dataset별로 명시해 다중 field를 함께 fitting할 수 있다.

Hamiltonian은 각 site에서
`H/2π = ΩN Nz + v1n Nx + ΩH Hz + v1Hx Hx + v1Hy Hy + J NzHz`다.
OC의 spin operator·column-major Liouvillian을 사용해 실제 NH 두 spin을 계산한다.
각 state에 16개의 real product operator, 두 state에 총 32차원을 사용하며,
모든 constant RF 구간에 RF·offset·J·relaxation·exchange를 동시에 포함한
`expm(L Δt)`를 적용한다. 완전 주기의 matrix power와 남은 구간을 정확히 전파한다.
작은 timestep으로 Hamiltonian을 나누는 근사나 sideband 위치 고정은 사용하지 않는다.

초기 상태는 `(pA Nz, pB Nz)`, 검출값은 `Nz(A,T)/pA`이다.
ONEST Matrix와 같은 longitudinal deviation convention으로 equilibrium recovery는 없다.
Relaxation은 product operator의 N/H rate 합을 사용한다. CSA–DD cross-correlation,
다른 spin과의 coupling, H chemical exchange, 시간에 따른 RF drift 등은 포함하지 않는다.
`J=0`에서 원래 ONEST Matrix 계산과 일치하는 것을 테스트한다.

90x–240y–90x 반복 decoupling과 잔여 sideband의 실험적 배경은
[Vallurupalli et al., JACS 2012, DOI 10.1021/ja3001419](https://doi.org/10.1021/ja3001419)의
Fig. 2–3을 참고할 수 있다. 이 문헌은 여기 구현의 검증을 대신하지 않는다.

## v1n fitting

| mode | fitted RF parameter | 설명 |
|---|---|---|
| `fixed` (기본값) | 없음 | 파일의 명목 v1 사용 |
| `scale` | `v1n_scale` | 모든 스펙트럼에서 `v1n = scale × nominal v1` |
| `per_dataset` | `v1n[0]`, `v1n[1]`, … | 파일별 actual RF (Hz), 파일 내 잔기들이 공유 |

`per_dataset`는 `initial: [25, 50]`, `bounds: [[20,30],[40,60]]`처럼 Hz 단위의
목록을 받는다. 생략하면 초기값은 nominal, 범위는 nominal의 0.5–1.5배다.
`scale`도 범위를 생략하면 0.5–1.5를 사용한다. RF는 양수로 제한한다.
`v1err`는 ONEST와 같은 Gaussian RF inhomogeneity sampling 폭이며,
RF를 fitting할 때 `v1err/v1` 비율을 유지한다. 이 값은 fitted v1n의 표준오차가 아니다.

함께 fitting하는 parameter는 전역 `kab`, `kba`, 잔기별 `peak_ppm`, `dw_ppm`,
`R1`, `R2a`, `R2b`이며 자동 H 모드에서는 peak별 `R1H`, `R2H`도 포함한다.
`kex=kab+kba`, `pB=kab/kex`로 보고한다.
여러 잔기·여러 B1·여러 T의 데이터를 같은 config로 fitting할 수 있다.
다중 field를 함께 사용하는 경우 위의 고정 H 모드 제한을 따른다.

초기값·범위·고정값은 다음처럼 지정한다. `vary` 생략 시 모든 parameter를 fitting한다.
`initial`의 값만 지정하면 고정되지 않는다. 고정하려면 `vary`에서 제외해야 한다.

```json
"init": {
  "Method": "Sideband",
  "kex": {"min": 100, "max": 500, "nsteps": 3},
  "pB": {"min": 0.02, "max": 0.08, "nsteps": 3},
  "initial": {"kab": 15, "kba": 285, "A1.peak_ppm": 120, "A1.R1": 1.5},
  "bounds": {"A1.dw_ppm": [1, 5]},
  "vary": ["v1n_scale"],
  "max_nfev": 500
}
```

위 예제의 `vary`는 **RF만** fitting한다. 교환·완화도 fitting할 때는 `vary`를
생략하거나 해당 parameter 이름을 모두 포함한다. kex/pB의 min/max는 초기 grid이며
최종 fit bounds가 아니다. grid는 시작점 선택용이고 global optimum을 보장하지 않는다.
나쁜 초기 dw나 좁은 offset 범위는 다른 최소점·비식별성으로 이어질 수 있다.

결과 JSON에는 RF Hz, 모든 parameter·표준오차, scaled Jacobian rank/condition,
RF와 다른 parameter의 correlation, bound 도달 여부, 사용한 config가 포함된다.
표준오차는 입력 intensity error를 절대 σ로 본 local covariance이며 reduced χ²로
재조정하지 않는다. rank가 부족하면 free parameter의 표준오차는 `null`이다.
고정된 parameter의 `stderr=0`은 데이터로 정밀하게 측정됐다는 뜻이 아니다.

**v1n을 fitting parameter로 사용할 수 있음은 합성 데이터에서 확인했다.**
실측에서 식별 가능한지는 RF 수·T·offset sampling·sideband SNR과 고정된 ¹H 조건에
달려 있다. 특히 ¹H B1/shift가 틀렸을 때 그 오차를 v1n이 흡수할 수 있다.
RF scale은 먼저 공통 parameter로 사용하고, 별도 calibration 근거가 있을 때
파일별 v1n을 고려한다. 다중 field에 공통 scale을 쓸 때도 동일 보정이 타당한지 확인한다.
현재 제공된 데이터는 합성 데이터이며 실측 fitting 결과는 없다.

## OC waveform 사용

기본 composite pulse 대신 OC `optimalcontrol.io.export_json`의 파일을 넣을 수 있다.
이 파일은 **실제로 반복할 전체 주기와 supercycle**을 포함해야 한다.

```json
"decoupling": {
  "h_larmor_mhz": 1200.0,
  "h_carrier_ppm": 8.5,
  "waveform_json": "my_oc_period.json",
  "rf_hz": 5000.0,
  "b1_scale": 1.0
}
```

Waveform은 x/y 두 channel, 0부터 시작하는 균등 시간축, `metadata.pulse_dt`가 필요하다.
`units="a.u."`이면 실제 amplitude로 변환할 `rf_hz`가 필수다.
`units="Hz"`이면 `rf_hz`를 생략한다. waveform에는 `p90_s`/`cycle`을 함께 주지 않는다.
서로 다른 waveform sample 수만큼 expm 비용이 증가한다. 정확히 같은 이웃 sample만 합친다.
`cestdec.Scheme` JSON과 OC Waveform JSON은 다른 형식이다.

## 검증 재실행

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg
.venv/bin/python test_sideband.py
.venv/bin/python demo_sideband.py --out results/demo_new
.venv/bin/python demo_sideband.py --offset-step 5 --out results/uniform_5Hz_new
.venv/bin/python plot_full_profile.py results/uniform_5Hz_new
.venv/bin/python test_performance.py
.venv/bin/python test_debugging.py
.venv/bin/python verify_3state.py
```

`test_sideband.py`는 H-first complex density-matrix 기준 계산과 N-first real OC 계산을
비교한다. pulse phase/부분주기, J=0 ONEST 환원, native OC waveform,
RF-only 두 fitting mode를 검증한다. `demo_sideband.py`는 독립 기준으로 합성 데이터를 만들고
실제 `run.py` CLI를 통해 모든 parameter를 fitting하며, 결과·log·입력·그림을 보존한다.

`--offset-step 5`는 −2400~+2400 Hz에 5 Hz 간격으로 **합성 입력 데이터 자체**를 생성해
다시 fitting한다(961점/RF). 기존 불균일 데이터를 보간해 fitting하는 옵션이 아니다.
간격은 유한한 양수여야 하고 전체 4800 Hz 범위를 정확히 나눠야 한다.
기존 결과는 유지되며 `--out`에 새 디렉터리를 지정한다.

전용 웹 UI와 sideband Monte Carlo 경로는 추가하지 않았다. 기존 ONEST 웹/MC는 원래
모델들에 사용한다. 새 Sideband config를 기존 MC에 전달하면 unsupported method로 실패한다.
GPL-3.0 ONEST 소스와 LICENSE를 유지했다. 원본 복사 기준 ONEST commit:
`6d178f3e5d6dc82b4cfc7a546d9004c17ff492d9`; OC source commit:
`139ff23aa53cb03902d725a818d847ad55895fc1`.
