# SBONEST 검증 결과 — 2026-10-02 및 2026-10-03

실제 `run.py` CLI를 사용한 독립 합성 검증이다. 원본 입력·config·결과 JSON·text·log는
[example/sideband_demo](../../example/sideband_demo)에, 비교 그림은
[PNG](../../example/sideband_demo/sideband_fit.png) / [PDF](../../example/sideband_demo/sideband_fit.pdf)에 있다.
아래 첫 세 절은 2026-10-02의 1.2 GHz 검증 기록이다. 마지막 절은
2026-10-03의 600/800 MHz 모델 비교를 추가한다. 모두 합성 데이터이며
실측 데이터는 검증하지 않았다.

## 실험과 참값

- ¹H 1200 MHz, ¹⁵N 121.5949416 MHz, T=0.4 s, nominal v1n=25/50/100 Hz.
- p90=70 µs, stock RR 90x–240y–90x, H carrier=8.5 ppm, H shifts A/B=6.2/6.5 ppm.
- J=92 Hz, R1H=2 s⁻¹, R2H=25 s⁻¹ (고정한 합성 모델 조건).
- kex=300 s⁻¹, pB=0.05, δA=120 ppm, ΔδN=3 ppm, R1=1.5, R2A=12, R2B=15 s⁻¹.
- 실제 v1n은 nominal의 **1.08배**: 27, 54, 108 Hz.
- Gaussian intensity σ=0.001, seed=20261002. 무잡음 입력에도 error=0.001을 지정했다.
- 기준 데이터는 H-first complex density matrix를 직접 구성한 `test_sideband.reference`에서 생성했다.
  production은 OC를 사용하는 N-first real product-operator 계산이다.
- main/minor 영역과 양쪽 sideband를 샘플링했다. 이 H-offset 조건의 주요 sideband는 약 ±1720 Hz로,
  on-resonance 근사 위치에 강제로 두지 않는다.

## 동시 fitting

scale 모드는 kab/kba, RF scale, δA/ΔδN/R1/R2A/R2B의 총 8개 parameter를 동시에 풀었다.
파일별 RF 모드는 10개 parameter다. ¹H 조건은 모든 fitting에서 참값으로 고정했다.

| 조건 | 데이터 수 | reduced χ² | kex (s⁻¹) | pB | v1n scale ± local σ | 시간 (s) |
|---|---:|---:|---:|---:|---|---:|
| 무잡음 · 전체 · scale | 399 | 2.27249e-15 | 300.00000 | 0.0500000 | 1.080000 ± 0.001077 | 3.26 |
| 잡음 · 전체 · scale | 399 | 1.02424 | 300.21630 | 0.0502497 | 1.078586 ± 0.001078 | 4.00 |
| 잡음 · 전체 · RF 고정 | 399 | 14.0269 | 294.77560 | 0.0567454 | 고정 1.0 | 3.33 |
| 잡음 · |offset|<600 Hz · scale | 153 | 1.05172 | 297.69937 | 0.0502435 | 1.080855 ± 0.003817 | 1.93 |
| 잡음 · 전체 · 파일별 RF | 399 | 1.02883 | 299.95668 | 0.0502280 | 파일별 값은 아래 | 3.57 |

파일별 RF fitting 결과(참값 27/54/108 Hz):

- `v1n[0]` = 26.973698 ± 0.037451 Hz
- `v1n[1]` = 53.964977 ± 0.089710 Hz
- `v1n[2]` = 107.834565 ± 0.124686 Hz

**판정: v1n을 fitting parameter로 사용할 수 있다.** 무잡음에서 모든 참값을 복원했고,
잡음 조건에서도 RF를 포함한 모델의 reduced χ²가 약 1이었다. 잘못된 nominal RF에 고정하면
sideband 깊이를 재현하지 못하고 reduced χ²가 약 14가 되었다.
이 예제에서 전체 데이터의 RF scale 표준오차는 main/minor 영역만 쓸 때보다
**3.54배 작았다**. 추가 데이터 수와 baseline 정보도 달라지므로 보편적인 sideband-only
정보 증가율로 해석해서는 안 된다.

RF와 나머지 변수의 correlation, rank, bounds 경고는 각 `_result.json`에 있다.
현재 예제는 full rank이고 bound 도달이 없다. 이 결과는 실제 실험의 정확도를 보장하지 않는다.
특히 ¹H B1/shift/relaxation과 J의 모델 오차는 이 covariance에 포함되지 않는다.

## 수치·프로그램 검사

- `tests/test_sideband.py`: 독립 complex 계산과 RR/RRbar/MLEV4, 0/짧은 T/정수 주기/부분 주기 비교,
  absolute 및 relative tolerance 2×10⁻¹¹ 통과.
- J=0에서 ONEST Matrix로 환원: tolerance 2×10⁻¹⁰ 통과.
- 독립 기준 데이터에서 scale 및 파일별 RF-only fitting 복원 통과.
- OC Waveform JSON의 x/y 순서, a.u.→Hz 단위, fractional propagation, 잘못된 RF 설정 거부 통과.
- 별도 코드 검토에서 다중 잔기·중간 inactive residue·grouped Jacobian·13개 parameter 복원 확인.
  peak position 미분 간격 피드백을 반영해 grouped relative step을 10⁻⁶으로 설정했다.
- 원래 ONEST `tests/test_performance.py`, `tests/test_debugging.py`, `tests/verify_3state.py` 모두 통과.
  CLI, 기존 웹 입력/출력, Monte Carlo와 3-state regression도 포함한다.
- Ruff F 검사, 신규 코드 E/F/I 검사, py_compile, git diff --check 통과.
- 기본 CLI PDF를 렌더링해 확인했고, 복사된 ONEST의 외부 범례가 잘리는 문제를
  `bbox_inches="tight"`로 수정한 뒤 재렌더링했다. 원본 ONEST/OC는 그대로 보존했다.

검증 환경: Python 3.14.7, NumPy 2.5.3, SciPy 1.18.1,
optimalcontrol 0.5.0. 단일 BLAS thread. 최초 OC 환경에는 Flask가 없어
웹 회귀검증 두 항목이 실패했으나, sbonest 전용 `.venv`에 원래 requirements를 모두 설치한 뒤
전체 검사를 다시 실행하여 통과했다. 이 절은 최초 검증 당시의 기록이다. 이후 추가된 자동 peakwise H fitting은
[전체 profile 재fitting 보고서](../../results/auto_H_refit/fits/REFIT_REPORT.md)를 참고한다.
현재 배포 상태와 검증 절차는 [HANDOFF.md](../../HANDOFF.md)에 기록한다.

재실행 명령과 물리 모델의 적용 범위는 [SIDEBAND.md](../guides/SIDEBAND.md)에 있다.

## 600/800 MHz 모델 비교 — 2026-10-03

[배포 결과](../../results/field_comparison_600_800_20261003_02/REPORT.txt)는 자기장별로
SBONEST, 실제 ONEST Matrix, signed nitrogen-only 대조를 같은 관측점과 오차로
피팅한 기록이다. 105–135 ppm, 25/50/100 Hz RF, σ = 0.001, R1H = 2와
R2H = 25 s⁻¹ 고정, RF scale을 포함한 8개 자유 파라미터를 사용했다.
관측점은 600 MHz에서 225개, 800 MHz에서 297개이며 두 자기장을 따로 피팅했다.

| ¹H 자기장 | 모델 | kex ± SE (s⁻¹) | pB ± SE (%) | Reduced χ² |
|---|---|---:|---:|---:|
| 600 MHz | SBONEST | 297.89 ± 1.63 | 5.0287 ± 0.0299 | 1.058 |
| 600 MHz | ONEST Matrix | 297.95 ± 1.64 | 5.0142 ± 0.0357 | 1.055 |
| 800 MHz | SBONEST | 299.85 ± 1.23 | 4.9994 ± 0.0313 | 1.035 |
| 800 MHz | ONEST Matrix | 301.41 ± 1.23 | 4.9702 ± 0.0330 | 1.058 |

참값은 kex = 300 s⁻¹, pB = 5%다. ±는 reduced χ²로 재조정하지 않은 국소
1 SE이며, pB의 SE는 퍼센트포인트다. ONEST의 모델 불일치 오차는 포함하지 않는다.
600 MHz 결과는 거의 같고, 800 MHz 무잡음 ONEST에는 kex의 +0.51% 변위가 남았다.
SBONEST는 두 자장에서 무잡음 참값을 복원했다. 기존 1250–1850 Hz 마스크는
두 측정 창에서 점을 제외하지 않는다. 이 결과는 자기장별 잡음 표본 하나에 한정된다.

확인한 항목:

- 모델 3개 × 잡음/무잡음 2조건 × 자기장 2개 × 초기값 3개 = 36회 피팅 재현.
- 독립 density-matrix 생성기와 보존된 무잡음 입력 일치.
- J = 0에서 signed nitrogen-only와 NH 전파 일치 및 ONEST zero-clamp 관계 확인.
- 기존 SBONEST noisy fit 재현, 무잡음 참값 복원, Jacobian rank 8/8 확인.
- 저장된 목적함수·예측값·공분산과 모든 소스 해시 확인, 두 비교 그림 렌더링 확인.
- 잡음 포함 ONEST의 초기값별 χ² 차이 보존: 600 MHz 0.153006, 800 MHz 0.065897.
  최저 χ²를 채택했으며, 모든 초기값의 수렴점 일치나 전역 최적성을 주장하지 않는다.

검증 기록은 [verification.json](../../results/field_comparison_600_800_20261003_02/verification.json),
상세 결과와 해시는 [summary.json](../../results/field_comparison_600_800_20261003_02/summary.json)에 있다.
재현 명령은 [매뉴얼 7.1절](../manual/SBONEST_MANUAL.ko.md)을 따른다. 이 36회 연구 계산은
CI의 기본 회귀 검사와 별도로 실행했다. 코드와 결과를 배포한
[PR #3](https://github.com/Gohyang-Matzip/sbonest/pull/3)의 병합 커밋 `8d484de`는
[전체 CI](https://github.com/Gohyang-Matzip/sbonest/actions/runs/37114856838)를 통과했다.
