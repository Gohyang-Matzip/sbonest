# 균일 offset 격자 재피팅

−2400~+2400 Hz, **5 Hz 간격**, RF당 961점 / 전체 2883점.
기존 adaptive sampling 데이터와 결과를 보존하고, 독립 complex density-matrix 기준으로
이 격자에서 합성 데이터를 새로 생성해 실제 `run.py`로 fitting했다.
기존 데이터의 보간이나 그림의 x축 간격만 바꾼 결과가 아니다.

조건: ¹H 1200 MHz, p90=70 μs, RR 90x–240y–90x, T=0.4 s,
nominal v1n=25/50/100 Hz, 참 scale=1.08, kex=300 s⁻¹, pB=5%,
intensity σ=0.001, seed=20261002. ¹H 조건과 J 등은 이전 예제와 동일하다.

| 조건 | 데이터 수 | kex (s⁻¹) | pB (%) | RF scale 또는 파일별 Hz | reduced χ² |
|---|---:|---:|---:|---|---:|
| 무잡음, v1n scale fit | 2883 | 300.0000 | 5.00000 | 1.080000 ± 0.000870 | 0.00000 |
| 잡음, v1n scale fit | 2883 | 299.2572 | 5.01239 | 1.079405 ± 0.000870 | 1.01393 |
| 잡음, v1n 고정 | 2883 | 300.5573 | 5.67124 | 1.0 (고정) | 3.80660 |
| 잡음, 중앙 영역만 | 717 | 299.3535 | 5.00063 | 1.081070 ± 0.001912 | 1.09705 |
| 잡음, 파일별 v1n fit | 2883 | 299.0624 | 5.01206 | 26.9882, 53.9801, 107.8998 Hz | 1.01424 |

전체 범위 fitting은 모든 점에 같은 σ=0.001을 적용한다. 중앙 영역 조건만 |offset|<600 Hz로 제한한다.
JSON의 rank, covariance, correlation과 경고도 확인했다. 무잡음의 모든 parameter 복원,
잡음의 RF scale 복원 및 reduced χ² 검증을 통과했다.

[전체 profile](sideband_full_profile.png) · [전체 PDF](sideband_full_profile.pdf) ·
[main/sideband 확대 비교](sideband_fit.png) · [전체 결과 JSON](validation.json).

이전 adaptive 결과는 RF당 133점이었고 이번은 961점이므로 표준오차 차이를
균일 spacing만의 효과로 해석할 수 없다. 같은 seed이지만 표본 수·위치가 달라 잡음 실현도 다르다.
본 결과는 합성 데이터 검증이며 실측 fitting이 아니다.

재현(저장 경로는 새 디렉터리):

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg \
  .venv/bin/python demo_sideband.py --offset-step 5 --out results/uniform_new
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg \
  .venv/bin/python plot_full_profile.py results/uniform_new
```
