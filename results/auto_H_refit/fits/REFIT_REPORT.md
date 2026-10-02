# R1H/R2H 입력 없이 전체 sideband profile 재fitting

1.2 GHz·105–135 ppm·약 25 Hz 간격의 기존 3-peak **합성 데이터**를 그대로 사용했다.
90°x–240°y–90°x decoupling sideband를 포함한 모든 점을 fitting했다. 새로운 실측 결과는 아니다.

## 프로그램 동작

사용자는 R1H/R2H를 입력하지 않아도 된다. 프로그램이 (2,25) s⁻¹에서 시작해
각 peak의 두 값을 질소·교환·RF 파라미터와 함께 최적화한다. 임의의 값을 고정하는 방식이 아니다.
H rates는 서로 다른 peak 사이에 공유하지 않으며, 같은 peak의 RF 파일과 A/B 상태 사이에서는 공유한다.
24개 자유 파라미터의 covariance에 H-rate 불확도를 포함한다. 경계에서의 대칭 국소 오차는 정식 confidence interval이 아니다.

## 기본 자동 fitting 결과

| RF (Hz) | 관측점 | 자유 파라미터 | kex (s⁻¹) | pB | v1n_scale | Reduced χ² |
|---|---:|---:|---:|---:|---:|---:|
| 25/100 | 882 | 24 | 299.68025 | 0.0498302 | 1.0822687 | 0.93569 |
| 25/50/100 | 1323 | 24 | 299.49797 | 0.0499607 | 1.0806102 | 0.94625 |

생성값은 kex=300 s⁻¹, pB=0.05, v1n_scale=1.08이다.

| 파라미터 | 생성값 | 25/100 Hz 추정값 ± SE | 25/50/100 Hz 추정값 ± SE |
|---|---:|---:|---:|
| v1n_scale | 1.08 | 1.0822687 ± 0.002242 | 1.0806102 ± 0.0018 |
| A1.peak_ppm | 120 | 120.00003 ± 0.0003871 | 120.00001 ± 0.0003614 |
| A1.dw_ppm | 3 | 3.000912 ± 0.001108 | 3.0014423 ± 0.0009472 |
| A1.R1 | 1.5 | 1.5003824 ± 0.000374 | 1.5001365 ± 0.0003038 |
| A1.R2a | 12 | 11.910548 ± 0.06159 | 11.969175 ± 0.05126 |
| A1.R2b | 15 | 15.148129 ± 0.4466 | 15.08502 ± 0.4106 |
| G2.peak_ppm | 118.5 | 118.49988 ± 0.0003472 | 118.50002 ± 0.0003331 |
| G2.dw_ppm | 2.2 | 2.2011034 ± 0.001058 | 2.2005838 ± 0.0009092 |
| G2.R1 | 1.3 | 1.29965 ± 0.0003319 | 1.2999311 ± 0.0002695 |
| G2.R2a | 10 | 9.9686827 ± 0.0528 | 10.000368 ± 0.04382 |
| G2.R2b | 20 | 19.645365 ± 0.3931 | 19.562122 ± 0.3569 |
| S3.peak_ppm | 121.5 | 121.50033 ± 0.0006526 | 121.50006 ± 0.0005984 |
| S3.dw_ppm | -2.5 | -2.502344 ± 0.0013 | -2.5024185 ± 0.00112 |
| S3.R1 | 1.7 | 1.7001173 ± 0.0004021 | 1.7001703 ± 0.0003263 |
| S3.R2a | 16 | 15.971953 ± 0.08149 | 16.038106 ± 0.06704 |
| S3.R2b | 25 | 24.538194 ± 0.5 | 24.595263 ± 0.4551 |

![전체 profile 재fitting](full_profile_fit.png)

## 초기값 재검사

R1H는 0.2–20, R2H는 5–100 s⁻¹에서 peak별로 독립 log-uniform 난수를 뽑아
세 번 재시작했다(seed 20261005). 두 RF 조합에 같은 초기값 세트를 사용했다.

| 조합 | 기본 fit과의 최대 절댓값 Δχ² | kex 범위 (s⁻¹) | v1n_scale 범위 |
|---|---:|---:|---:|
| two_RF | 4.33685e-08 | 299.680244–299.680248 | 1.08226872–1.08226876 |
| three_RF | 1.95342e-07 | 299.497971–299.497971 | 1.08061012–1.08061019 |

무잡음 대조군 두 개에서 전체 24개 생성값 복원도 확인했다.

## 내부 H-rate 진단

입력이 필요 없다는 것이 H rates가 중요하지 않거나 정확히 결정된다는 뜻은 아니다.
개별 H rates가 약하게 결정되거나 하한에 도달하면 진단 파일에 기록한다.

| 조합 | Peak | R1H ± SE (s⁻¹) | R2H ± SE (s⁻¹) |
|---|---|---:|---:|
| two_RF | A1 | 8.6647e-07 ± 17.18 | 25.903 ± 8.422 |
| two_RF | G2 | 0.78352 ± 5.839 | 19.64 ± 2.705 |
| two_RF | S3 | 3.1063e-16 ± 7.857 | 34.271 ± 3.852 |
| three_RF | A1 | 3.0283e-09 ± 16.36 | 26.382 ± 8.053 |
| three_RF | G2 | 0.54216 ± 5.61 | 19.698 ± 2.637 |
| three_RF | S3 | 0.24509 ± 7.645 | 34.41 ± 3.734 |

two_RF 진단: Weakly constrained proton nuisance parameters: A1.R1H, G2.R1H, S3.R1H; Strong peakwise R1H/R2H correlation; individual H rates are poorly separated; Parameters at bounds: S3.R1H


three_RF 진단: Weakly constrained proton nuisance parameters: A1.R1H, G2.R1H, S3.R1H; Strong peakwise R1H/R2H correlation; individual H rates are poorly separated

## 실행 및 산출물

기본 설정에는 R1H/R2H 입력이 없다. 동일 설정을 재실행하려면 Project Name을 새 출력 경로로 바꾼다.

```bash
cd /Users/donghanlee/work/projects/sbonest
.venv/bin/python refit_automatic_proton.py --out results/auto_H_repeat_01
```

[두 RF 설정](two_RF.json) · [두 RF 결과](two_RF_result.json) · [세 RF 설정](three_RF.json) · [세 RF 결과](three_RF_result.json) ·
[전체 결과](summary.json) · [측정점별 예측·잔차](predictions.csv) · [그림 PDF](full_profile_fit.pdf) · [입력·코드 hash](metadata.json)
