# 1.2 GHz: uniform 25 Hz sampling and sideband fitting

25 Hz 간격의 독립 합성 데이터를 새로 생성해 fitting했다. 기존 5 Hz 결과를 보간하거나 재표본화한 결과가 아니다. 전체 profile은 `sideband_full_profile.png` 및 `.pdf`, 확대 profile은 `sideband_fit.png` 및 `.pdf`에 있다.

## 간격 선택과 조건

- 대표 간격: **25 Hz**, 1.2 GHz 장비의 ¹⁵N 주파수 121.5949416 MHz에서 **0.205600658 ppm**.
- 범위: site A 기준 −2400…+2400 Hz, RF 세기당 193점, 3개 RF에서 총 579점. 같은 범위의 5 Hz 간격(961점/RF)보다 약 5배 적다. 측정점당 scan 수와 acquisition 조건이 같다는 비교이며, 실제 측정 시간을 재지는 않았다.
- 근거: Vallurupalli, Bouvignies & Kay, JACS 2012, DOI [10.1021/ja3001419](https://doi.org/10.1021/ja3001419), Methods, pp. 8158–8159. 500·800 MHz 실험에서 25 Hz 간격을 사용했고, 일부 RF 조건에서는 15 또는 30 Hz를 사용했다. 따라서 25 Hz를 대표값으로 택했으며, 1.2 GHz의 보편적 표준이라고 주장하지 않는다. [저자 제공 원문](https://pound.med.utoronto.ca/lek-publications/348.pdf)
- ¹H decoupling: 반복 90x–240y–90x, RR, p90=70 μs, nominal RF=3.571 kHz. H carrier 8.5 ppm, H shifts A/B=6.2/6.5 ppm, JNH=92 Hz.
- Saturation T=0.4 s; nominal ν₁N=25/50/100 Hz, 실제값=27/54/108 Hz. 실제 RF scale=1.08.
- 참값: kex=300 s⁻¹, pB=5%, ΔδN=3 ppm, R1=1.5 s⁻¹, R2A/B=12/15 s⁻¹. 독립 Gaussian intensity noise σ=0.001, seed=20261002.
- 모든 비교는 동일한 noisy data 또는 그 부분집합을 사용했다. `masked`는 |offset|=1250…1850 Hz를 제외해 두 주요 sideband 영역을 제거했다. `core`는 |offset|<600 Hz만 사용했다. 이 두 경우도 동일한 NH sideband 모델을 사용하므로, 비교는 fitting 모델 변경이 아닌 포함 데이터의 효과이다.

## 실제 fitting 결과

| Fitting | 점 수 | kex / s⁻¹ | pB / % | ν₁N scale ± local SE | χ²/dof |
|---|---:|---:|---:|---:|---:|
| 전체, RF scale fitting | 579 | 300.5143 | 4.98626 | 1.080604 ± 0.003295 | 1.10965 |
| Sideband 영역 제외, RF scale fitting | 429 | 300.5695 | 5.02012 | 1.075961 ± 0.003933 | 1.07222 |
| 중심 영역만, RF scale fitting | 141 | 299.6851 | 5.03931 | 1.073689 ± 0.004204 | 1.13150 |
| 전체, RF를 nominal 값에 고정 | 579 | 302.7128 | 5.65260 | 1.0, fixed | 2.10770 |

Noisefree 579점에서는 모든 8개 parameter를 참값으로 복원했고 χ²/dof=1.22×10⁻¹⁵였다. RF를 세 데이터셋별로 별도 fitting하면 27.0750, 54.2129, 108.5540 Hz, χ²/dof=1.11050이었다. 별도 RF fitting이 이 공통 scale 오차 예제에서 개선을 보이지 않아, 간단한 공통 scale을 우선한다.

전체 및 masked fit 모두 Jacobian rank=8/8, bound에 걸린 parameter는 없다. 전체 fit의 ν₁N scale–R2A 상관계수는 **−0.95166**이고, 프로그램이 strong-correlation 경고를 기록했다. SE는 입력 intensity error를 절대 오차로 취급한 local Jacobian covariance이며, 반복 실험이나 Monte Carlo 신뢰구간이 아니다. 서로 다른 데이터 부분집합의 χ²/dof 차이를 모델 우열의 검정으로 사용하지 않았다.

## 1.2 GHz에서 유리한 점과 한계

**이미 전체 sweep에 포함된 sideband를 버리지 않고 fitting하는 것은 유리하다.** 이 예제에서 sideband 영역을 포함하면 ν₁N scale SE가 0.003933에서 0.003295로 **16.2% 감소**하고, R2A SE도 0.10272에서 0.08370 s⁻¹로 감소한다. `v1n`은 fitting parameter로 사용 가능하지만 R2A와의 강한 상관관계 때문에 독립 RF calibration을 대체한다고 해석해서는 안 된다.

**Sideband를 얻기 위해 측정 시간을 늘리는 이득은 작다.** 전체/마스킹 비교의 점 수는 579/429=1.35배다. 동일 총 시간에서 마스킹한 지점의 시간을 나머지 지점의 추가 scan으로 돌리고, 오차가 scan 수의 제곱근에 반비례한다고 가정하면 masked RF SE는 약 0.003386이다. 전체의 0.003295는 이보다 약 **2.7%** 작다. 중심 영역만 측정하는 141점의 경우 같은 가정의 SE는 약 0.002075이므로, 이 단일 residue의 RF 정밀도만 목표로 하면 넓은 sweep이 효율적이지 않다. 이는 overhead를 무시한 local covariance 추정이며, 여러 residue가 전체 offset sweep을 공유하는 실제 실험 설계의 우열을 정하지 않는다.

**같은 decoupling RF 세기에서 고자장은 sideband를 ppm 기준으로 더 가까이 가져온다.** 위 논문의 반복 inversion pulse 관계 ±1/(2·pwINV)에 p90=70 μs, pwINV=(420/90)·p90를 적용하면 nominal sideband 거리는 1530.6 Hz이다. 이는 600 MHz에서 ¹⁵N 25.18 ppm, 1.2 GHz에서 12.59 ppm에 해당한다. 이 계산은 on-resonance 기준이며 실제 off-resonance NH 시뮬레이션에서는 주요 sideband가 약 ±1.72 kHz 부근에 나타난다. 따라서 1.2 GHz에서는 같은 ppm 관측 범위 안으로 sideband가 들어올 가능성이 커져 모델링의 필요성이 높아질 수 있다. 동시에 같은 proton ppm offset이 더 큰 Hz offset이 되므로 실제 proton carrier, pulse timing, RF calibration을 정확히 반영해야 한다. Decoupling 세기까지 자장에 맞춰 늘리는 경우에는 이 단순 비교가 달라진다.

**실무 판단:** 25 Hz uniform sampling으로 전체 sweep을 이미 측정한다면 sideband와 공통 ν₁N scale을 함께 fitting하는 것이 합리적이다. 25 Hz 간격에서 sideband의 추가 정보는 제한적이고 좁은 dip을 성기게 표본화하므로, sideband만을 위해 측정 범위를 크게 늘리는 것은 이 결과로 정당화되지 않는다. ν₁N을 실제보다 8% 작게 고정한 이번 경우에는 pB가 5%에서 5.653%로 치우쳤다. RF 오차를 다루는 효과와 sideband를 추가하는 효과를 구분해야 한다.

이는 단일 NH, 두 상태, 알려진 proton 조건, RF inhomogeneity=0의 합성 검증이다. 실제 1.2 GHz 데이터 또는 600 MHz와의 감도·완화 비교는 아니다. 실제 sample/probe의 RF inhomogeneity와 proton 조건 불확실성이 포함되면 이득과 bias가 달라질 수 있다.

## 재현

새 출력 폴더를 사용한다. 기존 파일을 덮어쓰지 않는다.

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg .venv/bin/python demo_sideband.py --offset-step 25 --out example/sideband_uniform_25Hz_new
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 MPLBACKEND=Agg .venv/bin/python plot_full_profile.py example/sideband_uniform_25Hz_new
```

이 폴더의 `noisy_masked.json`과 `masked_*Hz.txt`에는 위에 명시한 mask를 적용한 입력을 보존했다. 해당 config의 Project Name을 새 이름으로 바꾼 뒤 이 폴더에서 `../../.venv/bin/python ../../run.py noisy_masked.json --no-pdf`로 재실행할 수 있다. 각 fit의 JSON·text·log와 독립 생성 조건 `validation.json`을 함께 보존했다.

실제 입력 파일의 193점 및 25 Hz 등간격을 확인했고, 전체 profile을 만들 때 저장된 fit과 재계산한 χ²의 일치를 확인했다. PNG를 열어 전체 범위와 sampling을 시각적으로 검사했다.
