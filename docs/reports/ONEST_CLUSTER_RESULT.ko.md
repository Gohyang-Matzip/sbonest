# ONEST 클러스터 분석 검증 결과

**날짜:** 2026-08-08 · **도구:** ONEST (`python3.11 run.py config.json`) · **모델:** Baldwin 2-state
**Global 파라미터:** kab, kba (primary: kex = kab + kba) · **판정 절차:** clustering skill (red_chi2 게이트 → chi2_glob 랭킹 → kex_indiv gap → 그룹 fit + AICc)

## 1. 예제 데이터 (example/syn10.txt, A1–A5)

**결론: K = 1** — 5개 residue 전부 하나의 global 세트 공유. kab = 14.99 s⁻¹, kba = 298.36 s⁻¹ (kex ≈ 313.3 s⁻¹, pB ≈ 4.8%).

| 가설 | chi2 | k | AICc (n=1000) |
|---|---|---|---|
| **1-global** | 1023.93 | 27 | **1079.48** ← 최저 |
| K=2 {A1 \| A2–A5} | 1019.97 | 29 | 1079.76 |
| 전부 individual | 1017.04 | 35 | 1089.65 |

- A1의 red_chi2 1.28은 individual fit에서도 1.31로 남음(Δchi2 = 2.06) → kex 이질성이 아니라 residue-local 미스핏(오차 과소평가 가능성)
- A1 제외 refit: red_chi2 = 0.985, kex 313.3 → 314.9 (사실상 불변) → A1이 global 파라미터를 끌지 않음

## 2. 합성 데이터 검증 (ground truth K = 1, 2, 3, 4)

**설계:** 시나리오당 8 residue (B1–B8) × 200점 = n = 1600. 공통 pB = 0.05, R1 = 1, R2a = 10, R2b = 20 s⁻¹, B0 = 80.12 MHz, T = 0.4 s, v1 = 25 Hz, 노이즈 σ = 0.01. ONEST의 Baldwin 모델로 생성(시드 고정), full residue 포맷으로 초기 dw 제공.

### 결과: 4/4 시나리오 모두 정답 K 복원

| 시나리오 | 진짜 K (구성) | 복원 K | AICc: 1-global | AICc: K-global | AICc: individual |
|---|---|---|---|---|---|
| S1 | 1 (kex 300×8) | **1** ✓ | **1680.9** | — (gap 없음) | 1695.4 |
| S2 | 2 (150×4, 600×4) | **2** ✓ | 2998.5 | **1697.8** | 1709.5 |
| S3 | 3 (100×3, 350×3, 1000×2) | **3** ✓ | 3675.0 | **1591.9** | 1606.3 |
| S4 | 4 (100/300/700/1400 ×2) | **4** ✓ | 4124.1 | **1676.8** | 1689.2 |

AICc = chi2 + 2k + 2k(k+1)/(n−k−1), 모든 가설에서 n = 1600 (총 데이터 수) 동일.

### 클러스터별 kex 복원 (그룹 fit, 전부 red_chi2 ≈ 1로 내부 통과)

| 시나리오 | 복원 kex (진값) |
|---|---|
| S1 | 302.2 (300) |
| S2 | 146.3 (150) · 603.2 (600) |
| S3 | 107.6 (100) · 345.9 (350) · 1008.8 (1000) |
| S4 | 97.9 (100) · 299.0 (300) · 724.8 (700) · 1303.2 (1400) |

### 실증된 원리

1. **Global least-squares는 robust하지 않다.** K≥2 데이터에서 1-global fit의 kex는 항상 한 그룹으로 끌려감 (S2: 169, S3: 84, S4: 102) — 개별 `better` 플래그 대신 크기(chi2_glob 랭킹, red_chi2)를 봐야 하는 이유.
2. **chi2_glob 랭킹은 오염에도 살아남는다.** 오염된 global fit에서도 residue별 chi2 상승 경계가 진짜 클러스터 경계와 일치 (예: S2에서 B1–B4 ~200–233 vs B5–B8 384–703).
3. **kex_indiv gap이 곧 클러스터 가설.** 모든 시나리오에서 gap 수 = 진짜 그룹 수 − 1. S1은 gap이 없어 분할 가설 자체가 생성되지 않음 (과분할 방지).
4. **AICc는 진짜 K에서 최저.** K를 더 늘리면(individual = K=8 극한) 다시 상승.
5. **한계:** fast exchange (kex = 1400)에서는 복원값이 ~7% 낮게 나옴(1303) — 단, 클러스터 판정에는 영향 없음.

### 부수 발견: 초기 dw의 중요성

simple residue 포맷(초기 dw = 0.1)으로 fit하면 dw 부호가 반전된 국소 최소값에 빠질 수 있음. v1 검증에서 진짜 K=1 데이터인데 B5 하나가 dw +1.48 (진값 −2.5)로 수렴해 chi2 5351을 만들었고 red_chi2가 4.35까지 상승 — 클러스터 문제로 오판할 수 있는 패턴. **실데이터 분석 시 full residue 포맷(`# B5 R2a: 10 R2b: 20 dw: -2.5`)으로 초기 dw를 제공하고, 클러스터 판정 전 dw_fit 부호를 확인할 것.**

## 실행 정보

- 총 50 fit: global 4(+실패 v1 4), individual 32, group 9, 로컬 검증 1 — 병렬 computation agent로 실행
- 산출물: `session_artifacts/clustering/` (이 폴더에 보존 — indiv/group/global fit별 config, result.txt, PDF 전부 포함)
- 데이터 생성 스크립트: `scripts/gen_cluster_synth.py` (저장소 루트에서 실행, 시드 1–4로 재현 가능 — `python3.11 scripts/gen_cluster_synth.py` 실행 시 현재 폴더에 k1v2–k4v2.txt 생성)
