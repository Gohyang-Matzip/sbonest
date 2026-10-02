SBONEST sideband fitting은 [SBONEST 매뉴얼](SBONEST_MANUAL.ko.md)을 참고한다.
아래 매뉴얼은 원래 ONEST 모델을 설명한다.

# ONEST 단계별 사용 매뉴얼

ONEST는 CEST NMR 데이터를 다중 자기장 동시 피팅(Baldwin 해석해 모델, Matrix 및 3-state 변형 포함)으로 분석하여 단백질의 보이지 않는 들뜬 상태(excited state)를 규명하는 도구입니다.

## 1. 설치

```bash
git clone https://github.com/jhyeokchoi/ONEST/
cd ONEST
python3 -m venv venv
source venv/bin/activate   # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

의존성: `numpy`, `scipy`, `matplotlib` (CLI); `flask`, `werkzeug` (웹 서버 전용).

## 2. 입력 데이터 파일 준비

B0 자기장/실험당 텍스트 파일 하나. 형식 (`example/syn10.txt` 참고):

```
80.12           # B0 자기장 (MHz)
0.400000        # 포화 시간 T (s)
10.00  0.20     # B1(v1) 세기(Hz)와 오차
#offset(ppm)    Intensity    error
# A1 R2a: 10.0 R2b: 10.0 dw: -5.0
 100.000    0.672    0.006
 100.201    0.682    0.006
 ...
```

- **1행**: B0 자기장 (MHz). **2행**: 포화 시간 (초). **3행**: B1 세기(Hz)와 오차.
- **4행**: 컬럼 헤더 주석 (무시됨).
- **잔기(residue) 헤더**: 각 잔기 블록의 시작. 두 가지 형식:
  - 간단형: `# A1`
  - 전체형 (권장): `# A1 R2a: 10.0 R2b: 10.0 dw: -5.0` — 초기값 제공. 특히 초기 `dw`를 올바른 부호로 주면 부호가 반전된 국소 최소값에 빠지는 것을 막을 수 있습니다.
- **데이터 행**: `offset(ppm)  intensity  error`, 한 점당 한 줄.
- 다중 자기장: 자기장별로 파일을 나누고 `prepare.py`에 모두 전달. 잔기 이름은 파일 간에 일치해야 합니다.

원시 스펙트럼에서 입력 파일을 생성하려면 [input4onest](https://github.com/jhyeokchoi/input4onest)를 참고하세요.

## 3. 설정 파일 생성 (prepare.py)

```bash
python prepare.py example/syn10.txt > config.json
# 다중 자기장:
python prepare.py field1.txt field2.txt > config.json
# 계산 방법 선택 (기본: Baldwin):
python prepare.py --method Matrix example/syn10.txt > config.json
```

방법: `Baldwin` (2-state 해석해, 빠름 — 기본값), `Matrix` (2-state 수치해), `NoEx` (교환 없음), `Matrix_3st_Linear`, `Matrix_3st_Triangle` (3-state).

## 4. config.json 수정 (선택)

```json
{
    "Project Name": "default",
    "init": {
        "kex": {"min": 10.0, "max": 400.0, "nsteps": 6},
        "pB":  {"min": 0.01, "max": 0.1,  "nsteps": 6},
        "Method": "Baldwin"
    },
    "datasets": ["example/syn10.txt"],
    "residues": [{"name": "A1", "flag": "on"}, ...]
}
```

- `Project Name` — 모든 출력 파일의 접두어.
- `init.kex`, `init.pB` — 전역 교환 속도와 들뜬 상태 분율의 그리드 탐색 범위. kex > 400 s⁻¹이 예상되면 `max`를 늘리세요 (예: 2000).
- `residues[].flag` — `"off"`로 설정하면 해당 잔기를 피팅에서 제외.

## 5. 피팅 실행

```bash
python run.py config.json
```

Matrix 방법 사용 시에는 BLAS를 단일 스레드로 강제하여 멀티프로세싱 오버헤드를 피하세요:

```bash
export OMP_NUM_THREADS=1
python run.py config.json
```

출력 파일 (`Project Name` = `default`인 경우):

| 파일 | 내용 |
|------|------|
| `default_data.pdf` | 잔기별 원시 데이터 그래프 |
| `default_result.txt` | 피팅 파라미터, Chi2, dof, 잔기별 결과 |
| `default.pdf` | 데이터 + 피팅 곡선 |

`default_result.txt`에서 전역 `kex`, `pB`, 잔기별 `R2a`, `R2b`, `dw`, Chi2를 확인하세요. 각 잔기의 피팅된 `dw` 부호가 물리적으로 타당한지 검토하세요 — 부호가 반전되고 해당 잔기의 chi2가 유독 높으면 국소 최소값 신호입니다 (2단계의 전체형 잔기 헤더로 초기 `dw`를 제공하면 해결).

## 6. 몬테카를로 오차 추정 (mcrun.py)

```bash
python mcrun.py config.json 100          # 100회 MC, 모든 CPU 코어 사용
python mcrun.py config.json 100 4        # 4개 프로세스로 제한
```

각 실행은 오차 범위 내에서 데이터를 재추출하여 다시 피팅합니다. 추가 출력: `default_mc.txt` (파라미터 평균 ± 표준편차), `default_mcmean.pdf`.

참고 (macOS): 파이프된 stdin이 아닌 저장된 스크립트/파일로 실행하세요 — `spawn` 방식이 `<stdin>`을 재임포트하지 못합니다.

## 7. 웹 인터페이스 (선택)

```bash
python server_run.py
```

브라우저에서 `http://127.0.0.1:5001` 접속. 데이터 업로드, 설정, 피팅을 웹에서 수행할 수 있습니다.

## 8. 벤치마크 (선택)

```bash
export OMP_NUM_THREADS=1
python benchmark.py config.json           # 시간 측정
python benchmark.py config.json profile   # + cProfile → benchmark_profile.prof
```

## 문제 해결

- **`Field line parse error` / `V1 line parse error`** — 헤더 3행 확인: 숫자 1개, 숫자 1개, 숫자 2개 순서.
- **`Residue X not found in dataset`** — `config.json`의 잔기 이름이 데이터 파일 헤더와 불일치.
- **한 잔기가 Chi2를 지배** — `dw` 부호 반전 국소 최소값 가능성이 높음; 전체형 잔기 헤더로 초기 `dw`를 제공하세요.
- **kex가 그리드 경계에 도달** — `init.kex.max`를 넓히고 재실행.
