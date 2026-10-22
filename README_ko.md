# fptp-seat-audit

[![DOI](https://zenodo.org/badge/DOI/__CONCEPT_DOI__.svg)](https://doi.org/__CONCEPT_DOI__)

**요인 기반 FPTP 의석 분포의 수치 감사: 2024 영국 총선 사례.**
계산 노트 [`docs/SEATDIST_NOTE_v0_2_3.md`](docs/SEATDIST_NOTE_v0_2_3.md)(초안 v0.2.3, 동료 심사 전)의 코드, 자료, 저장된 출력입니다. 노트는 영어로 썼고, 이 문서는 영문 README 의 한국어판입니다.

> **범위.** 여기 모형은 고정된 계산 사례입니다(모형 v0.2.1, 2024 최종 여론조사). 예측 상품이 아니며, 2024 년 노동당의 의석 효율을 크게 잘못 봤습니다. 다루는 것은 주어진 요인 모형이 함의하는 의석 분포를 어떻게 계산하고 어떻게 점검하는가입니다.

## 노트가 보이는 것

공통 요인(전국 여론조사 오차, 지역 스윙, 선거구 유형 요인)이 고정되면 선거구 결과가 서로 독립이라는 모형이 있다고 합시다. 그러면 요인 값마다 각 정당의 의석 분포를 Poisson-binomial 재귀(신용 포트폴리오 손실 모형에서 쓰는 것)와 지역 사이의 합성곱으로 구할 수 있습니다. 요인이 주어지면 이 계산은 대수적으로 정확합니다. 구현에는 부동소수 오차와 상한이 있는 절단 오차가 더해지고, 66차원 요인 적분이 주된 근사입니다.

가장 많이 쓴 방식은 혼합 방식입니다. 전국·선거구 유형 요인 9개에는 scrambled Sobol 점을, 지역 요인에는 Gauss–Hermite 규칙을 씁니다. 이 방식의 통상적 오차 막대(독립 scramble 사이의 퍼짐)는 성격이 다른 오차 두 개를 놓쳤습니다.

1. **지역 규칙 편향.** P(Green ≥ 4석)이 약 3e-3 틀렸는데 scramble 표준오차는 3.5e-5 였습니다. 바깥 오차 막대는 안쪽 규칙의 평균 편향을 보지 못합니다.
2. **덜 풀린 드문 바깥 기여.** scramble 4개로는 P(Other ≥ 5)의 추정치와 그 표준오차가 모두 불안정했습니다.

오차 성분은 네 가지 계산 경로로 나눴습니다.
- 고정 바깥 점에서의 쌍 지역 규칙 비교
- 지역 구적을 쓰지 않는 조건부 RQMC 기준값
- 지역 전용 RQMC 편향 추정
- 독립된 생성형 몬테카를로

강화한 지역 규칙에서는 18개 사건 모두, 혼합 방식과 기준값의 차이에 대한 사건별 명목 95% 구간이 [−1e-3, 1e-3] 안에 들었습니다(편향 보정을 한 보조 비교에서는 18개 중 16개). 이는 사건별 탐색적 비교이며 동시적·결정적 오차 보증이 아닙니다.

재귀의 새로움이나 몬테카를로 대비 일반적 속도 우위는 주장하지 않습니다. 기여는 흔한 오차 막대가 놓친 두 실패를 각각 다른 계산 경로로 찾아낸 재현 가능한 사례라는 점입니다.

## 설치

```
python -m pip install -r requirements.txt
```

검증한 버전은 `requirements.txt` 에 고정돼 있습니다(Linux x86-64, Python 3.13.16). 아래 1·2단계는 numpy 와 scipy 만 있으면 됩니다. 3단계, 예제, 시험에는 pandas 와 Numba 도 필요합니다.

Windows 에서는 UTF-8 모드로 실행하세요(`python -X utf8 ...` 또는 `set PYTHONUTF8=1`). 이 경우에도 생성된 Markdown 은 줄바꿈(CRLF)만, JSON 은 부동소수 끝자리만 저장본과 다를 수 있습니다. Windows 에서 한 독립 실행에서는 1e-14 수준의 차이가 났습니다. 표와 모든 판정은 같습니다.

## 재현 세 단계

**1. 저장된 출력에서 노트와 표를 다시 만들기 (몇 초).**

```
python scripts/seatdist_note_tables_v0_2.py      # out/seatdist_note_tables_v0_2.md (표 M1, N1-N9)
python scripts/seatdist_note_tables_v0_2_2.py    # out/seatdist_note_tables_v0_2_2.md (표 N7d)
python scripts/build_seatdist_note_v0_2_3.py     # docs/SEATDIST_NOTE_v0_2_3.md
```

`scripts/build_seatdist_note_v0_2_2.py` 와 `scripts/build_seatdist_note_v0_2_1.py` 는 이전 판 `docs/SEATDIST_NOTE_v0_2_2.md`, `docs/SEATDIST_NOTE_v0_2_1.md` 를 그대로 다시 만듭니다.

**2. 저장된 감사 결과를 다시 집계하기 (몇 초).** 최종 비교(표 N7a–c)를 저장된 실행에서 다시 계산합니다.

```
python scripts/r1_note_tables_v02.py             # out/r1_note_tables_v02.json, 이어서 1단계
```

**3. 실험 다시 돌리기 (2코어에서 몇 시간).** 스크립트와 순서는 노트 부록 F 에 표별로 있습니다. 세 가지를 주의하세요.
- 저장된 지역 편향 실행은 난수열 수정 이전 것입니다. 수정된 스크립트는 새 파일에 쓰고 저장본을 덮어쓰지 않습니다(부록 C).
- 저장된 조건부 실행은 dense 재귀 커널(`UKSEAT_PB=dense`)로 돌렸습니다. 지금 기본 커널과는 1.4e-14 까지 같습니다.
- `scripts/r1_followup_v02.py` 는 커널 성능 비교와 Green 축 감사를 함께 다시 돌리고 저장 파일을 덮어씁니다. 성능 비교만 다시 확인하려면 `scripts/pb_kernel_bench_v0_1.py` 를 쓰세요.

CI 워크플로(`.github/workflows/tests.yml`)는 단위 시험을 돌리고, 1·2단계가 커밋된 노트와 표를 그대로 재현하는지 확인합니다.

## 감사한 지역 규칙으로 실행하는 예제

```
python scripts/example_audited_v0_2.py 8 4       # 바깥 2^8 x 4 scrambled Sobol, 지역 규칙 (3,3,3,9,5,3)
```

엔진 생성자의 기본 지역 규칙 `(3,3,2,2,2,2)` 는 저장된 실행이 재현되도록 남겨 둔 것입니다. 위에서 말한 Green ≥ 4 편향(약 +3e-3)이 있습니다. 예제처럼 `ukseat/engine2.py` 의 `REG_NS_AUDITED` 를 명시해서 넘기세요.

`Hung (no party>=326)` 으로 저장된 사건은 "어떤 모형 범주도 326석에 이르지 못함"이라는 뜻입니다. Oth(무소속과 군소 정당을 합친 것)를 하나의 범주로 셉니다. 노트는 이 이름으로 보고합니다. 정부 구성 예측이 아닙니다.

같은 모형과 바깥 점에서의 선거구별 승리 확률은 `QEngine2(...).seat_marginals(...)` 로 구합니다. 인자는 `run(...)` 과 같습니다. 지역 규칙이 크면 `run()` 보다 훨씬 느립니다. (v0.2.1 에서는 이 메서드가 상속된 것이라 `QEngine2` 와 `FastEngine` 에서 틀린 값을 냈습니다. 노트 부록 F 참고.)

## 시험

```
python tests/test_engine_v0_1.py
python tests/test_engine2_v0_2.py
```

첫 파일은 작은 선거의 전수 열거와 재귀의 일치, 절단 상한, 구적 모멘트, 닫힌 형식 혼합과 Gumbel-max 표본의 일치를 확인합니다. 둘째 파일은 선거구 유형 항을 포함한 v0.2 핵심 계산을 전수 열거와 비교하고, `QEngine2` 와 `FastEngine` 의 선거구별 확률 합이 같은 인자로 돌린 `run()` 의 기대 의석과 같은지 확인합니다.

## 구성

| 경로 | 내용 |
|---|---|
| `docs/SEATDIST_NOTE_v0_2_3.md` | 노트(현재 판)와 템플릿. `v0_1`, `v0_2`, `v0_2_1`, `v0_2_2` 는 이전 판으로 보관 |
| `ukseat/` | 모형(`model.py`, `model2.py`), 엔진(`engine.py`, `engine2.py`, Numba 핵심 `core.py`), 구적, 2024 입력(`config_2024.py`), 사건(`report.py`) |
| `scripts/` | 자료 준비, 실험, 표와 노트 생성 |
| `data/` | 입력(출처와 라이선스는 `DATA_LICENSES.md`) |
| `out/` | 저장된 출력(JSON, 로그, 생성된 표) |
| `tests/` | 단위 시험 |

## 인용

모든 판을 인용할 때는 concept DOI [__CONCEPT_DOI__](https://doi.org/__CONCEPT_DOI__) 를 쓰세요. 이 DOI 는 늘 최신 판으로 연결됩니다. v0.2.2 는 [__V022_DOI__](https://doi.org/__V022_DOI__) 이고, 이 판(v0.2.3)의 DOI 는 Zenodo 기록에 있습니다. GitHub 의 "Cite this repository" 버튼은 [`CITATION.cff`](CITATION.cff) 를 씁니다.

> Han, N. (2026). *fptp-seat-audit: a numerical audit of factor-based FPTP seat distributions (Great Britain 2024 case study)*, version 0.2.3. Risk & Return Co., Ltd. Zenodo. https://doi.org/__CONCEPT_DOI__

## 라이선스

- **코드**(`ukseat/`, `scripts/`, `tests/`): MIT(`LICENSE`)
- **노트와 이 프로젝트가 만든 저장 출력**(`docs/`, `out/`): CC BY 4.0(귀속 표시 Nari Han)
- **자료:** 원 출처의 조건을 따릅니다. 영국 의회 선거 결과는 Open Parliament Licence v3.0(귀속 문구: Contains Parliamentary information licensed under the Open Parliament Licence v3.0), BES 2024 선거구 결과 파일은 CC BY 4.0 입니다. 자세한 내용과 귀속 표시는 [`DATA_LICENSES.md`](DATA_LICENSES.md) 에 있습니다.

## AI 사용 공개

생성형 AI 도구를 많이 썼습니다. Claude(Anthropic; Claude Opus 5.5, 2026년 10월)가 저자의 지시에 따라 코드 대부분을 쓰고, 수치 실험을 돌리고, 노트 초안을 썼습니다. GPT 6.1 Sol(OpenAI)은 계산과 초안을 여러 차례 검토했고 감사의 일부를 독립적으로 다시 실행했습니다. 그에 따른 수정은 감사 기록의 일부입니다. 저자는 질문을 정하고, 작업을 지시·검토하고, 보고할 결과와 표현을 결정했으며, 내용에 대한 모든 책임을 집니다.

## 연락

의견과 정정은 GitHub Issues 로 보내 주세요.
