# AlphaFold3 재현 및 변이 감별력 검정

**원논문** Yao, M. *et al.* (2024) *Deciphering proteins in Alzheimer's disease: A new Mendelian
randomization method integrated with AlphaFold3 for 3D structure prediction.*
**Cell Genomics 4, 100700.** [`10.1016/j.xgen.2024.100700`](https://doi.org/10.1016/j.xgen.2024.100700)

---

## 한 줄 요약

> 논문은 알츠하이머 인과 단백질의 미스센스 변이가 **3D 구조를 바꾼다**고 보고했다.
> 동일 조건으로 재현한 결과 **그 구조 변화는 모델 자체의 변동과 구분되지 않았고**,
> 결합 예측 모델(Boltz-2)로 확장해도 **변이를 감별하지 못했다.**
> → 구조 예측으로 변이 효과를 논하기 전에 **모델의 변이 감별력을 먼저 검정해야 한다.**

---

## 가설과 검정

| | 내용 |
|---|---|
| **귀무가설 H₀** | 야생형과 변이형의 예측 구조 차이는 **동일 서열 반복 예측의 변동과 같다** |
| **대립가설 H₁** | 변이형의 구조 차이가 그 변동을 **초과한다** |
| **음성 대조군** | **동일 서열을 반복 예측한 구조 쌍** (변이 없음 → 차이는 모델 변동뿐) |
| **처리군** | 야생형 × 변이형 구조 쌍 |
| **검정** | Mann–Whitney U (단측), α = 0.05 |
| **통제 변수** | 동일 도구·동일 서열 구간·시드 고정·동일 샘플 수 |

---

## 결과

### 1. 논문 재현 — 성공

| 단백질 | 변이 | 재현 pTM | 논문 보고 pTM |
|---|---|---|---|
| **CD33** | R69G | **0.59–0.60** | **0.6** ✅ |
| PILRA | R78G | 0.44–0.46 | 미보고 |
| TREM2 | R62H | 0.52–0.54 | 미보고 |

논문은 AlphaFold 방법을 STAR Methods에 기술하지 않았다. 입력 구간이 명시되지 않아
두 근거로 **전장 서열**임을 복원했다 — (i) 논문 Fig 4B·S5 원본에 무질서 꼬리가 존재,
(ii) 전장일 때만 pTM 0.6이 재현됨.

### 2. 구조 축 — 3개 중 1개는 노이즈와 구분 불가

변이 잔기 8 Å 이내 Cα 국소 RMSD, Ig 도메인 정렬 후 재정렬 없이 측정.

| 단백질 | 음성 대조군 (동일 서열, 20쌍) | 처리군 (WT×변이형, 25쌍) | p | 판정 |
|---|---|---|---|---|
| **CD33** | 0.124 Å | **0.120 Å** | **0.54** | **H₀ 기각 실패** |
| TREM2 | 0.168 Å | 0.257 Å | 0.0054 | H₀ 기각 (차이 0.09 Å) |
| PILRA | 0.100 Å | 0.995 Å | 6×10⁻⁹ | H₀ 기각 |

**CD33은 논문이 본문 Figure 4B에 대표로 실은 사례다.**

### 3. 결합 축 — Boltz-2도 감별 실패

TREM2 Ig 도메인(19–130) + Aβ42(APP 672–713) 복합체, 시드 3 × 모델 5 = 30 구조.

| 지표 | WT | R62H | p | 판정 |
|---|---|---|---|---|
| **ipTM** (결합면 신뢰도) | 0.851 | 0.824 | **0.30** | 구분 안 됨 |
| 인터페이스 pLDDT | 0.791 | 0.787 | 0.90 | 구분 안 됨 |
| pTM | 0.935 | 0.923 | 0.65 | 구분 안 됨 |

복합체 예측 자체는 성공했다(ipTM 0.85는 높은 값). **변이 감별력만 없다.**

---

## 선행 연구의 한계와 본 연구의 차별성

| | 선행 연구 | 본 연구 |
|---|---|---|
| 원논문 (Yao 2024) | WT 1개 vs 변이형 1개를 겹쳐 "구조가 변했다"고 기술 | **음성 대조군을 두고 통계 검정** |
| AlphaFold 한계 (Buel 2022) | AF2가 미스센스 변이에 둔감함을 지적 | **정량화**하여 단백질별 판정 |
| Boltz-2 한계 (Bret 2026) | affinity가 결합부위 변이에 둔감함을 보고 | **AD 표적에서 독립 확인** |
| 실험 문헌 (Zhao 2018 등) | R62H가 Aβ 결합을 감소시킴 | 두 예측 모델 모두 이를 **못 잡음**을 확인 |

**본 연구의 기여**: 구조·결합 두 축에서 동일한 검정 프로토콜을 적용해,
예측 모델의 **변이 감별력 자체를 측정하는 절차**를 제시했다.

---

## 재현 방법

### 환경

```bash
python3 -m pip install numpy scipy matplotlib pandas requests
conda create -n pymol -c conda-forge python=3.11 pymol-open-source -y   # RMSD·렌더링용
```

### 실행 순서

```bash
# 1. 입력 서열 생성 (UniProt 조회 + 잔기 번호 assert 검증)
python3 src/sequences.py

# 2. AlphaFold Server 업로드용 JSON 생성
python3 src/make_af_jobs.py
#    → alphafoldserver.com 에서 inputs/af_jobs/af_jobs_A_full_seed1.json 업로드
#    → 결과를 folds_*/ 로 내려받는다 (본 저장소는 data/alphafold3/ 에 정리본 포함)

# 3. 신뢰도 집계
python3 src/collect_confidence.py

# 4. RMSD 쌍 계산 (PyMOL 환경)
/opt/anaconda3/envs/pymol/bin/python src/rmsd_analysis.py

# 5. 논문 스타일 그림 + 슬라이드 표/차트
/opt/anaconda3/envs/pymol/bin/python src/pymol_render.py
python3 src/compose_figure.py
python3 src/make_slide_tables.py
python3 src/make_slide_charts.py

# 6. Boltz-2 확장 (Colab T4) — notebooks/02_boltz_trem2_ab42.ipynb
python3 src/make_boltz_inputs.py
python3 src/analyze_boltz.py
```

### 재현성 장치

- **잔기 번호 `assert` 강제 검증** — UniProt 서열 버전이 바뀌면 즉시 중단 (`src/sequences.py`)
- **시드 고정** — AF Server / Boltz-2 모두 시드를 명시하고 기록
- **매니페스트 CSV** — 모든 입력 서열의 구간·길이·변이 위치를 기록
- **이어달리기** — 이미 끝난 작업은 건너뛰므로 세션이 끊겨도 재개 가능

---

## 저장소 구조

```
src/
  sequences.py            UniProt 조회 · 변이 적용 · 잔기 번호 검증
  make_af_jobs.py         AlphaFold Server 업로드용 JSON 생성
  collect_confidence.py   pTM · pLDDT 집계 (mmCIF 파서 자체 구현)
  rmsd_analysis.py        모델 쌍별 RMSD 135개 계산 (PyMOL)
  pymol_render.py         논문 Fig 4B · S5 스타일 렌더링
  compose_figure.py       논문 레이아웃으로 합성
  make_slide_tables.py    슬라이드용 표 생성 + 통계 검정
  make_slide_charts.py    노이즈/신호 분포 차트
  make_boltz_inputs.py    Boltz-2 입력 YAML 생성
  analyze_boltz.py        ipTM 집계 + 검정

notebooks/
  02_boltz_trem2_ab42.ipynb   Colab T4용 Boltz-2 실행 노트북

inputs/     FASTA · AF Server JSON · Boltz YAML
data/       예측 구조 60개 + 신뢰도 JSON (MSA·PAE는 재생성 가능하므로 제외)
results/    집계 CSV · 통계 요약 · 슬라이드 표
figures/    논문 재현 그림 · 분석 차트
```

---

## 사용 도구

| 도구 | 용도 | 라이선스 |
|---|---|---|
| [AlphaFold Server (AF3)](https://alphafoldserver.com) | 구조 예측 (논문과 동일) | 비영리 한정 |
| [Boltz-2](https://github.com/jwohlwend/boltz) 2.2.1 | 복합체·결합 예측 | MIT |
| PyMOL 3.1.0 (open-source) | 정렬 · RMSD · 렌더링 | 오픈소스 |
| SciPy · pandas · matplotlib | 통계 · 시각화 | BSD |

> **라이선스 주의** — AlphaFold Server 약관은 개인·비영리 조직만 허용한다.
> 본 저장소는 학술 재현 목적이다. 상업적 확장은 MIT 라이선스인 Boltz-2 기반이어야 한다.

---

## 참고 문헌

1. Yao, M. *et al.* (2024) *Cell Genomics* **4**, 100700. — 재현 대상
2. Abramson, J. *et al.* (2024) *Nature* **630**, 493. — AlphaFold3
3. Passaro, S. *et al.* (2025) — Boltz-2
4. **Buel, G. & Walters, K.** (2022) *Nat Struct Mol Biol* — *Can AlphaFold2 predict the impact of missense mutations on structure?*
5. **Bret, G.** *et al.* (2026) *J Chem Inf Model* — Boltz-2의 결합부위 변이 둔감성
6. Zhao, Y. *et al.* (2018) *Neuron* **97**, 1023. — TREM2가 Aβ 올리고머에 결합, AD 변이가 결합 감소
7. Zhong, L. *et al.* (2018) *Mol Neurodegener* — oAβ1-42 결합, 필수 잔기 31–91
8. Yeh, F. *et al.* (2016) *Neuron* **91**, 328. — TREM2–APOE/CLU/LDL 결합
9. Pillai, J. *et al.* (2025) *Comput Struct Biotechnol J* — TREM2 R62H의 구조 영향이 R47H보다 작음
10. Jansen, I. *et al.* (2019) *Nat Genet* **51**, 404. — AD GWAS

데이터베이스: UniProt `P20138` `Q9UKJ1` `Q9NZC2` `P05067` · RCSB PDB `4NFB` `3WUZ` `5UD7` `5UD8`

---

## 한계

- 예측 대 예측 비교이며 **실험적 검증은 없다**
- "차이가 없다"는 **결합이 변하지 않는다는 뜻이 아니라 모델이 감별하지 못했다**는 뜻이다
- ipTM은 결합력(affinity)이 아니라 인터페이스 신뢰도다 — Boltz-2 affinity head는 저분자 전용
- MR-SPI 통계 파트는 재현하지 않고 논문 값을 인용했다 (UKB-PPP는 통제접근 자원)
