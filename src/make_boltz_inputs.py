"""Boltz-2 입력 YAML을 만든다 — TREM2 ± Aβ42 복합체.

목적
  논문은 "변이 → 구조 변화"만 봤다. 우리 재현 결과 구조 변화는 노이즈 수준이었다.
  그런데 실험 문헌은 TREM2 변이가 **Aβ 결합**을 떨어뜨린다고 보고한다.
  → 구조가 아니라 결합을 봐야 한다. AlphaFold3는 결합을 못 재고, Boltz-2는 잰다.

근거
  · Zhao et al. 2018, Neuron  — TREM2가 Aβ 올리고머에 나노몰 친화도로 결합.
                                AD 관련 변이가 결합을 감소시킴
  · Zhong et al. 2018, Mol Neurodegener — oAβ1-42가 TREM2에 고친화도 결합.
                                결합에 필수인 잔기가 31–91 구간 (R62가 여기 포함)
  · Yeh et al. 2016, Neuron   — TREM2가 APOE·CLU·LDL 결합, 질병 변이가 결합 저해

⚠️ 한계 (반드시 함께 보고할 것)
  · Boltz-2의 affinity head는 **저분자 전용**이다. Aβ42는 펩타이드이므로
    결합력 수치(affinity)는 나오지 않는다. 우리가 얻는 건 **ipTM(결합 신뢰도)** 다.
    → `properties: affinity` 블록을 넣지 않는다.
  · Bret et al. 2026, JCIM — Boltz-2 affinity가 결합부위 변이에 둔감하다는 보고.
    → ipTM 차이도 노이즈 바닥 검정을 거쳐야 한다. 그냥 믿으면 안 된다.

구간 선택
  TREM2는 **Ig 도메인(19–130)** 만 쓴다. 전장을 쓰면 막관통·세포질 꼬리가
  복합체 인터페이스 지표(ipTM)를 오염시킨다. 19–130은 결정구조 5UD8·5ELI가
  덮는 구간이고, 결합 필수 잔기 31–91을 모두 포함한다.
  (논문 재현 파트는 전장을 썼다 — 그건 논문 조건을 맞추기 위함이고,
   이 확장 파트는 우리 설계이므로 목적에 맞는 구간을 쓴다.)

실행: python3 src/make_boltz_inputs.py
출력: inputs/boltz/*.yaml, inputs/boltz/README.md
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from sequences import TARGETS, fetch_uniprot, apply_mutation  # noqa: E402

OUT = ROOT / "inputs" / "boltz"

# Aβ42 = APP(P05067) 672–713. UniProt에서 직접 잘라 검증한다.
APP_ACC = "P05067"
AB42_RANGE = (672, 713)
AB42_EXPECTED = "DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA"

SEEDS = [1, 2, 3]          # 노이즈 바닥용. 실행 여유 보고 늘린다
DIFFUSION_SAMPLES = 5      # AF3의 job당 모델 5개와 맞춤


def yaml_for(trem2_seq: str, ab42: str, with_peptide: bool) -> str:
    """Boltz-2 입력 YAML. affinity 블록은 넣지 않는다 (펩타이드라 대상 아님)."""
    lines = [
        "version: 1",
        "sequences:",
        "  - protein:",
        "      id: A",
        f'      sequence: "{trem2_seq}"',
    ]
    if with_peptide:
        lines += [
            "  - protein:",
            "      id: B",
            f'      sequence: "{ab42}"',
        ]
    return "\n".join(lines) + "\n"


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)

    # --- Aβ42 확보 및 검증 -------------------------------------------------
    app = fetch_uniprot(APP_ACC)
    a, b = AB42_RANGE
    ab42 = app[a - 1 : b]
    assert len(app) == 770, f"APP 길이 {len(app)} — 770을 기대"
    assert ab42 == AB42_EXPECTED, "Aβ42 서열 불일치 — UniProt 버전 확인 필요"
    print(f"Aβ42  APP {a}-{b}  {len(ab42)} aa  ✓ 검증 통과")

    # --- TREM2 도메인 WT / R62H -------------------------------------------
    t = next(x for x in TARGETS if x["gene"] == "TREM2")
    s, e = t["construct"]                       # 19–130
    full_wt = fetch_uniprot(t["acc"])
    full_mut = apply_mutation(full_wt, t["pos"], t["wt"], t["mut"], t["gene"])
    dom_wt, dom_mut = full_wt[s - 1 : e], full_mut[s - 1 : e]
    label = f"{t['wt']}{t['pos']}{t['mut']}"

    assert len(dom_wt) == len(dom_mut) == e - s + 1
    assert sum(x != y for x, y in zip(dom_wt, dom_mut)) == 1, "변이가 1개가 아님"
    print(f"TREM2 {s}-{e}  {len(dom_wt)} aa  WT / {label}  ✓ 검증 통과")
    print(f"      결합 필수 구간 31-91 포함 여부: {s <= 31 and e >= 91}")

    # --- YAML 생성 ---------------------------------------------------------
    jobs = []
    for allele, seq in (("WT", dom_wt), (label, dom_mut)):
        for with_pep, tag in ((True, "AB42"), (False, "alone")):
            name = f"TREM2_{allele}_{tag}"
            (OUT / f"{name}.yaml").write_text(yaml_for(seq, ab42, with_pep))
            jobs.append((name, len(seq) + (len(ab42) if with_pep else 0)))

    for n, ln in jobs:
        print(f"  {n:28s} {ln:>4d} aa")

    # --- 실행 안내 ---------------------------------------------------------
    cmds = "\n".join(
        f"boltz predict inputs/boltz/{n}.yaml \\\n"
        f"    --out_dir outputs/boltz --use_msa_server \\\n"
        f"    --diffusion_samples {DIFFUSION_SAMPLES} --seed {seed} \\\n"
        f"    --output_format mmcif"
        for n, _ in jobs if n.endswith("AB42")
        for seed in SEEDS
    )

    (OUT / "README.md").write_text(f"""# Boltz-2 실행 — TREM2 ± Aβ42

## 무엇을 왜 하는가

논문은 변이의 **구조** 변화를 봤다. 우리 재현 결과 그 변화는 노이즈 수준이었다.
실험 문헌은 TREM2 변이가 **Aβ 결합**을 떨어뜨린다고 보고한다.
→ 결합을 재야 한다. AlphaFold3는 못 재고 Boltz-2는 잰다.

## 입력

| 파일 | 내용 | 길이 |
|---|---|---|
| `TREM2_WT_AB42.yaml` | TREM2 야생형 + Aβ42 | {len(dom_wt) + len(ab42)} aa |
| `TREM2_{label}_AB42.yaml` | TREM2 {label} + Aβ42 | {len(dom_mut) + len(ab42)} aa |
| `TREM2_WT_alone.yaml` | TREM2 야생형 단독 (대조) | {len(dom_wt)} aa |
| `TREM2_{label}_alone.yaml` | TREM2 {label} 단독 (대조) | {len(dom_mut)} aa |

TREM2는 Ig 도메인 {s}–{e} (결정구조 5UD8·5ELI 구간, 결합 필수 잔기 31–91 포함).
Aβ42는 APP P05067 {a}–{b}.

## 읽을 지표

- **ipTM** (0~1) — 두 사슬 사이 인터페이스 신뢰도. **이게 결합 지표다**
- **interface pLDDT** — 접촉면 잔기들의 확신도
- ❌ affinity 수치는 안 나온다. Boltz-2 affinity head는 **저분자 전용**이고
  Aβ42는 펩타이드다. `properties: affinity` 블록을 일부러 넣지 않았다.

## 실행 (Colab T4 권장)

```bash
!pip install -q "boltz[cuda]==2.2.1"
# 런타임 재시작 후

{cmds}
```

## ⚠️ 해석 시 반드시 함께 보고할 것

Bret et al. 2026 (J Chem Inf Model)은 Boltz-2 affinity가 **결합부위 변이에 둔감**하다고
보고했다. 따라서 WT와 {label}의 ipTM 차이도 **노이즈 바닥 검정을 거쳐야** 한다.
시드 {SEEDS} × diffusion sample {DIFFUSION_SAMPLES}개로 분포를 만들어,
구조 파트와 동일한 Mann–Whitney U 검정을 적용한다.

## 근거 문헌

- Zhao et al. 2018, *Neuron* — TREM2가 Aβ 올리고머에 나노몰 결합, AD 변이가 결합 감소
- Zhong et al. 2018, *Mol Neurodegener* — oAβ1-42 고친화도 결합, 필수 잔기 31–91
- Yeh et al. 2016, *Neuron* — TREM2–APOE/CLU/LDL 결합, 질병 변이가 저해
- Passaro et al. 2025, *bioRxiv* — Boltz-2 (MIT 라이선스, FEP 근접 성능)
- Bret et al. 2026, *JCIM* — Boltz-2의 결합부위 변이 둔감성 (한계)
""")

    print(f"\n→ {OUT.relative_to(ROOT)}/  YAML {len(jobs)}개 + README.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
