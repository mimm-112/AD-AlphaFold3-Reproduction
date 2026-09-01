"""Boltz-2 결과 집계 — TREM2 ± Aβ42 결합 신뢰도 비교.

읽는 지표
  iptm            두 사슬(TREM2·Aβ42) 사이 인터페이스 신뢰도. **이게 결합 지표다**
  complex_iplddt  인터페이스 잔기들의 pLDDT
  ptm             복합체 전체
  complex_plddt   복합체 전체 pLDDT

검정 설계
  ipTM은 구조당 스칼라 하나이므로 RMSD처럼 쌍으로 만들 필요가 없다.
  같은 서열 15개(시드 3 × 모델 5)의 흩어짐이 곧 노이즈 바닥이고,
  WT 15개 vs R62H 15개를 Mann–Whitney U로 직접 비교한다.

⚠️ affinity 수치는 없다. Boltz-2 affinity head는 저분자 전용이고 Aβ42는 펩타이드다.
⚠️ Bret et al. 2026(JCIM)은 Boltz-2가 결합부위 변이에 둔감하다고 보고했다.
   "구분 안 됨"이 나와도 그 자체가 보고할 결과다.

실행: python3 src/analyze_boltz.py
출력: results/boltz_iptm.csv, results/boltz_summary.md
"""

from __future__ import annotations

import csv
import json
import re
import statistics as st
import sys
from pathlib import Path

from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parent.parent
SEARCH = ROOT / "notebooks" / "boltz_results" / "outputs"
RESULTS = ROOT / "results"

METRICS = ["iptm", "complex_iplddt", "ptm", "complex_plddt", "confidence_score"]


def collect() -> list[dict]:
    rows = []
    for f in sorted(SEARCH.rglob("confidence_*.json")):
        d = json.loads(f.read_text())
        # 경로 예: outputs/TREM2_WT_AB42_seed1/boltz_results_.../predictions/.../confidence_..._model_0.json
        m = re.search(r"TREM2_(WT|R62H)_AB42_seed(\d+)", str(f))
        if not m:
            continue
        allele, seed = m.group(1), int(m.group(2))
        model = int(f.stem.rsplit("_", 1)[-1])
        row = {"allele": allele, "seed": seed, "model": model}
        row.update({k: d.get(k) for k in METRICS})
        rows.append(row)
    return rows


def describe(vals: list[float]) -> str:
    return (f"{st.median(vals):.3f}  "
            f"(min {min(vals):.3f} / max {max(vals):.3f}, n={len(vals)})")


def main() -> int:
    if not SEARCH.exists():
        print(f"결과 폴더 없음: {SEARCH}")
        return 1

    rows = collect()
    if not rows:
        print("confidence json 을 못 찾았다.")
        return 1

    RESULTS.mkdir(exist_ok=True)
    out_csv = RESULTS / "boltz_iptm.csv"
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    wt = [r for r in rows if r["allele"] == "WT"]
    mu = [r for r in rows if r["allele"] == "R62H"]

    lines: list[str] = []
    A = lines.append

    A("# Boltz-2 결과 · TREM2 ± Aβ42 결합 신뢰도\n")
    A(f"구조 {len(rows)}개 (WT {len(wt)} + R62H {len(mu)}) · "
      f"시드 {sorted({r['seed'] for r in rows})} × 모델 5개")
    A("입력: TREM2 Ig 도메인 19–130 (112 aa) + Aβ42 (42 aa) = 154 aa\n")

    A("## 지표별 비교\n")
    A("| 지표 | WT | R62H | 중앙값 차 | Mann–Whitney p | 판정 |")
    A("|---|---|---|---|---|---|")

    verdicts = {}
    for k in METRICS:
        a = [r[k] for r in wt if r[k] is not None]
        b = [r[k] for r in mu if r[k] is not None]
        if not a or not b:
            continue
        # 양측: 어느 쪽으로든 다르면 잡는다 (방향은 중앙값 차로 본다)
        _, p = mannwhitneyu(a, b, alternative="two-sided")
        diff = st.median(b) - st.median(a)
        sep = p < 0.05
        verdicts[k] = (st.median(a), st.median(b), diff, p, sep)
        A(f"| **{k}** | {describe(a)} | {describe(b)} | {diff:+.3f} | {p:.3g} | "
          f"{'다름' if sep else '**구분 안 됨**'} |")

    A("\n## 해석\n")
    mi = verdicts.get("iptm")
    if mi:
        wt_m, mu_m, diff, p, sep = mi
        A(f"- **ipTM** — 결합 지표. WT {wt_m:.3f} vs R62H {mu_m:.3f} (차이 {diff:+.3f}), p = {p:.3g}")
        if sep:
            direction = "낮아졌다" if diff < 0 else "높아졌다"
            A(f"  → 변이형에서 결합 신뢰도가 통계적으로 유의하게 **{direction}**.")
            if diff < 0:
                A("  → 실험 보고(Zhao 2018 · Zhong 2018: AD 변이가 Aβ 결합 감소)와 **방향 일치**.")
            else:
                A("  → ⚠️ 실험 보고와 **방향이 반대**다. 해석에 주의.")
        else:
            A("  → **Boltz-2가 이 변이를 구분하지 못했다.**")
            A("  → Bret et al. 2026(JCIM)이 보고한 '결합부위 변이 둔감성'과 일치하는 결과.")

    A("\n## 반드시 함께 보고할 한계\n")
    A("1. **affinity 수치가 아니다.** Boltz-2 affinity head는 저분자 전용이고 "
      "Aβ42는 펩타이드다. ipTM은 결합력이 아니라 **인터페이스 신뢰도**다.")
    A("2. **Bret et al. 2026** (*J Chem Inf Model*) — Boltz-2 affinity가 결합부위 변이에 둔감. "
      "타겟을 바꿔도 분류가 잘 안 바뀌는 사례 보고.")
    A("3. **King et al. 2025** — Boltz-2를 단백질–단백질 친화도로 미세조정해도 "
      "서열 기반 모델보다 성능이 낮았다.")
    A("4. 실험적 검증 없음. 예측 대 예측 비교다.")

    A("\n## 근거 문헌\n")
    A("- Zhao et al. 2018, *Neuron* — TREM2가 Aβ 올리고머에 나노몰 결합, AD 변이가 결합 감소")
    A("- Zhong et al. 2018, *Mol Neurodegener* — oAβ1-42 고친화도 결합, 결합 필수 잔기 31–91 (R62 포함)")
    A("- Yeh et al. 2016, *Neuron* — TREM2–APOE/CLU/LDL 결합, 질병 변이가 저해")
    A("- Passaro et al. 2025 — Boltz-2 (MIT, FEP 근접 성능)")
    A("- Bret et al. 2026, *JCIM* — Boltz-2 결합부위 변이 둔감성")

    (RESULTS / "boltz_summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    print(f"\n저장 → {out_csv.relative_to(ROOT)}, results/boltz_summary.md")
    return 0


if __name__ == "__main__":
    sys.exit(main())
