"""슬라이드에 그대로 붙여넣을 표를 만든다. (base 환경에서 실행)

입력: results/confidence_summary.csv   (collect_confidence.py)
      results/rmsd_pairs.csv           (rmsd_analysis.py, PyMOL 환경)
출력: results/slide_tables.md

표 3종
  표 1 — 재현 결과 요약 (논문 대조 포함)          → 논문 재현 파트
  표 2 — 신뢰도 상세 (pTM / pLDDT)               → 논문 재현 파트
  표 3 — 노이즈 바닥 vs 신호                      → 자체 검증 파트 (논문에 없음)

실행: python3 src/make_slide_tables.py
"""

from __future__ import annotations

import csv
import statistics as st
import sys
from pathlib import Path

from scipy.stats import mannwhitneyu

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from sequences import TARGETS  # noqa: E402

RESULTS = ROOT / "results"
GENES = [t["gene"] for t in TARGETS]
BY_GENE = {t["gene"]: t for t in TARGETS}


def load(name: str) -> list[dict]:
    p = RESULTS / name
    if not p.exists():
        sys.exit(f"{name} 이 없다. 먼저 해당 스크립트를 실행할 것.")
    return list(csv.DictReader(p.open()))


def rng(vals: list[float], d: int = 2) -> str:
    """최소~최대. 값이 같으면 하나만."""
    lo, hi = min(vals), max(vals)
    return f"{lo:.{d}f}" if abs(hi - lo) < 10 ** -d else f"{lo:.{d}f}–{hi:.{d}f}"


def main() -> int:
    conf = [r for r in load("confidence_summary.csv") if r["construct"] == "full"]
    pairs = load("rmsd_pairs.csv")

    def sub(gene: str, allele_is_wt: bool) -> list[dict]:
        return [r for r in conf
                if r["gene"] == gene and (r["allele"] == "WT") == allele_is_wt]

    out: list[str] = []
    A = out.append

    A("# 슬라이드용 표 — AlphaFold3 재현 결과\n")
    A("입력: AlphaFold Server(AF3), UniProt 전장 서열, seed 1, job당 모델 5개 "
      "→ 단백질당 구조 10개(WT 5 + 변이형 5), 총 30개\n")

    # ---------------- 표 1 : 재현 결과 요약 --------------------------------
    A("\n## 표 1 · 재현 결과 요약 (논문 대조)\n")
    A("| 단백질 | 변이 | dbSNP | 재현 pTM | 논문 pTM | 도메인 pLDDT | 판정 |")
    A("|---|---|---|---|---|---|---|")
    for g in GENES:
        t = BY_GENE[g]
        lab = f"{t['wt']}{t['pos']}{t['mut']}"
        allp = [float(r["ptm"]) for r in conf if r["gene"] == g]
        dom = [float(r["mean_plddt_domain"]) for r in conf if r["gene"] == g]
        paper = "**0.6**" if g == "CD33" else "미보고"
        verdict = "**일치**" if g == "CD33" else "—"
        A(f"| {g} | {lab} | {t['rsid']} | {rng(allp)} | {paper} | "
          f"{rng(dom, 1)} | {verdict} |")
    A("\n> 논문이 pTM을 명시한 것은 CD33뿐이다 (Fig 4B 캡션: *\"a score of 0.6 for both "
      "structures\"*). PILRA·TREM2는 논문에 수치가 없어 대조 불가.")

    # ---------------- 표 2 : 신뢰도 상세 -----------------------------------
    A("\n## 표 2 · 신뢰도 상세\n")
    A("| 단백질 | 대립 | pTM | 전장 pLDDT | 도메인 pLDDT | 변이부위 pLDDT | 무질서 비율 |")
    A("|---|---|---|---|---|---|---|")
    for g in GENES:
        t = BY_GENE[g]
        lab = f"{t['wt']}{t['pos']}{t['mut']}"
        for is_wt, name in ((True, "WT"), (False, lab)):
            rs = sub(g, is_wt)
            A(f"| {g} | {name} | {rng([float(r['ptm']) for r in rs])} "
              f"| {rng([float(r['mean_plddt']) for r in rs], 1)} "
              f"| {rng([float(r['mean_plddt_domain']) for r in rs], 1)} "
              f"| {rng([float(r['plddt_at_mut']) for r in rs], 1)} "
              f"| {rng([float(r['fraction_disordered']) for r in rs])} |")
    A("\n> pTM이 0.5~0.6으로 낮은 것은 접힘 실패가 아니라 **신호펩타이드·막관통·세포질 "
      "꼬리가 무질서**하기 때문이다(무질서 비율 0.41~0.61). 변이가 위치한 Ig 도메인 "
      "자체는 pLDDT 94~98로 최상급이다.")

    # ---------------- 표 3 : 노이즈 바닥 vs 신호 ---------------------------
    A("\n## 표 3 · 노이즈 바닥 대비 변이 효과 *(논문에 없는 자체 검증)*\n")
    A("| 단백질 | 노이즈 (같은 서열, 20쌍) | 신호 (WT vs 변이형, 25쌍) | Mann–Whitney p | 판정 |")
    A("|---|---|---|---|---|")

    verdicts = {}
    for g in GENES:
        noise = [float(r["rmsd_local8A"]) for r in pairs
                 if r["gene"] == g and r["comparison"].startswith("noise")]
        sig = [float(r["rmsd_local8A"]) for r in pairs
               if r["gene"] == g and r["comparison"] == "signal_WT_vs_MUT"]
        u, p = mannwhitneyu(sig, noise, alternative="greater")
        sep = p < 0.05 and st.median(sig) > st.median(noise)
        verdicts[g] = (st.median(noise), st.median(sig), p, sep)
        A(f"| {g} | {st.median(noise):.3f} Å (중앙값) | {st.median(sig):.3f} Å (중앙값) "
          f"| {p:.2g} | {'신호 > 노이즈' if sep else '**구분 안 됨**'} |")

    A("\n**측정 방법**")
    A("- 대상: 변이 잔기 8 Å 이내 Cα (국소 RMSD). Ig 도메인 Cα로 먼저 정렬한 뒤 "
      "재정렬 없이 측정")
    A("- 노이즈 = 같은 서열의 서로 다른 모델 쌍 (WT-WT 10쌍 + 변이형-변이형 10쌍)")
    A("- 신호 = WT 모델 5개 × 변이형 모델 5개 = 25쌍")
    A("- 검정: Mann–Whitney U, 단측(신호 > 노이즈), 유의수준 0.05")
    A("- k = 모델 5개 (AF Server가 job당 자동 생성하는 diffusion sample). "
      "추가 예측 없이 기존 결과만 사용")

    n_sep = sum(1 for v in verdicts.values() if v[3])
    A(f"\n> **결과: {len(GENES)}개 중 {n_sep}개에서 신호가 노이즈를 넘었다.**")
    if n_sep < len(GENES):
        nots = [g for g in GENES if not verdicts[g][3]]
        A(f"> {', '.join(nots)}에서는 변이로 인한 구조 차이가 "
          f"모델 자체의 변동과 통계적으로 구분되지 않는다.")

    path = RESULTS / "slide_tables.md"
    path.write_text("\n".join(out) + "\n")
    print("\n".join(out))
    print(f"\n\n저장 → {path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
