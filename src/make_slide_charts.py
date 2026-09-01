"""슬라이드 11·13용 시각자료를 만든다.

  figures/slide11_targets.png  — 7개 인과 단백질 중 구조 예측 가능한 3개 (표 형식)
  figures/slide13_noise.png    — 노이즈 바닥 vs 변이 신호 분포 (작은 배수 3패널)

색은 검증된 카테고리 팔레트 슬롯 1·2를 쓴다 (validate_palette.js 전 항목 PASS).
  노이즈 = 파랑 #2a78d6 (배경 역할) · 신호 = 주황 #eb6834 (주목 대상)

실행: python3 src/make_slide_charts.py
"""

from __future__ import annotations

import csv
import statistics as st
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402
from matplotlib.patches import Rectangle            # noqa: E402
import numpy as np                                  # noqa: E402
from scipy.stats import mannwhitneyu                # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from sequences import TARGETS                       # noqa: E402

FIG_DIR = ROOT / "figures"
RESULTS = ROOT / "results"

# --- 색 역할 (검증된 팔레트) ------------------------------------------------
SURFACE = "#ffffff"        # 슬라이드 배경에 맞춤
INK = "#0b0b0b"
INK_2 = "#52514e"
INK_MUTED = "#8a8a86"
GRID = "#e6e5e1"
C_NOISE = "#2a78d6"        # 슬롯 1 파랑
C_SIGNAL = "#eb6834"       # 슬롯 2 주황
C_NAVY = "#1e3a5f"         # 발표 타이틀색

plt.rcParams.update({
    "font.family": "AppleGothic",
    "axes.unicode_minus": False,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
})


# ===========================================================================
# 슬라이드 11 — 구조 예측 대상 표
# ===========================================================================
# 논문 Fig 4A(방향) + Fig 4B·S5(구조 유무) 기준
SEVEN = [
    # (단백질, 방향, 미스센스 변이 있음?, 변이 표기, 논문 그림)
    ("CD33",  "위험 ↑", True,  "R69G",  "Fig 4B"),
    ("PILRA", "위험 ↑", True,  "R78G",  "Fig S5(a)"),
    ("TREM2", "위험 ↓", True,  "R62H",  "Fig S5(b)"),
    ("PILRB", "위험 ↑", False, "—",     "없음"),
    ("RET",   "위험 ↑", False, "—",     "없음"),
    ("CD55",  "위험 ↓", False, "—",     "없음"),
    ("EPHA1", "위험 ↓", False, "—",     "없음"),
]


def slide11() -> Path:
    fig, ax = plt.subplots(figsize=(10, 5.4), dpi=200)
    ax.axis("off")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 9.4)

    ax.text(0.15, 9.0, "MR-SPI가 지목한 인과 단백질 7개 중, 구조를 예측할 수 있는 것은 3개",
            fontsize=13.5, color=C_NAVY, fontweight="bold", va="center")
    ax.text(0.15, 8.45,
            "나머지 4개는 미스센스 변이 자체가 없어 변이형 서열을 만들 수 없다",
            fontsize=10, color=INK_2, va="center")

    cols = [0.35, 2.35, 4.35, 6.6, 8.4]
    heads = ["단백질", "AD 위험 방향", "미스센스 변이", "변이 표기", "논문 그림"]

    y0, dy = 7.5, 0.86
    for x, h in zip(cols, heads):
        ax.text(x, y0, h, fontsize=10, color=INK_2, fontweight="bold", va="center")
    ax.plot([0.15, 9.85], [y0 - 0.34] * 2, color=INK_2, lw=1.2)

    for i, (gene, direction, has_mis, mut, figref) in enumerate(SEVEN):
        y = y0 - dy * (i + 1)
        if has_mis:
            ax.add_patch(Rectangle((0.15, y - 0.35), 9.7, 0.7,
                                   facecolor="#fdf0e9", edgecolor="none", zorder=0))
        ink = INK if has_mis else INK_MUTED
        weight = "bold" if has_mis else "normal"

        ax.text(cols[0], y, gene, fontsize=11.5, color=ink, fontweight=weight, va="center")
        ax.text(cols[1], y, direction, fontsize=10.5, color=ink, va="center")
        ax.text(cols[2], y, "○" if has_mis else "×", fontsize=13,
                color=C_SIGNAL if has_mis else INK_MUTED,
                fontweight="bold", va="center")
        ax.text(cols[3], y, mut, fontsize=10.5, color=ink, fontweight=weight, va="center")
        ax.text(cols[4], y, figref, fontsize=10, color=ink, va="center")

        if i == 2:  # 3개와 4개 사이 구분선
            ax.plot([0.15, 9.85], [y - dy / 2] * 2, color=GRID, lw=1.4)

    ax.text(0.15, 0.55,
            "초록은 \"seven proteins with structural alterations\"로 기술되어 있어 "
            "7개 모두 구조가 보고된 것처럼 읽힌다.",
            fontsize=9, color=INK_2, va="center")
    ax.text(0.15, 0.15,
            "출처: Yao et al. 2024, Figure 4A · 4B · S5",
            fontsize=8, color=INK_MUTED, va="center")

    out = FIG_DIR / "slide11_targets.png"
    fig.savefig(out, dpi=200, bbox_inches="tight", pad_inches=0.25)
    plt.close(fig)
    return out


# ===========================================================================
# 슬라이드 13 — 노이즈 바닥 vs 신호
# ===========================================================================
def slide13() -> Path:
    rows = list((RESULTS / "rmsd_pairs.csv").open() and
                csv.DictReader((RESULTS / "rmsd_pairs.csv").open()))
    genes = [t["gene"] for t in TARGETS]

    data = {}
    for g in genes:
        noise = [float(r["rmsd_local8A"]) for r in rows
                 if r["gene"] == g and r["comparison"].startswith("noise")]
        sig = [float(r["rmsd_local8A"]) for r in rows
               if r["gene"] == g and r["comparison"] == "signal_WT_vs_MUT"]
        _, p = mannwhitneyu(sig, noise, alternative="greater")
        data[g] = (noise, sig, p)

    xmax = max(max(v[0] + v[1]) for v in data.values()) * 1.18

    fig, axes = plt.subplots(len(genes), 1, figsize=(10, 6.4), dpi=200, sharex=True)
    fig.subplots_adjust(left=0.15, right=0.83, top=0.80, bottom=0.13, hspace=0.42)

    rng = np.random.default_rng(0)  # 지터 재현성 고정

    for ax, g in zip(axes, genes):
        noise, sig, p = data[g]
        sep = p < 0.05

        for vals, y, color, lab in ((noise, 1.0, C_NOISE, "노이즈"),
                                    (sig, 0.0, C_SIGNAL, "신호")):
            jit = rng.uniform(-0.17, 0.17, len(vals))
            ax.scatter(vals, np.full(len(vals), y) + jit, s=34, color=color,
                       alpha=0.55, linewidths=1.2, edgecolors=SURFACE, zorder=3)
            m = st.median(vals)
            ax.plot([m, m], [y - 0.32, y + 0.32], color=color, lw=2.6, zorder=4)
            ax.text(m, y + 0.42, f"{m:.3f}", fontsize=8.5, color=color,
                    ha="center", va="bottom", fontweight="bold")

        ax.set_ylim(-0.62, 1.72)
        ax.set_xlim(0, xmax)
        ax.set_yticks([1.0, 0.0])
        ax.set_yticklabels(["노이즈\n(같은 서열)", "신호\n(WT vs 변이형)"],
                           fontsize=9, color=INK_2)
        ax.tick_params(axis="y", length=0)
        ax.tick_params(axis="x", colors=INK_2, labelsize=9)
        ax.grid(axis="x", color=GRID, lw=0.8, zorder=0)
        ax.set_axisbelow(True)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.spines["bottom"].set_color(GRID)

        mut = next(f"{t['wt']}{t['pos']}{t['mut']}" for t in TARGETS if t["gene"] == g)
        ax.text(-0.005, 1.98, f"{g}  {mut}", transform=ax.get_yaxis_transform(),
                fontsize=11.5, color=INK, fontweight="bold", va="top", ha="right")

        verdict = "분리됨" if sep else "구분 안 됨"
        vcolor = C_SIGNAL if sep else "#b02020"
        ax.text(1.02, 0.72, verdict, transform=ax.transAxes, fontsize=10,
                color=vcolor, fontweight="bold", va="center")
        ax.text(1.02, 0.30, f"p = {p:.2g}", transform=ax.transAxes, fontsize=9,
                color=INK_2, va="center")

    axes[-1].set_xlabel("변이 잔기 8 Å 이내 Cα 국소 RMSD (Å)",
                        fontsize=10, color=INK_2, labelpad=8)

    fig.text(0.15, 0.955,
             "변이로 생긴 구조 차이가 모델 자체의 변동을 넘는가",
             fontsize=13.5, color=C_NAVY, fontweight="bold", ha="left")
    fig.text(0.15, 0.905,
             "노이즈 = 같은 서열의 모델 쌍 20개 · 신호 = WT × 변이형 모델 쌍 25개 · "
             "세로선 = 중앙값 · Mann–Whitney U 단측",
             fontsize=9, color=INK_2, ha="left")
    fig.text(0.15, 0.032,
             "CD33은 논문이 본문 Figure 4B에 대표로 실은 사례다. "
             "AlphaFold Server(AF3), 전장 서열, seed 1, job당 모델 5개.",
             fontsize=8.5, color=INK_MUTED, ha="left")

    out = FIG_DIR / "slide13_noise.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


# ===========================================================================
# 슬라이드 14 — Boltz-2 결합 축 (ipTM)
# ===========================================================================
def slide14() -> Path:
    import csv as _csv
    rows = list(_csv.DictReader((RESULTS / "boltz_iptm.csv").open()))
    wt  = [float(r["iptm"]) for r in rows if r["allele"] == "WT"]
    mut = [float(r["iptm"]) for r in rows if r["allele"] == "R62H"]
    _, p = mannwhitneyu(mut, wt, alternative="two-sided")

    fig, ax = plt.subplots(figsize=(10, 4.6), dpi=200)
    fig.subplots_adjust(left=0.17, right=0.80, top=0.68, bottom=0.22)
    rng = np.random.default_rng(0)

    for vals, y, color, lab in ((wt, 1.0, C_NOISE, "WT"),
                                (mut, 0.0, C_SIGNAL, "R62H")):
        jit = rng.uniform(-0.17, 0.17, len(vals))
        ax.scatter(vals, np.full(len(vals), y) + jit, s=42, color=color,
                   alpha=0.6, linewidths=1.2, edgecolors=SURFACE, zorder=3)
        m = st.median(vals)
        ax.plot([m, m], [y - 0.30, y + 0.30], color=color, lw=2.8, zorder=4)
        ax.text(m, y + 0.40, f"{m:.3f}", fontsize=9.5, color=color,
                ha="center", va="bottom", fontweight="bold")

    ax.set_ylim(-0.62, 1.66)
    ax.set_yticks([1.0, 0.0])
    ax.set_yticklabels(["TREM2 WT\n+ Aβ42", "TREM2 R62H\n+ Aβ42"],
                       fontsize=10, color=INK_2)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", colors=INK_2, labelsize=9)
    ax.grid(axis="x", color=GRID, lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.set_xlabel("ipTM — 두 사슬 인터페이스 신뢰도 (0~1, 높을수록 결합 확신)",
                  fontsize=10, color=INK_2, labelpad=8)

    ax.text(1.03, 0.70, "구분 안 됨", transform=ax.transAxes, fontsize=11,
            color="#b02020", fontweight="bold", va="center")
    ax.text(1.03, 0.32, f"p = {p:.2f}", transform=ax.transAxes, fontsize=9.5,
            color=INK_2, va="center")

    fig.text(0.17, 0.93, "Boltz-2도 이 변이를 감별하지 못했다",
             fontsize=13.5, color=C_NAVY, fontweight="bold", ha="left")
    fig.text(0.17, 0.855,
             f"시드 3 × 모델 5 = 각 {len(wt)}개 · 세로선 = 중앙값 · Mann–Whitney U 양측",
             fontsize=9, color=INK_2, ha="left")
    fig.text(0.17, 0.055,
             "복합체 예측 자체는 성공했다 (ipTM 0.85는 높은 값). 변이 감별력만 없다.  "
             "Bret et al. 2026의 결합부위 변이 둔감성 보고와 일치.",
             fontsize=8.5, color=INK_MUTED, ha="left")

    out = FIG_DIR / "slide14_boltz_iptm.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def main() -> int:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for f in (slide11(), slide13(), slide14()):
        print(f"  → {f.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
