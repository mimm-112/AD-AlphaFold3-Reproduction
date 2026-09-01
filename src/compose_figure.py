"""PyMOL 렌더 결과를 논문 Figure 4B / S5 레이아웃으로 합성한다.

pymol_render.py가 만든 두 이미지(전체 오버레이 + 변이부위 확대)를
논문과 같은 구도로 붙인다: 왼쪽 큰 구조, 오른쪽 아래 검은 테두리 인셋, 왼쪽 위 주석.

논문 주석 문구는 sequences.py의 TARGETS[*]["paper"] 에 그대로 옮겨 두었다.
  CD33  : pQTL rs2455069-A>G  (chr19:51225385)  Arginine [AGG] > Glycine [GGG]
  PILRA : SNP rs1859788-A>G   (chr7:100374211)  Arginine [AGG] > Glycine [GGG]
  TREM2 : SNP rs143332484-C>T (chr6:41161469)   Arginine [CGT] > Histidine [CAT]

⚠️ pymol 환경에는 matplotlib이 없으므로 이 스크립트는 base 환경에서 돌린다.
   렌더:  /opt/anaconda3/envs/pymol/bin/python src/pymol_render.py
   합성:  python3 src/compose_figure.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg          # noqa: E402
import matplotlib.pyplot as plt           # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from sequences import TARGETS  # noqa: E402

FIG_DIR = ROOT / "figures"

# 논문이 쓴 색 (본문 색 강조용)
C_WT = "#00A000"    # Arginine — 초록
C_MUT = "#C8C800"   # Glycine/Histidine — 노랑(가독성 위해 약간 어둡게)
C_REF = "#0000FF"   # 참조 대립 — 파랑
C_ALT = "#FF0000"   # 변이 대립 — 빨강


def trim_white(img, pad: int = 8):
    """흰 여백을 잘라 구조가 화면을 채우게 한다."""
    import numpy as np
    a = img[..., :3] if img.ndim == 3 else img
    mask = (a < 0.98).any(axis=2) if a.ndim == 3 else (a < 0.98)
    if not mask.any():
        return img
    ys, xs = np.where(mask)
    y0, y1 = max(ys.min() - pad, 0), min(ys.max() + pad + 1, img.shape[0])
    x0, x1 = max(xs.min() - pad, 0), min(xs.max() + pad + 1, img.shape[1])
    return img[y0:y1, x0:x1]


def compose(target: dict) -> Path | None:
    gene = target["gene"]
    p = target["paper"]
    label = f"{target['wt']}{target['pos']}{target['mut']}"

    overlay_p = FIG_DIR / f"{gene}_overlay.png"
    zoom_p = FIG_DIR / f"{gene}_zoom.png"
    if not (overlay_p.exists() and zoom_p.exists()):
        print(f"  건너뜀 — 렌더 결과 없음: {gene}")
        return None

    overlay = trim_white(mpimg.imread(overlay_p))
    zoom = trim_white(mpimg.imread(zoom_p))

    fig = plt.figure(figsize=(10, 7.5), dpi=200)
    fig.patch.set_facecolor("white")

    # --- 본 구조 ----------------------------------------------------------
    ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
    ax.imshow(overlay)
    ax.axis("off")

    # --- 주석 (논문 왼쪽 위 텍스트 블록) -----------------------------------
    def colored_line(x, y, segments, size=10.5):
        """(문자열, 색, 굵기) 조각들을 겹치지 않게 이어 그린다.

        글자 폭을 렌더러로 실제 측정해서 다음 조각의 x를 정한다.
        (문자수 × 상수로 어림하면 폰트마다 어긋난다.)
        """
        renderer = fig.canvas.get_renderer()
        cx = x
        for text, color, bold in segments:
            t = ax.text(cx, y, text, transform=ax.transAxes, fontsize=size,
                        color=color, va="top",
                        fontweight="bold" if bold else "normal")
            bb = t.get_window_extent(renderer=renderer)
            cx += bb.width / fig.bbox.width

    fig.canvas.draw()  # 렌더러 준비
    tx, ty = 0.035, 0.965

    ax.text(tx, ty, gene, transform=ax.transAxes, fontsize=17,
            fontweight="bold", fontstyle="italic", va="top")

    ref, alt = p["snp_label"].rsplit("-", 1)[1].split(">")
    colored_line(tx, ty - 0.075, [
        (p["snp_label"].rsplit("-", 1)[0] + "-", "black", False),
        (ref, C_REF, True),
        (">", "black", False),
        (alt, C_ALT, True),
    ])

    ax.text(tx, ty - 0.115, f"({p['locus']})", transform=ax.transAxes,
            fontsize=10.5, va="top")

    ax.text(tx, ty - 0.185, "Amino Acid [Codon]:", transform=ax.transAxes,
            fontsize=10.5, va="top")
    colored_line(tx, ty - 0.225, [
        (f"{p['wt_name']} [{p['wt_codon']}]", C_WT, True),
        (" > ", "black", False),
        (f"{p['mut_name']} [{p['mut_codon']}]", C_MUT, True),
    ])

    # --- 확대 인셋 (논문처럼 검은 테두리 상자) -----------------------------
    ins = fig.add_axes([0.50, 0.045, 0.47, 0.44])
    ins.imshow(zoom)
    ins.set_xticks([]); ins.set_yticks([])
    for s in ins.spines.values():
        s.set_linewidth(2.5)
        s.set_edgecolor("black")

    # --- 재현 표기 (인셋과 겹치지 않게 왼쪽 아래) --------------------------
    ax.text(0.02, 0.018,
            f"reproduction of {p['figure']}\n"
            f"AlphaFold Server (AF3) · full-length UniProt {target['acc']} · model 0",
            transform=ax.transAxes, fontsize=7.5, color="#555555",
            ha="left", va="bottom", linespacing=1.5)

    out = FIG_DIR / f"{gene}_{label}_paper_style.png"
    fig.savefig(out, dpi=200, facecolor="white")
    plt.close(fig)
    print(f"  {gene:6s} {label:6s} → {out.name}")
    return out


def main() -> int:
    if not FIG_DIR.exists():
        print("figures/ 가 없다. 먼저 pymol_render.py 를 돌릴 것.")
        return 1

    print("논문 레이아웃으로 합성 (오버레이 + 확대 인셋 + 주석)\n")
    made = [compose(t) for t in TARGETS]
    made = [m for m in made if m]

    if not made:
        return 1
    print(f"\n완성 {len(made)}장 → figures/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
