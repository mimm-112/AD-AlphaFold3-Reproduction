"""논문 Figure 4B / Figure S5 를 우리 예측 결과로 재현한다.

Yao et al., Cell Genomics 4:100700 (2024) 의 그림 설명을 그대로 따른다.

  Fig 4B 캡션 원문:
    "The ribbon representation of 3D structures of CD33 with arginine and glycine
     at position 69 are colored in blue and red, respectively. The amino acids at
     position 69 are displayed in stick representation, with arginine and glycine
     colored in green and yellow, respectively."

  → 야생형(Arg) 리본 = 파랑 / 변이형(Gly·His) 리본 = 빨강
  → 변이 위치 잔기만 stick, 야생형 잔기 = 초록 / 변이 잔기 = 노랑

논문은 어느 모델을 썼는지 밝히지 않았다. 우리는 AF Server가 1순위로 내놓는
model_0 을 쓴다 (ranking_score 최상위).

실행:
  /opt/anaconda3/envs/pymol/bin/python src/pymol_render.py
출력:
  figures/<GENE>_overlay.png        전체 오버레이 (논문 Fig 4B 구도)
  figures/<GENE>_zoom.png           변이 부위 확대
  results/rmsd_overlay.csv          정렬 시 나온 RMSD (참고용)
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from sequences import TARGETS  # noqa: E402

FOLDS = sorted(ROOT.glob("folds_*"), reverse=True)
FIG_DIR = ROOT / "figures"
RESULTS_DIR = ROOT / "results"

# 논문 Fig 4B 색 지정
COLOR_WT_RIBBON = "blue"
COLOR_MUT_RIBBON = "red"
COLOR_WT_RESIDUE = "green"    # arginine
COLOR_MUT_RESIDUE = "yellow"  # glycine / histidine

MODEL_INDEX = 0  # AF Server 1순위 모델

from pymol import cmd  # noqa: E402


def cif_for(job: str, model: int = MODEL_INDEX) -> Path:
    return FOLDS[0] / job / f"fold_{job}_model_{model}.cif"


def render_pair(target: dict) -> dict:
    gene = target["gene"]
    pos = target["pos"]
    wt, mut = target["wt"], target["mut"]
    label = f"{wt}{pos}{mut}"

    wt_job = f"{gene.lower()}_full_wt"
    mut_job = f"{gene.lower()}_full_{label.lower()}"

    wt_cif, mut_cif = cif_for(wt_job), cif_for(mut_job)
    if not (wt_cif.exists() and mut_cif.exists()):
        print(f"  건너뜀 — 파일 없음: {wt_cif.name} / {mut_cif.name}")
        return {}

    cmd.reinitialize()
    cmd.load(str(wt_cif), "WT")
    cmd.load(str(mut_cif), "MUT")

    # --- 정렬 --------------------------------------------------------------
    # 전장에는 흐물흐물한 꼬리가 있어 그걸 포함해 맞추면 도메인이 어긋난다.
    # 논문은 정렬 범위를 밝히지 않았으므로, 구조가 잡힌 Ig 도메인 기준으로 맞춘다.
    d0, d1 = target["domain"]
    sel = f"resi {d0}-{d1} and name CA and polymer"
    rms_dom = cmd.align(f"MUT and {sel}", f"WT and {sel}", cycles=0)[0]
    rms_all = cmd.rms_cur("MUT and name CA", "WT and name CA", matchmaker=-1)

    # --- 논문 색 지정 ------------------------------------------------------
    cmd.hide("everything")
    cmd.show("cartoon", "polymer")
    cmd.color(COLOR_WT_RIBBON, "WT")
    cmd.color(COLOR_MUT_RIBBON, "MUT")

    wt_res = f"WT and resi {pos}"
    mut_res = f"MUT and resi {pos}"
    cmd.show("sticks", f"({wt_res} or {mut_res}) and not hydro")
    cmd.color(COLOR_WT_RESIDUE, wt_res)
    cmd.color(COLOR_MUT_RESIDUE, mut_res)

    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("cartoon_transparency", 0.1)
    cmd.set("ray_shadows", 0)
    cmd.set("antialias", 2)
    cmd.set("stick_radius", 0.22)

    FIG_DIR.mkdir(parents=True, exist_ok=True)

    # --- 전체 오버레이 (논문 Fig 4B 구도) ----------------------------------
    cmd.orient("polymer")
    cmd.zoom("polymer", 3)
    overlay = FIG_DIR / f"{gene}_overlay.png"
    cmd.png(str(overlay), width=1600, height=1200, dpi=300, ray=1)

    # --- 변이 부위 확대 ----------------------------------------------------
    cmd.orient(f"({wt_res} or {mut_res})")
    cmd.zoom(f"({wt_res} or {mut_res})", 6)
    zoom = FIG_DIR / f"{gene}_zoom.png"
    cmd.png(str(zoom), width=1600, height=1200, dpi=300, ray=1)

    print(f"  {gene:6s} {label:6s}  도메인 RMSD {rms_dom:.3f} Å  "
          f"전장 RMSD {rms_all:.3f} Å  →  {overlay.name}, {zoom.name}")

    return {
        "gene": gene,
        "mutation": label,
        "wt_model": wt_cif.name,
        "mut_model": mut_cif.name,
        "align_range": f"{d0}-{d1}",
        "rmsd_domain_CA": round(rms_dom, 3),
        "rmsd_fulllength_CA": round(rms_all, 3),
        "overlay_png": str(overlay.relative_to(ROOT)),
        "zoom_png": str(zoom.relative_to(ROOT)),
    }


def main() -> int:
    if not FOLDS:
        print("folds_* 폴더가 없다. AlphaFold Server 결과를 먼저 넣을 것.")
        return 1

    print(f"입력: {FOLDS[0].name}   (모델 {MODEL_INDEX}번 사용)")
    print("색: 야생형 리본 파랑 / 변이형 리본 빨강, 변이잔기 초록·노랑 (논문 Fig 4B와 동일)\n")

    rows = [r for r in (render_pair(t) for t in TARGETS) if r]
    if not rows:
        return 1

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out = RESULTS_DIR / "rmsd_overlay.csv"
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    print(f"\n그림 {len(rows) * 2}장 → figures/")
    print(f"수치       → {out.relative_to(ROOT)}")
    print("\n⚠️ RMSD는 논문에 없는 값이다. 정렬 품질 확인용으로만 기록했다.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
