"""모델 쌍별 RMSD를 전부 계산해 CSV로 낸다. (PyMOL 환경에서 실행)

AF Server는 job당 모델 5개를 준다. 이걸 이용해 두 가지를 잰다.

  노이즈(noise)  : 같은 서열끼리의 모델 쌍  WT-WT 10쌍, MUT-MUT 10쌍
  신호(signal)   : 다른 서열끼리의 모델 쌍  WT-MUT 25쌍

같은 서열인데도 모델마다 구조가 조금씩 다르다. 그 차이가 "노이즈 바닥"이고,
변이 때문에 생긴 차이(신호)가 그 바닥을 넘는지 보면 된다.

⚠️ 이 분석은 논문에 없다. 논문은 WT 1개 vs 변이체 1개를 겹쳐 보기만 했다.
   다만 추가 예측 없이 이미 받은 결과 30개만으로 계산하므로 비용은 0이다.

정렬 방식
  · 전장에는 흐물흐물한 꼬리가 있어 그걸 포함해 맞추면 도메인이 어긋난다.
  · 그래서 Ig 도메인 CA로 정렬한 뒤, 그 상태에서 RMSD를 잰다.
  · local RMSD는 정렬을 다시 하지 않는다. 다시 맞추면 국소 차이가 정의상 0이 된다.

실행: /opt/anaconda3/envs/pymol/bin/python src/rmsd_analysis.py
출력: results/rmsd_pairs.csv
"""

from __future__ import annotations

import csv
import itertools
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "src"))
from sequences import TARGETS  # noqa: E402

from pymol import cmd  # noqa: E402

FOLDS = sorted(ROOT.glob("folds_*"), reverse=True)
N_MODELS = 5
LOCAL_CUTOFF = 8.0


def cif(job: str, model: int) -> Path:
    return FOLDS[0] / job / f"fold_{job}_model_{model}.cif"


def rmsd(cif_a: Path, cif_b: Path, domain: tuple[int, int], pos: int) -> tuple[float, float]:
    """도메인 정렬 후 (도메인 RMSD, 변이부위 국소 RMSD)."""
    # ⚠️ 객체 이름을 "A"/"B"로 두면 체인 ID A와 충돌해 align이 실패한다.
    cmd.reinitialize()
    cmd.load(str(cif_a), "ref")
    cmd.load(str(cif_b), "mob")

    d0, d1 = domain
    sel = f"resi {d0}-{d1} and name CA and polymer"
    dom = cmd.align(f"mob and {sel}", f"ref and {sel}", cycles=0)[0]

    # 정렬을 유지한 채 변이부위 주변만 다시 잰다 (재정렬 금지)
    cmd.select("loc", f"ref and byres (polymer within {LOCAL_CUTOFF} of "
                      f"(ref and resi {pos})) and name CA")
    ids = sorted({at.resi for at in cmd.get_model("loc").atom}, key=int)
    if ids:
        s = "resi " + "+".join(ids) + " and name CA"
        loc = cmd.rms_cur(f"mob and {s}", f"ref and {s}", matchmaker=-1)
    else:
        loc = float("nan")
    return dom, loc


def main() -> int:
    if not FOLDS:
        print("folds_* 폴더가 없다.")
        return 1

    rows = []
    for t in TARGETS:
        gene = t["gene"]
        label = f"{t['wt']}{t['pos']}{t['mut']}"
        wt_job = f"{gene.lower()}_full_wt"
        mut_job = f"{gene.lower()}_full_{label.lower()}"

        if not cif(wt_job, 0).exists():
            print(f"  건너뜀: {gene}")
            continue

        pairs = []
        # 노이즈: 같은 서열 안에서 모델끼리
        for i, j in itertools.combinations(range(N_MODELS), 2):
            pairs.append(("noise_WT", cif(wt_job, i), cif(wt_job, j), i, j))
            pairs.append(("noise_MUT", cif(mut_job, i), cif(mut_job, j), i, j))
        # 신호: WT 모델 × 변이체 모델 전조합
        for i in range(N_MODELS):
            for j in range(N_MODELS):
                pairs.append(("signal_WT_vs_MUT", cif(wt_job, i), cif(mut_job, j), i, j))

        for kind, a, b, i, j in pairs:
            dom, loc = rmsd(a, b, t["domain"], t["pos"])
            rows.append({
                "gene": gene, "mutation": label, "comparison": kind,
                "model_a": i, "model_b": j,
                "rmsd_domain": round(dom, 4),
                "rmsd_local8A": round(loc, 4),
            })

        n_noise = sum(1 for r in rows if r["gene"] == gene and r["comparison"].startswith("noise"))
        n_sig = sum(1 for r in rows if r["gene"] == gene and r["comparison"] == "signal_WT_vs_MUT")
        print(f"  {gene:6s} 노이즈 {n_noise}쌍 · 신호 {n_sig}쌍 계산 완료")

    out = ROOT / "results" / "rmsd_pairs.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"\n총 {len(rows)}쌍 → {out.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
