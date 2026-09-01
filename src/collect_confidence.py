"""AlphaFold Server 결과에서 신뢰도 지표를 모아 표로 만든다.

AF Server 출력 폴더(folds_YYYY_MM_DD_HH_MM/)를 읽어
job × model 별로 pTM, pLDDT, 변이부위 pLDDT 등을 뽑는다.

논문 대조 포인트:
  Fig 4B 캡션 — CD33 두 구조 모두 "predicted template modeling yields a score of 0.6"
  → CD33_full_WT / CD33_full_R69G 의 ptm이 0.6 근처면 논문 재현 성공.

⚠️ AF Server의 ptm은 0~1 스케일, CIF의 B-factor(pLDDT)는 0~100 스케일이다.

실행:
  python3 src/collect_confidence.py                      # 최신 folds_* 폴더 자동 탐색
  python3 src/collect_confidence.py folds_2026_08_31_16_46
"""

from __future__ import annotations

import csv
import json
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sequences import TARGETS  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR = ROOT / "results"

# gene 이름(소문자) → 타겟 정의
BY_GENE = {t["gene"].lower(): t for t in TARGETS}


def parse_cif_ca_plddt(cif_path: Path) -> dict[int, float]:
    """AF3 mmCIF에서 잔기별 pLDDT를 뽑는다.

    AlphaFold는 pLDDT를 B-factor 칼럼에 넣는다. CA 원자만 취하면 잔기당 하나가 된다.
    biopython 없이 atom_site 루프만 직접 읽는다.
    """
    cols: list[str] = []
    in_loop = False
    out: dict[int, float] = {}

    for line in cif_path.read_text().splitlines():
        s = line.strip()

        if s.startswith("_atom_site."):
            cols.append(s.split(".", 1)[1])
            in_loop = True
            continue

        if in_loop and (s.startswith("_") or s == "loop_" or s == "#"):
            if cols and out:
                break          # atom_site 루프가 끝났다
            in_loop = bool(cols)
            continue

        if not in_loop or not cols or not s or s.startswith("#"):
            continue

        f = s.split()
        if len(f) < len(cols):
            continue

        rec = dict(zip(cols, f))
        if rec.get("group_PDB") != "ATOM":
            continue
        if rec.get("label_atom_id") != "CA":
            continue

        try:
            resid = int(rec.get("label_seq_id", rec.get("auth_seq_id", "0")))
            out[resid] = float(rec["B_iso_or_equiv"])
        except (ValueError, KeyError):
            continue

    return out


def find_latest_folds_dir() -> Path | None:
    dirs = sorted(ROOT.glob("folds_*"), reverse=True)
    return dirs[0] if dirs else None


def parse_job_name(job: str) -> tuple[str, str, str]:
    """'cd33_full_r69g' → ('CD33', 'full', 'R69G')"""
    parts = job.split("_")
    gene = parts[0].upper()
    construct = parts[1] if len(parts) > 1 else "?"
    allele = "WT" if parts[-1].lower() == "wt" else parts[-1].upper()
    return gene, construct, allele


def main(folds_dir: Path) -> int:
    if not folds_dir.exists():
        print(f"폴더 없음: {folds_dir}")
        return 1

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    rows = []

    for job_dir in sorted(p for p in folds_dir.iterdir() if p.is_dir()):
        job = job_dir.name
        gene, construct, allele = parse_job_name(job)
        target = BY_GENE.get(gene.lower())

        for conf_path in sorted(job_dir.glob("*_summary_confidences_*.json")):
            model = int(conf_path.stem.rsplit("_", 1)[1])
            conf = json.loads(conf_path.read_text())

            cif = job_dir / f"fold_{job}_model_{model}.cif"
            plddt = parse_cif_ca_plddt(cif) if cif.exists() else {}

            row = {
                "job": job,
                "gene": gene,
                "construct": construct,
                "allele": allele,
                "model": model,
                "ptm": conf.get("ptm"),
                "ranking_score": conf.get("ranking_score"),
                "fraction_disordered": conf.get("fraction_disordered"),
                "has_clash": conf.get("has_clash"),
                "n_res": len(plddt),
                "mean_plddt": round(statistics.mean(plddt.values()), 2) if plddt else None,
            }

            if target and plddt:
                pos = target["pos"]
                # 전장이면 UniProt 번호 그대로, 도메인 절단이면 오프셋을 뺀다
                idx = pos if construct == "full" else pos - target["construct"][0] + 1
                row["mut_pos_uniprot"] = pos
                row["plddt_at_mut"] = plddt.get(idx)

                d0, d1 = target["domain"]
                if construct == "dom":
                    d0, d1 = d0 - target["construct"][0] + 1, d1 - target["construct"][0] + 1
                dom_vals = [v for k, v in plddt.items() if d0 <= k <= d1]
                row["mean_plddt_domain"] = round(statistics.mean(dom_vals), 2) if dom_vals else None

            rows.append(row)

    if not rows:
        print("결과 파일을 못 찾았다.")
        return 1

    out_csv = RESULTS_DIR / "confidence_summary.csv"
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # ---- 화면 출력 -------------------------------------------------------
    print(f"입력: {folds_dir.name}   ({len(rows)}개 모델)\n")
    hdr = (f"{'job':22s} {'model':>5s} {'pTM':>5s} {'rank':>5s} "
           f"{'disord':>6s} {'pLDDT평균':>9s} {'도메인':>7s} {'변이부위':>8s}")
    print(hdr)
    print("-" * len(hdr))

    for r in rows:
        print(f"{r['job']:22s} {r['model']:>5d} {r['ptm']:>5.2f} "
              f"{r['ranking_score']:>5.2f} {r['fraction_disordered']:>6.2f} "
              f"{r['mean_plddt'] or 0:>9.1f} "
              f"{r.get('mean_plddt_domain') or 0:>7.1f} "
              f"{r.get('plddt_at_mut') or 0:>8.1f}")

    # ---- 논문 대조 -------------------------------------------------------
    print("\n" + "=" * 62)
    print("논문 대조: Fig 4B — CD33 두 구조 모두 pTM = 0.6")
    for allele in ("WT", "R69G"):
        vals = [r["ptm"] for r in rows
                if r["gene"] == "CD33" and r["construct"] == "full" and r["allele"] == allele]
        if vals:
            print(f"  CD33_full_{allele:5s} pTM = {min(vals):.2f} ~ {max(vals):.2f} "
                  f"(모델 {len(vals)}개, 중앙값 {statistics.median(vals):.2f})")

    print(f"\n저장: {out_csv.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    d = Path(sys.argv[1]) if len(sys.argv) > 1 else find_latest_folds_dir()
    if d and not d.is_absolute():
        d = ROOT / d
    sys.exit(main(d) if d else 1)
