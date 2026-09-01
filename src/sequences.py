"""GPX4... 가 아니라 AD 논문 재현용 서열 생성기.

Yao et al., Cell Genomics 4:100700 (2024) 의 AlphaFold 파트 재현 입력을 만든다.

하는 일:
  1. UniProt에서 CD33 / PILRA / TREM2 정본 서열을 받는다
  2. 변이 위치의 아미노산이 실제로 기대값인지 assert로 강제 확인한다
  3. 야생형 / 변이형 서열을 만든다
  4. 두 가지 구간으로 저장한다
       full : UniProt 전장  → 논문 그대로 재현 (A안, 본체)
       dom  : Ig-like V 도메인만 → 우리 개선판 (B안, 결정구조 비교용)

⚠️ 논문에는 AlphaFold 방법 서술이 없다 (STAR Methods에 항목 자체가 없음).
   어떤 구간을 넣었는지 명시되지 않았으므로 "전장"은 가장 자연스러운 해석이며,
   논문의 pTM=0.6 이 전장 막단백질에 걸맞은 값이라는 정황증거가 이를 뒷받침한다.

실행:  python3 src/sequences.py
출력:  inputs/seqs/*.fasta  (12개)
       inputs/seqs/ALL_for_alphafold_server.txt  (붙여넣기용 모음)
       results/sequences_manifest.csv
"""

from __future__ import annotations

import csv
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SEQ_DIR = ROOT / "inputs" / "seqs"
CACHE_DIR = SEQ_DIR / "_uniprot"
RESULTS_DIR = ROOT / "results"


# ---------------------------------------------------------------------------
# 타겟 정의
#
# 잔기 번호는 모두 UniProt 전구체(precursor) 번호 = dbSNP·논문과 같은 체계다.
# 아래 값은 UniProt 서열에 직접 인덱싱해 확인했다 (2026-08-31).
#
# domain  : UniProt이 정의한 Ig-like V-type 도메인 경계
# crystal : 참조 결정구조가 실제로 덮는 UniProt 구간 (RCSB 정렬 정보에서 확인)
# ---------------------------------------------------------------------------
TARGETS = [
    {
        "gene": "CD33",
        "acc": "P20138",
        "length": 364,
        "pos": 69, "wt": "R", "mut": "G",
        "rsid": "rs2455069",
        # 논문 Fig 4B 주석 문구 그대로
        "paper": {
            "snp_label": "pQTL rs2455069-A>G",
            "locus": "chr19:51225385",
            "wt_name": "Arginine", "wt_codon": "AGG", "wt_codon_hl": 0,
            "mut_name": "Glycine", "mut_codon": "GGG", "mut_codon_hl": 0,
            "figure": "Figure 4B",
        },
        "signal": (1, 17),
        "domain": (19, 135),        # Ig-like V-type
        "crystal": (21, 232),       # 5IHB/5J06 — Ig V + Ig C2 두 도메인을 다 덮는다
        "construct": (19, 135),     # ← 도메인 경계 채택. 이유는 아래 NOTE 참조
        "note": (
            "결정구조 5IHB/5J06은 21-232로 Ig 도메인 2개를 모두 포함한다. "
            "그 구간을 그대로 예측하면 두 도메인 사이 경첩 운동이 시드마다 흔들려 "
            "69번 국소 신호를 덮어버린다. 그래서 V 도메인만 쓰고, "
            "결정구조와의 비교는 겹치는 구간에서만 한다."
        ),
    },
    {
        "gene": "PILRA",
        "acc": "Q9UKJ1",
        "length": 303,
        # ⚠️ 문헌은 'G78R'이라 부르지만 UniProt 정본이 이미 R이다.
        #    rs1859788은 UniProt에 R->G 로 등재돼 있다. 만들 변이체는 78R -> 78G.
        "pos": 78, "wt": "R", "mut": "G",
        "rsid": "rs1859788",
        # 논문 Fig S5(a) 주석. 논문도 Arg를 기준, Gly를 변이로 적었다 → 우리 R78G 방향과 일치
        "paper": {
            "snp_label": "SNP rs1859788-A>G",
            "locus": "chr7:100374211",
            "wt_name": "Arginine", "wt_codon": "AGG", "wt_codon_hl": 0,
            "mut_name": "Glycine", "mut_codon": "GGG", "mut_codon_hl": 0,
            "figure": "Figure S5(a)",
        },
        "signal": (1, 19),
        "domain": (32, 150),
        "crystal": (32, 150),       # 4NFB(78R) / 3WUZ(78G) — 도메인 경계와 정확히 일치
        "construct": (32, 150),
        "note": (
            "결정구조 구간과 도메인 경계가 정확히 일치한다. "
            "4NFB(78R) vs 3WUZ(78G)가 이 변이의 실험적 정답이므로 직접 대조 가능."
        ),
    },
    {
        "gene": "TREM2",
        "acc": "Q9NZC2",
        "length": 230,
        "pos": 62, "wt": "R", "mut": "H",
        "rsid": "rs143332484",
        # 논문 Fig S5(b) 주석
        "paper": {
            "snp_label": "SNP rs143332484-C>T",
            "locus": "chr6:41161469",
            "wt_name": "Arginine", "wt_codon": "CGT", "wt_codon_hl": 1,
            "mut_name": "Histidine", "mut_codon": "CAT", "mut_codon_hl": 1,
            "figure": "Figure S5(b)",
        },
        "signal": (1, 18),
        "domain": (29, 112),
        "crystal": (19, 130),       # 5UD8/5ELI 구간. 5UD7은 19-174로 더 길다
        "construct": (19, 130),     # ← 결정구조 구간 채택 (단일 도메인 + 결정화된 양옆)
        "note": (
            "5UD8(R47H)/5ELI가 19-130을 덮는다. 단일 Ig 도메인이라 이 구간을 그대로 써도 "
            "경첩 문제가 없고, 5UD7(WT) vs 5UD8(R47H)을 점돌연변이 양성 대조군으로 쓸 수 있다. "
            "⚠️ UniProt은 R62H를 'does not affect protein structure'로 등재하고 있다."
        ),
    },
]


def fetch_uniprot(acc: str) -> str:
    """UniProt 정본 서열을 받는다. 한 번 받으면 캐시해서 재사용."""
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cached = CACHE_DIR / f"{acc}.fasta"

    if not cached.exists():
        url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
        print(f"  UniProt 다운로드: {acc}")
        with urllib.request.urlopen(url, timeout=30) as r:
            cached.write_text(r.read().decode())
    else:
        print(f"  캐시 사용: {acc}")

    return "".join(
        line.strip()
        for line in cached.read_text().splitlines()
        if not line.startswith(">")
    )


def apply_mutation(seq: str, pos: int, wt: str, mut: str, gene: str) -> str:
    """1-based 전구체 번호로 치환. WT가 실제로 그 자리에 있는지 강제 확인.

    이 assert가 깨지면 isoform이나 번호 체계가 틀린 것이다.
    (PILRA에 'G78R'을 적용하려 하면 여기서 걸린다 — 78번은 이미 R이다.)
    """
    actual = seq[pos - 1]
    if actual != wt:
        raise AssertionError(
            f"{gene} {pos}번이 '{actual}'인데 '{wt}'를 기대함. "
            f"isoform 또는 번호 체계 오류 — CLAUDE.md의 잔기 번호 표를 확인할 것."
        )
    return seq[: pos - 1] + mut + seq[pos:]


def write_fasta(path: Path, header: str, seq: str, width: int = 60) -> None:
    lines = [f">{header}"]
    lines += [seq[i : i + width] for i in range(0, len(seq), width)]
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    SEQ_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    rows = []
    bundle = []

    for t in TARGETS:
        gene, acc = t["gene"], t["acc"]
        pos, wt, mut = t["pos"], t["wt"], t["mut"]
        label = f"{wt}{pos}{mut}"

        print(f"\n[{gene}] {acc}  {label}  ({t['rsid']})")
        seq_wt = fetch_uniprot(acc)

        # --- 검증 ---------------------------------------------------------
        if len(seq_wt) != t["length"]:
            raise AssertionError(
                f"{gene} 길이가 {len(seq_wt)}인데 {t['length']}를 기대함. "
                f"UniProt 서열 버전이 바뀌었을 수 있으니 확인할 것."
            )
        seq_mut = apply_mutation(seq_wt, pos, wt, mut, gene)
        print(f"  ✓ 길이 {len(seq_wt)} / {pos}번 = {wt} 확인 → {mut}로 치환")

        start, end = t["construct"]
        # 변이 위치가 잘라낸 구간 안에 있는지 확인
        if not (start <= pos <= end):
            raise AssertionError(f"{gene} 변이 {pos}번이 구간 {start}-{end} 밖에 있다")

        variants = [
            # (구간 태그, 설명, WT 서열, 변이 서열, 구간)
            ("full", "논문 그대로 재현 (A안)", seq_wt, seq_mut, (1, len(seq_wt))),
            ("dom", "Ig 도메인 절단 (B안)", seq_wt[start - 1 : end], seq_mut[start - 1 : end], (start, end)),
        ]

        for tag, desc, s_wt, s_mut, (a, b) in variants:
            for allele, s in (("WT", s_wt), (label, s_mut)):
                name = f"{gene}_{tag}_{allele}"
                header = (
                    f"{name} | {acc} {a}-{b} ({len(s)} aa) | "
                    f"{t['rsid']} {label} | {desc}"
                )
                path = SEQ_DIR / f"{name}.fasta"
                write_fasta(path, header, s)

                bundle.append(f"### {name}   [{desc}]  {len(s)} aa\n{s}\n")
                rows.append({
                    "name": name,
                    "gene": gene,
                    "uniprot": acc,
                    "rsid": t["rsid"],
                    "mutation": label,
                    "allele": allele,
                    "construct": tag,
                    "range_start": a,
                    "range_end": b,
                    "length": len(s),
                    # 잘라낸 서열 안에서 변이가 몇 번째인지 (1-based) — PyMOL/RMSD용
                    "mut_index_in_construct": pos - a + 1,
                    "fasta": str(path.relative_to(ROOT)),
                })

            print(f"  → {tag}: {a}-{b} ({b - a + 1} aa) WT + {label}")

        print(f"  NOTE: {t['note']}")

    # --- 매니페스트 -------------------------------------------------------
    manifest = RESULTS_DIR / "sequences_manifest.csv"
    with manifest.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    # --- 붙여넣기용 모음 --------------------------------------------------
    bundle_path = SEQ_DIR / "ALL_for_alphafold_server.txt"
    bundle_path.write_text(
        "AlphaFold Server 제출용 서열 모음\n"
        "https://alphafoldserver.com — 구글 계정 로그인, 하루 30 job, job당 모델 5개\n"
        "\n"
        "제출 순서 권장:\n"
        "  1) _full_ 6개 먼저 = 논문 그대로 재현 (A안)\n"
        "  2) _dom_  6개 그다음 = 결정구조 비교용 (B안)\n"
        "  ※ 시드는 자동생성 대신 직접 입력해 기록해 둘 것 (노이즈 바닥 분석에 필요)\n"
        "\n" + "=" * 70 + "\n\n" + "\n".join(bundle)
    )

    print("\n" + "=" * 62)
    print(f"FASTA {len(rows)}개 → {SEQ_DIR.relative_to(ROOT)}/")
    print(f"붙여넣기용 모음 → {bundle_path.relative_to(ROOT)}")
    print(f"매니페스트     → {manifest.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
