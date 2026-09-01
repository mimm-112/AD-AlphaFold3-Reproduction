"""AlphaFold Server 업로드용 JSON 생성기.

sequences.py가 만든 매니페스트를 읽어 AF Server가 받는 JSON 배치 파일을 만든다.
웹페이지의 [Upload JSON] 버튼으로 올리면 여러 job이 한 번에 등록된다.

AF Server JSON 형식 (dialect: "alphafoldserver", version: 1):
  최상위가 리스트여야 여러 job을 한 파일에 담을 수 있다.

  [
    {
      "name": "job 이름",
      "modelSeeds": [1],
      "sequences": [{"proteinChain": {"sequence": "...", "count": 1}}],
      "dialect": "alphafoldserver",
      "version": 1
    }
  ]

⚠️ modelSeeds를 빈 리스트로 두면 서버가 랜덤 시드를 자동 배정한다.
   우리는 노이즈 바닥 분석을 해야 하므로 시드를 명시해 기록을 남긴다.

실행:
  python3 src/make_af_jobs.py            # A안(전장) + B안(도메인), 시드 1
  python3 src/make_af_jobs.py 1 2 3 4    # 시드 4개짜리 앙상블용까지 생성
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MANIFEST = ROOT / "results" / "sequences_manifest.csv"
OUT_DIR = ROOT / "inputs" / "af_jobs"

# AF Server 하루 job 한도 (2026-01 기준). 넘으면 경고만 띄운다.
DAILY_JOB_LIMIT = 30


def read_seq(rel_path: str) -> str:
    p = ROOT / rel_path
    return "".join(
        line.strip()
        for line in p.read_text().splitlines()
        if not line.startswith(">")
    )


def make_job(name: str, sequence: str, seed: int) -> dict:
    return {
        "name": name,
        "modelSeeds": [seed],
        "sequences": [{"proteinChain": {"sequence": sequence, "count": 1}}],
        "dialect": "alphafoldserver",
        "version": 1,
    }


def main(seeds: list[int]) -> int:
    if not MANIFEST.exists():
        print("매니페스트가 없다. 먼저 실행:  python3 src/sequences.py")
        return 1

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = list(csv.DictReader(MANIFEST.open()))

    plans = [
        ("A_full", "full", "논문 그대로 재현 (전장 서열)"),
        ("B_dom", "dom", "결정구조 비교용 (Ig 도메인)"),
    ]

    total = 0
    for tag, construct, desc in plans:
        subset = [r for r in rows if r["construct"] == construct]

        jobs = []
        for seed in seeds:
            for r in subset:
                # 시드가 여러 개면 job 이름에 시드를 붙여 결과를 구분한다
                name = r["name"] if len(seeds) == 1 else f"{r['name']}_s{seed}"
                jobs.append(make_job(name, read_seq(r["fasta"]), seed))

        suffix = f"seed{seeds[0]}" if len(seeds) == 1 else f"seeds{min(seeds)}-{max(seeds)}"
        out = OUT_DIR / f"af_jobs_{tag}_{suffix}.json"
        out.write_text(json.dumps(jobs, indent=2) + "\n")

        total += len(jobs)
        print(f"{out.relative_to(ROOT)}")
        print(f"   {desc}")
        print(f"   job {len(jobs)}개  →  구조 {len(jobs) * 5}개 (job당 모델 5개)")
        for j in jobs[:3]:
            n = len(j["sequences"][0]["proteinChain"]["sequence"])
            print(f"     · {j['name']}  ({n} aa, seed {j['modelSeeds'][0]})")
        if len(jobs) > 3:
            print(f"     · ... 외 {len(jobs) - 3}개")
        print()

    print(f"합계 job {total}개")
    if total > DAILY_JOB_LIMIT:
        print(
            f"⚠️ 하루 한도 {DAILY_JOB_LIMIT}개를 넘는다. "
            f"파일을 나눠서 며칠에 걸쳐 올릴 것."
        )
    else:
        print(f"하루 한도 {DAILY_JOB_LIMIT}개 이내 — 한 번에 올려도 된다.")

    print("\n올리는 법:")
    print("  1. https://alphafoldserver.com 접속")
    print("  2. 오른쪽 위 [Upload JSON] 클릭")
    print(f"  3. {OUT_DIR.relative_to(ROOT)}/af_jobs_A_full_*.json 선택  ← 논문 재현부터")
    print("  4. job 목록이 뜨면 제출")
    return 0


if __name__ == "__main__":
    args = [int(a) for a in sys.argv[1:]] or [1]
    sys.exit(main(args))
