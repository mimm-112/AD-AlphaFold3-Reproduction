# Boltz-2 실행 — TREM2 ± Aβ42

## 무엇을 왜 하는가

논문은 변이의 **구조** 변화를 봤다. 우리 재현 결과 그 변화는 노이즈 수준이었다.
실험 문헌은 TREM2 변이가 **Aβ 결합**을 떨어뜨린다고 보고한다.
→ 결합을 재야 한다. AlphaFold3는 못 재고 Boltz-2는 잰다.

## 입력

| 파일 | 내용 | 길이 |
|---|---|---|
| `TREM2_WT_AB42.yaml` | TREM2 야생형 + Aβ42 | 154 aa |
| `TREM2_R62H_AB42.yaml` | TREM2 R62H + Aβ42 | 154 aa |
| `TREM2_WT_alone.yaml` | TREM2 야생형 단독 (대조) | 112 aa |
| `TREM2_R62H_alone.yaml` | TREM2 R62H 단독 (대조) | 112 aa |

TREM2는 Ig 도메인 19–130 (결정구조 5UD8·5ELI 구간, 결합 필수 잔기 31–91 포함).
Aβ42는 APP P05067 672–713.

## 읽을 지표

- **ipTM** (0~1) — 두 사슬 사이 인터페이스 신뢰도. **이게 결합 지표다**
- **interface pLDDT** — 접촉면 잔기들의 확신도
- ❌ affinity 수치는 안 나온다. Boltz-2 affinity head는 **저분자 전용**이고
  Aβ42는 펩타이드다. `properties: affinity` 블록을 일부러 넣지 않았다.

## 실행 (Colab T4 권장)

```bash
!pip install -q "boltz[cuda]==2.2.1"
# 런타임 재시작 후

boltz predict inputs/boltz/TREM2_WT_AB42.yaml \
    --out_dir outputs/boltz --use_msa_server \
    --diffusion_samples 5 --seed 1 \
    --output_format mmcif
boltz predict inputs/boltz/TREM2_WT_AB42.yaml \
    --out_dir outputs/boltz --use_msa_server \
    --diffusion_samples 5 --seed 2 \
    --output_format mmcif
boltz predict inputs/boltz/TREM2_WT_AB42.yaml \
    --out_dir outputs/boltz --use_msa_server \
    --diffusion_samples 5 --seed 3 \
    --output_format mmcif
boltz predict inputs/boltz/TREM2_R62H_AB42.yaml \
    --out_dir outputs/boltz --use_msa_server \
    --diffusion_samples 5 --seed 1 \
    --output_format mmcif
boltz predict inputs/boltz/TREM2_R62H_AB42.yaml \
    --out_dir outputs/boltz --use_msa_server \
    --diffusion_samples 5 --seed 2 \
    --output_format mmcif
boltz predict inputs/boltz/TREM2_R62H_AB42.yaml \
    --out_dir outputs/boltz --use_msa_server \
    --diffusion_samples 5 --seed 3 \
    --output_format mmcif
```

## ⚠️ 해석 시 반드시 함께 보고할 것

Bret et al. 2026 (J Chem Inf Model)은 Boltz-2 affinity가 **결합부위 변이에 둔감**하다고
보고했다. 따라서 WT와 R62H의 ipTM 차이도 **노이즈 바닥 검정을 거쳐야** 한다.
시드 [1, 2, 3] × diffusion sample 5개로 분포를 만들어,
구조 파트와 동일한 Mann–Whitney U 검정을 적용한다.

## 근거 문헌

- Zhao et al. 2018, *Neuron* — TREM2가 Aβ 올리고머에 나노몰 결합, AD 변이가 결합 감소
- Zhong et al. 2018, *Mol Neurodegener* — oAβ1-42 고친화도 결합, 필수 잔기 31–91
- Yeh et al. 2016, *Neuron* — TREM2–APOE/CLU/LDL 결합, 질병 변이가 저해
- Passaro et al. 2025, *bioRxiv* — Boltz-2 (MIT 라이선스, FEP 근접 성능)
- Bret et al. 2026, *JCIM* — Boltz-2의 결합부위 변이 둔감성 (한계)
