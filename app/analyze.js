/* 구조 분석 모듈 — 업로드된 예측 결과를 브라우저에서 직접 계산한다.
 *
 *   mmCIF 파싱 → Cα 추출 (pLDDT = B-factor)
 *   Kabsch/Horn 최적 중첩 → 전역 정렬
 *   변이 잔기 8 Å 이내 국소 RMSD (재정렬 없이)
 *
 * src/rmsd_analysis.py 와 동일한 프로토콜을 이식한 것이다.
 */

/* ── mmCIF 파싱 ───────────────────────────────────────────────
   AlphaFold·Boltz 출력의 atom_site 루프만 읽는다.
   AlphaFold 는 pLDDT 를 B-factor 칼럼에 넣는다. */
function parseCifCA(text){
  const lines = text.split('\n');
  let cols = [], inLoop = false;
  const out = {};                       // resid → {xyz:[x,y,z], plddt}

  for(const raw of lines){
    const s = raw.trim();
    if(s.startsWith('_atom_site.')){ cols.push(s.split('.')[1]); inLoop = true; continue; }
    if(inLoop && (s === '#' || s === 'loop_' || s.startsWith('_'))){
      if(cols.length && Object.keys(out).length) break;
      inLoop = cols.length > 0; continue;
    }
    if(!inLoop || !cols.length || !s || s.startsWith('#')) continue;

    const f = s.split(/\s+/);
    if(f.length < cols.length) continue;
    const rec = {};
    cols.forEach((c,i)=>rec[c]=f[i]);
    if(rec.group_PDB !== 'ATOM') continue;
    if(rec.label_atom_id !== 'CA') continue;

    const resid = parseInt(rec.label_seq_id ?? rec.auth_seq_id, 10);
    if(!Number.isFinite(resid)) continue;
    out[resid] = {
      xyz: [+rec.Cartn_x, +rec.Cartn_y, +rec.Cartn_z],
      plddt: +rec.B_iso_or_equiv
    };
  }
  return out;
}

/* ── 4×4 대칭행렬 고유분해 (Jacobi) ────────────────────────── */
function jacobi4(A){
  const n = 4;
  let a = A.map(r=>r.slice());
  let v = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]];
  for(let sweep=0; sweep<60; sweep++){
    let off = 0;
    for(let p=0;p<n;p++) for(let q=p+1;q<n;q++) off += a[p][q]*a[p][q];
    if(off < 1e-14) break;
    for(let p=0;p<n;p++) for(let q=p+1;q<n;q++){
      if(Math.abs(a[p][q]) < 1e-16) continue;
      const theta = (a[q][q]-a[p][p]) / (2*a[p][q]);
      const t = Math.sign(theta||1) / (Math.abs(theta)+Math.sqrt(theta*theta+1));
      const c = 1/Math.sqrt(t*t+1), s = t*c;
      for(let k=0;k<n;k++){
        const akp=a[k][p], akq=a[k][q];
        a[k][p]=c*akp-s*akq; a[k][q]=s*akp+c*akq;
      }
      for(let k=0;k<n;k++){
        const apk=a[p][k], aqk=a[q][k];
        a[p][k]=c*apk-s*aqk; a[q][k]=s*apk+c*aqk;
      }
      for(let k=0;k<n;k++){
        const vkp=v[k][p], vkq=v[k][q];
        v[k][p]=c*vkp-s*vkq; v[k][q]=s*vkp+c*vkq;
      }
    }
  }
  return {values: [a[0][0],a[1][1],a[2][2],a[3][3]], vectors: v};
}

/* ── 최적 중첩 (Horn 사원수법) ──────────────────────────────
   mob 을 ref 에 맞추는 회전 R 과 평행이동을 구한다. */
function superpose(mob, ref){
  const n = ref.length;
  const cm = [0,0,0], cr = [0,0,0];
  for(let i=0;i<n;i++) for(let k=0;k<3;k++){ cm[k]+=mob[i][k]/n; cr[k]+=ref[i][k]/n; }
  const P = mob.map(p=>[p[0]-cm[0],p[1]-cm[1],p[2]-cm[2]]);
  const Q = ref.map(p=>[p[0]-cr[0],p[1]-cr[1],p[2]-cr[2]]);

  // 공분산 행렬
  const M = [[0,0,0],[0,0,0],[0,0,0]];
  for(let i=0;i<n;i++) for(let r=0;r<3;r++) for(let c=0;c<3;c++) M[r][c]+=P[i][r]*Q[i][c];
  const [[Sxx,Sxy,Sxz],[Syx,Syy,Syz],[Szx,Szy,Szz]] = M;

  const N = [
    [Sxx+Syy+Szz, Syz-Szy,     Szx-Sxz,     Sxy-Syx],
    [Syz-Szy,     Sxx-Syy-Szz, Sxy+Syx,     Szx+Sxz],
    [Szx-Sxz,     Sxy+Syx,     -Sxx+Syy-Szz,Syz+Szy],
    [Sxy-Syx,     Szx+Sxz,     Syz+Szy,     -Sxx-Syy+Szz]
  ];
  const {values, vectors} = jacobi4(N);
  let bi = 0; for(let i=1;i<4;i++) if(values[i] > values[bi]) bi = i;
  const q = [vectors[0][bi], vectors[1][bi], vectors[2][bi], vectors[3][bi]];
  const nq = Math.hypot(...q); const [w,x,y,z] = q.map(v=>v/nq);

  const R = [
    [w*w+x*x-y*y-z*z, 2*(x*y-w*z),     2*(x*z+w*y)],
    [2*(x*y+w*z),     w*w-x*x+y*y-z*z, 2*(y*z-w*x)],
    [2*(x*z-w*y),     2*(y*z+w*x),     w*w-x*x-y*y+z*z]
  ];
  return {R, cm, cr};
}

const applyRT = (p, {R, cm, cr}) => {
  const d = [p[0]-cm[0], p[1]-cm[1], p[2]-cm[2]];
  return [0,1,2].map(r => R[r][0]*d[0]+R[r][1]*d[1]+R[r][2]*d[2] + cr[r]);
};
const dist2 = (a,b) => (a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2;

/* ── 두 구조의 편차 계산 ────────────────────────────────────
   정렬 구간: 양쪽 모두 pLDDT ≥ plddtMin 인 잔기
     (전장 입력에서 무질서 꼬리가 정렬을 오염시키는 것을 막는다.
      임의 단백질에 적용 가능한 일반 규칙이다.)
   국소 RMSD: 변이 잔기 8 Å 이내 Cα — 정렬을 유지한 채 측정 */
function comparePair(A, B, mutPos, cutoff=8.0, plddtMin=70){
  const common = Object.keys(A).map(Number).filter(r => B[r] !== undefined);
  if(!common.length) return null;

  const core = common.filter(r => A[r].plddt >= plddtMin && B[r].plddt >= plddtMin);
  const fit = core.length >= 20 ? core : common;      // 잘 접힌 부분이 너무 적으면 전체 사용

  const T = superpose(fit.map(r=>B[r].xyz), fit.map(r=>A[r].xyz));
  const Bt = {}; common.forEach(r => Bt[r] = applyRT(B[r].xyz, T));

  const rmsdOf = ids => Math.sqrt(ids.reduce((s,r)=>s+dist2(A[r].xyz, Bt[r]), 0) / ids.length);

  let local = null;
  if(A[mutPos]){
    const c = A[mutPos].xyz;
    const near = common.filter(r => dist2(A[r].xyz, c) <= cutoff*cutoff);
    if(near.length) local = rmsdOf(near);
  }
  return {global: rmsdOf(fit), local, nFit: fit.length, usedCore: core.length >= 20};
}

/* ── 그룹 간 검정 ───────────────────────────────────────────
   기준 반복 = 같은 그룹 내 쌍 · 변이 비교 = 그룹 간 쌍 */
function runDetectability(wtSet, mutSet, mutPos){
  const pairsIn = set => {
    const out = [];
    for(let i=0;i<set.length;i++) for(let j=i+1;j<set.length;j++) out.push([set[i],set[j]]);
    return out;
  };
  const val = pr => { const r = comparePair(pr[0], pr[1], mutPos); return r && r.local; };

  const baseline = [...pairsIn(wtSet), ...pairsIn(mutSet)].map(val).filter(v=>v!=null);
  const signal = [];
  wtSet.forEach(a => mutSet.forEach(b => { const v = val([a,b]); if(v!=null) signal.push(v); }));

  if(baseline.length < 3 || signal.length < 3) return {error:'비교 가능한 쌍이 부족합니다 (각 조건에 구조 3개 이상 필요)'};

  const {p} = mannWhitneyU(signal, baseline, 'greater');
  const meta = comparePair(wtSet[0], mutSet[0], mutPos);
  return {baseline, signal, p, detected: p < 0.05,
          nFit: meta?.nFit, usedCore: meta?.usedCore};
}

/* ── ZIP 파싱 ───────────────────────────────────────────────
   AlphaFold Server / Boltz 출력 구조를 폴더 또는 파일명으로 묶는다. */
async function readZip(file){
  const zip = await JSZip.loadAsync(file);
  const groups = {};                    // 그룹명 → [{name, struct}]
  const names = Object.keys(zip.files).filter(n => n.toLowerCase().endsWith('.cif') && !zip.files[n].dir);

  for(const n of names){
    const text = await zip.files[n].async('string');
    const st = parseCifCA(text);
    if(!Object.keys(st).length) continue;
    const parts = n.split('/').filter(Boolean);
    // 폴더가 있으면 폴더명, 없으면 model_N 을 뗀 파일명으로 묶는다
    const g = parts.length > 1
      ? parts[parts.length-2]
      : parts[0].replace(/_?model_\d+\.cif$/i,'').replace(/\.cif$/i,'');
    (groups[g] ||= []).push({name: parts[parts.length-1], struct: st});
  }
  return groups;
}
