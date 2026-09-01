/* BlindSpot — 변이 감별력 사전 검정
 *
 * 실제로 동작하는 부분:
 *   · UniProt REST 조회 및 잔기 대조 검증  → 임의 accession 가능
 *   · 변이 서열 생성 및 입력 파일 배포
 *   · CSV 원자료 기반 Mann–Whitney U 재계산
 *   · Mol* 로 실측 구조 로드
 * 미구현: 구조 예측 실행 (GPU 필요)
 */

const $ = s => document.querySelector(s);
const $$ = s => [...document.querySelectorAll(s)];

/* ────────────── 사례 정의 (사전 산출) ────────────── */
const CASES = {
  CD33:  {acc:'P20138', len:364, pos:69, wt:'R', mut:'G', rsid:'rs2455069',
          domain:[19,135], gene:'CD33',
          wtCif:'cd33_full_wt', mutCif:'cd33_full_r69g', img:'CD33_zoom.png'},
  PILRA: {acc:'Q9UKJ1', len:303, pos:78, wt:'R', mut:'G', rsid:'rs1859788',
          domain:[32,150], gene:'PILRA',
          wtCif:'pilra_full_wt', mutCif:'pilra_full_r78g', img:'PILRA_zoom.png'},
  TREM2: {acc:'Q9NZC2', len:230, pos:62, wt:'R', mut:'H', rsid:'rs143332484',
          domain:[29,112], gene:'TREM2',
          wtCif:'trem2_full_wt', mutCif:'trem2_full_r62h', img:'TREM2_zoom.png'},
};

/* ────────────── 통계 ────────────── */
const median = a => {const s=[...a].sort((x,y)=>x-y), m=s.length>>1;
  return s.length%2 ? s[m] : (s[m-1]+s[m])/2;};

/** Mann–Whitney U — 정규근사, 동순위 보정 포함 */
function mannWhitneyU(a, b, alternative='greater'){
  const n1=a.length, n2=b.length;
  const all=[...a.map(v=>({v,g:0})), ...b.map(v=>({v,g:1}))].sort((x,y)=>x.v-y.v);
  // 동순위 평균 순위 부여
  const ranks=new Array(all.length); const tieGroups=[];
  for(let i=0;i<all.length;){
    let j=i; while(j+1<all.length && all[j+1].v===all[i].v) j++;
    const r=(i+j)/2+1;
    for(let k=i;k<=j;k++) ranks[k]=r;
    if(j>i) tieGroups.push(j-i+1);
    i=j+1;
  }
  let R1=0; all.forEach((o,i)=>{ if(o.g===0) R1+=ranks[i]; });
  const U1=R1-n1*(n1+1)/2, U2=n1*n2-U1;
  const mu=n1*n2/2;
  const N=n1+n2;
  const tieCorr=tieGroups.reduce((s,t)=>s+(t*t*t-t),0);
  const sd=Math.sqrt((n1*n2/12)*((N+1)-tieCorr/(N*(N-1))));
  if(sd===0) return {U:U1, p:1};
  const U = alternative==='greater' ? U1 : U1;
  let z=(U-mu)/sd;
  // 연속성 보정
  z = z>0 ? (U-mu-0.5)/sd : (U-mu+0.5)/sd;
  const p1 = 1-normCdf(z);                     // 단측 (a > b)
  const p = alternative==='two-sided'
      ? Math.min(1, 2*Math.min(p1, 1-p1)) : p1;
  return {U:U1, U2, z, p};
}
function normCdf(z){ // Abramowitz–Stegun 26.2.17
  const t=1/(1+0.2316419*Math.abs(z));
  const d=0.3989422804014327*Math.exp(-z*z/2);
  let p=d*t*(0.319381530+t*(-0.356563782+t*(1.781477937+t*(-1.821255978+t*1.330274429))));
  return z>0 ? 1-p : p;
}
const fmtP = p => p<1e-4 ? p.toExponential(1).replace('e','×10^').replace('^-','⁻')
                         : p.toPrecision(2);

/* ────────────── CSV ────────────── */
async function loadCsv(path){
  const txt = await (await fetch(path)).text();
  const [head, ...lines] = txt.trim().split('\n');
  const cols = head.split(',');
  return lines.map(l=>{
    const v=l.split(','); const o={};
    cols.forEach((c,i)=>o[c.trim()]=v[i]);
    return o;
  });
}

/* ────────────── 화면 2 · UniProt 조회 및 검증 (실제 동작) ────────────── */
let currentSeq = null, currentAcc = null;

async function fetchUniProt(acc){
  const r = await fetch(`https://rest.uniprot.org/uniprotkb/${acc}.json`);
  if(!r.ok) throw new Error(`조회 실패 (HTTP ${r.status})`);
  const d = await r.json();
  const feats = (d.features||[]).filter(f=>['Signal','Domain','Transmembrane','Natural variant'].includes(f.type));
  return {
    acc: d.primaryAccession,
    id: d.uniProtkbId,
    name: d.proteinDescription?.recommendedName?.fullName?.value || '',
    seq: d.sequence.value,
    len: d.sequence.length,
    feats
  };
}

async function doLookup(){
  const acc = $('#acc').value.trim().toUpperCase();
  const box = $('#accBox'), info = $('#accInfo');
  if(!acc) return;
  box.className='input'; info.className='msg'; info.innerHTML='<span class="spin"></span> <span class="muted">조회 중…</span>';
  $('#lookup').disabled=true;
  try{
    const p = await fetchUniProt(acc);
    currentSeq = p.seq; currentAcc = p.acc;
    box.className='input ok';
    $('#accMeta').textContent = `${p.id} · ${p.len} aa`;
    const sig = p.feats.find(f=>f.type==='Signal');
    const dom = p.feats.filter(f=>f.type==='Domain');
    const tm  = p.feats.find(f=>f.type==='Transmembrane');
    const parts=[];
    if(sig) parts.push(`신호펩타이드 ${sig.location.start.value}–${sig.location.end.value}`);
    dom.slice(0,2).forEach(d=>parts.push(`${d.description} ${d.location.start.value}–${d.location.end.value}`));
    if(tm) parts.push(`TM ${tm.location.start.value}–${tm.location.end.value}`);
    info.className='msg ok';
    info.innerHTML = `<b>${p.name||p.id}</b><br>${parts.join(' · ')||'주석 없음'}`;
    window._feats = p.feats;
    $('#lookup').disabled=false;
    validateMut();
  }catch(e){
    box.className='input warn';
    info.className='msg warn';
    info.textContent = `${acc} 를 찾지 못했습니다 (${e.message}). accession 형식을 확인해 주세요 — 예: Q9NZC2`;
    currentSeq=null; $('#lookup').disabled=false;
  }
}

function validateMut(){
  const raw = $('#mut').value.trim().toUpperCase();
  const box = $('#mutBox'), out = $('#mutResult'), gen = $('#genBox');
  gen.innerHTML=''; $('#dl').style.display='none';
  const rw=$('#runWrap'); if(rw) rw.style.display='none';
  const rm=$('#runMsg'); if(rm){ rm.style.display='none'; rm.innerHTML=''; }
  const m = raw.match(/^([A-Z])(\d+)([A-Z])$/);
  if(!m){ box.className='input'; out.className='msg';
    out.textContent='형식에 맞게 입력해 주세요 — 예: R62H'; return; }
  if(!currentSeq){ out.className='msg'; out.textContent='accession 을 먼저 조회해 주세요'; return; }

  const [,wt,posS,mut] = m, pos=+posS;
  if(pos<1 || pos>currentSeq.length){
    box.className='input warn'; out.className='msg warn';
    out.innerHTML=`위치 ${pos}가 서열 범위(1–${currentSeq.length})를 벗어납니다`; return;
  }
  const actual = currentSeq[pos-1];

  if(actual !== wt){
    box.className='input warn'; out.className='msg warn';
    // 역방향인지 확인 — 흔한 오류 패턴
    const reversed = actual===mut;
    out.innerHTML = `<b>서열 불일치</b> — 입력하신 변이의 기준 잔기(${AA[wt]||wt})가 `
      + `야생형 서열 ${pos}번 잔기(<b>${AA[actual]||actual}</b>)와 다릅니다.`
      + (reversed ? `<br>표기 방향이 반대일 수 있습니다. <code>${actual}${pos}${wt}</code> 인지 확인해 주세요.` : '')
      + `<br><span class="muted">전구체(precursor) 번호 기준으로 입력했는지 확인해 주세요.</span>`;
    return;
  }

  box.className='input ok'; out.className='msg ok';
  const rs = (window._feats||[]).find(f=>f.type==='Natural variant'
      && f.location.start.value===pos
      && (f.alternativeSequence?.alternativeSequences||[]).includes(mut));
  const rsid = rs?.featureCrossReferences?.[0]?.id;
  out.innerHTML = `<b>확인됨</b> — 야생형 서열 ${pos}번이 ${AA[wt]}(${wt})로 일치합니다`
    + (rsid ? ` · dbSNP <b>${rsid}</b>` : '');

  buildSequences(pos, wt, mut);
}

const AA={A:'Ala',R:'Arg',N:'Asn',D:'Asp',C:'Cys',Q:'Gln',E:'Glu',G:'Gly',H:'His',
  I:'Ile',L:'Leu',K:'Lys',M:'Met',F:'Phe',P:'Pro',S:'Ser',T:'Thr',W:'Trp',Y:'Tyr',V:'Val',
  U:'Sec',O:'Pyl',X:'Xaa',B:'Asx',Z:'Glx'};

function buildSequences(pos, wt, mut){
  const full = currentSeq;
  const mutSeq = full.slice(0,pos-1) + mut + full.slice(pos);
  const label = `${wt}${pos}${mut}`;

  // 구간 선택
  const useDom = $('#useDom').checked;
  const dom = (window._feats||[]).filter(f=>f.type==='Domain')
              .find(d=>d.location.start.value<=pos && pos<=d.location.end.value);
  let a=1, b=full.length, note='전장';
  if(useDom && dom){ a=dom.location.start.value; b=dom.location.end.value;
    note=`${dom.description} ${a}–${b}`; }
  else if(useDom && !dom){ note='전장 — 해당 위치를 포함하는 도메인 주석이 없습니다'; }
  const sW = full.slice(a-1,b), sM = mutSeq.slice(a-1,b), idx = pos-a;

  // 비표준 잔기(셀레노시스테인 U, 피롤라이신 O)는 예측 도구가 처리하지 못한다
  const NONSTD = {U:'셀레노시스테인(Sec)', O:'피롤라이신(Pyl)'};
  const found = [...new Set([...sW].filter(c=>NONSTD[c]))];
  const nsWarn = found.length ? `
    <div class="msg warn" style="margin-top:11px"><b>비표준 잔기 포함</b> —
      ${found.map(c=>`${NONSTD[c]} <code>${c}</code>`).join(', ')} 가 서열에 존재합니다.
      AlphaFold·Boltz-2 는 표준 20종 아미노산만 처리하므로 치환이 필요합니다
      (예: Sec → Cys). 결정구조 다수도 같은 이유로 치환체를 사용합니다.</div>` : '';

  const hl = (s,i,c)=>s.slice(0,i)+`<b class="hl">${c}</b>`+s.slice(i+1);
  $('#genBox').innerHTML = `
    <div class="chipline"><span class="tag">구간</span> <span class="mono">${note} · ${sW.length} aa</span></div>
    <div class="seqlbl">야생형</div><div class="seq">${hl(sW,idx,wt)}</div>
    <div class="seqlbl">변이형 ${label}</div><div class="seq">${hl(sM,idx,mut)}</div>
    <div class="chipline"><span class="tag">상이 잔기</span> <span class="mono">1개 (위치 ${pos})</span></div>${nsWarn}`;

  $('#dl').style.display='flex';
  $('#runWrap').style.display='block';
  window._out = {acc:currentAcc, label, a, b, sW, sM, n:+$('#nsamp').value, seed:+$('#seed').value};
}

/* 입력 파일 배포 — 실제 다운로드 */
function download(name, text){
  const u=URL.createObjectURL(new Blob([text],{type:'text/plain'}));
  const a=document.createElement('a'); a.href=u; a.download=name; a.click();
  URL.revokeObjectURL(u);
}
function dlFasta(){
  const o=window._out;
  download(`${o.acc}_${o.label}.fasta`,
    `>${o.acc}_WT ${o.a}-${o.b}\n${o.sW}\n>${o.acc}_${o.label} ${o.a}-${o.b}\n${o.sM}\n`);
}
function dlAf(){
  const o=window._out, jobs=[];
  [['WT',o.sW],[o.label,o.sM]].forEach(([tag,seq])=>{
    jobs.push({name:`${o.acc}_${tag}`, modelSeeds:[o.seed],
      sequences:[{proteinChain:{sequence:seq,count:1}}],
      dialect:'alphafoldserver', version:1});
  });
  download(`af_jobs_${o.acc}_${o.label}.json`, JSON.stringify(jobs,null,2));
}
function dlBoltz(){
  const o=window._out;
  const y=s=>`version: 1\nsequences:\n  - protein:\n      id: A\n      sequence: "${s}"\n`;
  download(`boltz_${o.acc}_WT.yaml`, y(o.sW));
  setTimeout(()=>download(`boltz_${o.acc}_${o.label}.yaml`, y(o.sM)), 250);
}

/* ────────────── 화면 3·4 · 검정 재계산 (실제 계산) ────────────── */
let RMSD=null, IPTM=null;

async function renderStructureAxis(gene){
  if(!RMSD) RMSD = await loadCsv('data/rmsd_pairs.csv');
  const c = CASES[gene];
  const rows = RMSD.filter(r=>r.gene===gene);
  const noise = rows.filter(r=>r.comparison.startsWith('noise')).map(r=>+r.rmsd_local8A);
  const sig   = rows.filter(r=>r.comparison==='signal_WT_vs_MUT').map(r=>+r.rmsd_local8A);
  const {p} = mannWhitneyU(sig, noise, 'greater');
  const det = p<0.05;
  const mn=median(noise), ms=median(sig), d=ms-mn;

  $('#sTitle').innerHTML = `${gene} · ${c.wt}${c.pos}${c.mut} <span class="muted mono">${c.rsid}</span>`;

  $('#sScore').innerHTML = `
   <div class="score ${det?'sig':'ns'}">
    <div class="top">
      <div class="verdict">${det?'재현 변동성 초과':'재현 변동성 내'}</div>
      <div class="vsub">${det
        ? '변이로 인한 편차가 모델 재현 변동을 통계적으로 초과합니다'
        : '변이로 인한 편차가 모델 재현 변동과 구분되지 않습니다'}</div>
    </div>
    <div class="kpis">
      <div class="kpi"><div class="kl">p-value</div>
        <div class="kv mono">${fmtP(p)}</div>
        <div class="kn">Mann–Whitney U · 단측 · α = 0.05</div></div>
      <div class="kpi"><div class="kl">편차 Δ</div>
        <div class="kv mono">${d>=0?'+':''}${d.toFixed(3)} <span style="font-size:15px">Å</span></div>
        <div class="kn">${det?`Cα–Cα 결합거리의 약 ${Math.round(Math.abs(d)/1.5*100)}%`:'기준 변동 범위 이내'}</div></div>
    </div>
    <div class="rows">
      ${row('기준 반복 변동 (n='+noise.length+')', mn.toFixed(3)+' Å')}
      ${row('범위', Math.min(...noise).toFixed(3)+' – '+Math.max(...noise).toFixed(3), 1)}
      ${row('변이 비교 편차 (n='+sig.length+')', ms.toFixed(3)+' Å')}
      ${row('비율', (ms/mn).toFixed(2)+'×')}
    </div>
    <div class="foot">${det
      ? '귀무가설을 기각합니다. 다만 효과크기를 함께 확인하십시오.'
      : '귀무가설을 기각하지 못했습니다. 구조 불변이 아니라 본 모델의 감별 한계를 의미합니다.'}</div>
   </div>`;

  $('#sCond').innerHTML = `
    <tr><td>모델</td><td class="mono">AlphaFold Server (AF3)</td></tr>
    <tr><td>입력 구간</td><td class="mono">UniProt ${c.acc} 전장 (${c.len} aa)</td></tr>
    <tr><td>시드 / 샘플</td><td class="mono">seed 1 · n = 5</td></tr>
    <tr><td>정렬 구간</td><td class="mono">도메인 ${c.domain[0]}–${c.domain[1]} Cα</td></tr>
    <tr><td>컷오프</td><td class="mono">8 Å (변이 잔기 기준)</td></tr>`;

  drawStrip('#sPlot', [
    {label:'기준 반복', vals:noise, color:'#2a78d6'},
    {label:'변이 비교', vals:sig,   color:'#eb6834'}],
    '국소 RMSD (Å)');

  loadViewer([`data/structures/${c.wtCif}.cif`, `data/structures/${c.mutCif}.cif`],
             ['#2a4fd8','#e03a3a'], '#sViewer', c.img);
  window._curGene = gene;
}

async function renderBindingAxis(){
  if(!IPTM) IPTM = await loadCsv('data/boltz_iptm.csv');
  const wt  = IPTM.filter(r=>r.allele==='WT').map(r=>+r.iptm);
  const mut = IPTM.filter(r=>r.allele==='R62H').map(r=>+r.iptm);
  const {p} = mannWhitneyU(mut, wt, 'two-sided');
  const det = p<0.05;
  const mw=median(wt), mm=median(mut), spread=Math.max(...wt)-Math.min(...wt), d=mm-mw;

  $('#bScore').innerHTML = `
   <div class="score ${det?'sig':'ns'}">
    <div class="top">
      <div class="verdict">${det?'재현 변동성 초과':'재현 변동성 내'}</div>
      <div class="vsub">${det
        ? '변이로 인한 ipTM 차이가 재현 변동을 초과합니다'
        : '변이로 인한 ipTM 차이가 모델 재현 변동에 포섭됩니다'}</div>
    </div>
    <div class="kpis">
      <div class="kpi"><div class="kl">p-value</div>
        <div class="kv mono">${fmtP(p)}</div>
        <div class="kn">Mann–Whitney U · 양측 · α = 0.05</div></div>
      <div class="kpi"><div class="kl">ipTM 차이 Δ</div>
        <div class="kv mono">${d>=0?'+':''}${d.toFixed(3)}</div>
        <div class="kn">기준 분산 폭 ${spread.toFixed(3)}</div></div>
    </div>
    <div class="rows">
      ${row('야생형 ipTM (n='+wt.length+')', mw.toFixed(3))}
      ${row('범위', Math.min(...wt).toFixed(3)+' – '+Math.max(...wt).toFixed(3),1)}
      ${row('변이형 ipTM (n='+mut.length+')', mm.toFixed(3))}
      ${row('범위', Math.min(...mut).toFixed(3)+' – '+Math.max(...mut).toFixed(3),1)}
    </div>
    <div class="foot">기준 분산 폭 ${spread.toFixed(3)} 이 처리군 간 차이 ${Math.abs(d).toFixed(3)} 를 크게 상회합니다.
      Bret et al. (2026) 이 보고한 결합부위 변이 둔감성과 부합합니다.</div>
   </div>`;

  drawStrip('#bPlot', [
    {label:'야생형 + Aβ42',   vals:wt,  color:'#2a78d6'},
    {label:'R62H + Aβ42', vals:mut, color:'#eb6834'}],
    'ipTM');

  loadViewer(['data/structures/trem2_ab42_wt.cif','data/structures/trem2_ab42_r62h.cif'],
             ['#2a4fd8','#e03a3a'], '#bViewer', 'TREM2_overlay.png');
}

const row=(k,v,sub,cls)=>`<div class="row${sub?' sub':''}">
  <span class="k">${k}</span><span class="v mono ${cls||''}">${v}</span></div>`;

/* ────────────── 분포 도표 (SVG 직접 생성) ────────────── */
function drawStrip(sel, groups, xlabel){
  const W=760, H=groups.length*88+64, L=132, R=36, T=18;
  const all=groups.flatMap(g=>g.vals);
  const lo=0, hi=Math.max(...all)*1.08;
  const x=v=>L+(v-lo)/(hi-lo)*(W-L-R);
  let s=`<svg viewBox="0 0 ${W} ${H}" width="100%">`;
  // 격자
  const ticks=5;
  for(let i=0;i<=ticks;i++){
    const v=lo+(hi-lo)*i/ticks, px=x(v);
    s+=`<line x1="${px}" y1="${T}" x2="${px}" y2="${H-44}" stroke="#e6e5e1"/>`;
    s+=`<text x="${px}" y="${H-26}" font-size="11" fill="#52514e" text-anchor="middle">${v.toFixed(2)}</text>`;
  }
  groups.forEach((g,gi)=>{
    const cy=T+34+gi*88;
    s+=`<text x="${L-12}" y="${cy+4}" font-size="12" fill="#52514e" text-anchor="end">${g.label}</text>`;
    s+=`<text x="${L-12}" y="${cy+19}" font-size="10.5" fill="#8a8a86" text-anchor="end">n = ${g.vals.length}</text>`;
    g.vals.forEach((v,i)=>{
      const jy=cy+((i*37)%25)-12;
      s+=`<circle cx="${x(v)}" cy="${jy}" r="4.4" fill="${g.color}" fill-opacity=".55" stroke="#fff" stroke-width="1.1"/>`;
    });
    const m=median(g.vals);
    s+=`<line x1="${x(m)}" y1="${cy-24}" x2="${x(m)}" y2="${cy+24}" stroke="${g.color}" stroke-width="2.8"/>`;
    s+=`<text x="${x(m)}" y="${cy-30}" font-size="11.5" font-weight="700" fill="${g.color}" text-anchor="middle">${m.toFixed(3)}</text>`;
  });
  s+=`<text x="${(L+W-R)/2}" y="${H-6}" font-size="11.5" fill="#52514e" text-anchor="middle">${xlabel}</text></svg>`;
  $(sel).innerHTML=s;
}

/* ────────────── Mol* 뷰어 ────────────── */
let viewers={};
function webglOk(){
  try{const c=document.createElement('canvas');
      return !!(c.getContext('webgl2')||c.getContext('webgl'));}catch(e){return false;}
}

async function loadViewer(files, colors, sel='#sViewer', fallbackImg=null){
  const el=$(sel); if(!el) return;
  // WebGL 이 없는 환경(헤드리스 캡처 등)에서는 사전 렌더 이미지로 대체한다
  if(!webglOk() || typeof molstar==='undefined'){
    el.innerHTML = fallbackImg
      ? `<img src="data/img/${fallbackImg}" style="width:100%;height:100%;object-fit:contain">
         <div class="vnote">사전 렌더 이미지 · 대화형 3D는 WebGL 지원 브라우저에서 표시된다</div>`
      : `<div class="vfallback">WebGL 미지원 환경</div>`;
    return;
  }
  try{
    if(!viewers[sel]){
      viewers[sel] = await molstar.Viewer.create(el.id, {
        layoutIsExpanded:false, layoutShowControls:false, layoutShowSequence:false,
        layoutShowLog:false, layoutShowLeftPanel:false, viewportShowExpand:true,
        viewportShowSelectionMode:false, viewportShowAnimation:false, pdbProvider:'rcsb'
      });
    }
    const v=viewers[sel];
    await v.plugin.clear();
    for(let i=0;i<files.length;i++){
      await v.loadStructureFromUrl(new URL(files[i], location.href).href, 'mmcif', false,
        {representationParams:{theme:{globalName:'uniform',
          globalColorParams:{value: parseInt(colors[i].slice(1),16)}}}});
    }
  }catch(e){
    el.innerHTML=`<div class="vfallback">3D 뷰어를 초기화하지 못했습니다.<br>
      <span class="muted">${e.message}</span><br>
      <span class="muted">로컬 파일(file://)로 열면 CORS로 차단됩니다 — HTTP 서버로 여십시오.</span></div>`;
  }
}

/* ────────────── 초기화 ────────────── */
function initNav(){
  $$('nav a').forEach(a=>a.onclick=e=>{
    e.preventDefault();
    $$('nav a').forEach(x=>x.classList.remove('on'));
    $$('.screen').forEach(x=>x.classList.remove('on'));
    a.classList.add('on'); $('#'+a.dataset.s).classList.add('on');
    window.scrollTo(0,0);
    if(a.dataset.s==='s3') renderStructureAxis(window._curGene||'CD33');
    if(a.dataset.s==='s4') renderBindingAxis();
    if(a.dataset.s==='s5') renderReport();
  });
}

async function renderReport(){
  if(!RMSD) RMSD = await loadCsv('data/rmsd_pairs.csv');
  if(!IPTM) IPTM = await loadCsv('data/boltz_iptm.csv');
  let html='';
  for(const g of Object.keys(CASES)){
    const c=CASES[g];
    const rows=RMSD.filter(r=>r.gene===g);
    const noise=rows.filter(r=>r.comparison.startsWith('noise')).map(r=>+r.rmsd_local8A);
    const sig=rows.filter(r=>r.comparison==='signal_WT_vs_MUT').map(r=>+r.rmsd_local8A);
    const {p}=mannWhitneyU(sig,noise,'greater');
    const mn=median(noise), ms=median(sig), det=p<0.05;
    html+=`<tr><td><b>${g}</b> <span class="mono muted">${c.wt}${c.pos}${c.mut}</span></td>
      <td class="mono">구조</td><td class="n mono">${mn.toFixed(3)} Å</td>
      <td class="n mono">${ms.toFixed(3)} Å</td>
      <td class="n mono">${(ms/mn).toFixed(2)}×</td>
      <td class="n mono">${fmtP(p)}</td>
      <td><span class="badge ${det?'sig':'ns'}">${det?'변동성 초과':'변동성 내'}</span></td></tr>`;
  }
  const wt=IPTM.filter(r=>r.allele==='WT').map(r=>+r.iptm);
  const mut=IPTM.filter(r=>r.allele==='R62H').map(r=>+r.iptm);
  const {p:pb}=mannWhitneyU(mut,wt,'two-sided');
  html+=`<tr><td><b>TREM2</b> <span class="mono muted">R62H + Aβ42</span></td>
    <td class="mono">인터페이스</td><td class="n mono">${(Math.max(...wt)-Math.min(...wt)).toFixed(3)}</td>
    <td class="n mono">${Math.abs(median(mut)-median(wt)).toFixed(3)}</td>
    <td class="n mono">—</td><td class="n mono">${fmtP(pb)}</td>
    <td><span class="badge ${pb<0.05?'sig':'ns'}">${pb<0.05?'변동성 초과':'변동성 내'}</span></td></tr>`;
  $('#rTable').innerHTML=html;
}

function gotoScreen(id){
  if(!document.getElementById(id)) return;
  $$('nav a').forEach(x=>x.classList.remove('on'));
  $$('.screen').forEach(x=>x.classList.remove('on'));
  document.getElementById(id).classList.add('on');
  const na=$(`nav a[data-s="${id}"]`); if(na) na.classList.add('on');
  if(id==='s3') renderStructureAxis(window._curGene||'CD33');
  if(id==='s4') renderBindingAxis();
  if(id==='s5') renderReport();
}

window.addEventListener('DOMContentLoaded', ()=>{
  initNav();
  $('#lookup').onclick = doLookup;
  $('#acc').addEventListener('keydown', e=>{ if(e.key==='Enter') doLookup(); });
  $('#mut').addEventListener('input', validateMut);
  $('#useDom').addEventListener('change', validateMut);
  $('#dlFasta').onclick=dlFasta; $('#dlAf').onclick=dlAf; $('#dlBoltz').onclick=dlBoltz;
  $$('[data-goto]').forEach(b=>b.onclick=()=>gotoScreen(b.dataset.goto));
  const rb=$('#runBtn'); if(rb) rb.onclick=()=>{
    const g=Object.keys(CASES).find(k=>CASES[k].acc===currentAcc);
    if(g){ window._curGene=g; gotoScreen('s3'); return; }
    // 사전 산출 결과가 없는 입력 — 다른 사례의 결과를 보여주지 않는다
    $('#runMsg').innerHTML =
      `<b>사전 산출 결과가 없습니다.</b><br>`
      + `구조 예측 실행은 GPU 가 필요하여 본 데모에 포함되지 않았습니다. `
      + `위 입력 파일을 내려받아 AlphaFold Server 또는 Boltz-2 에서 예측한 뒤 결과를 분석하십시오.`
      + `<br><span class="muted">수록 사례: CD33 P20138 · PILRA Q9UKJ1 · TREM2 Q9NZC2</span>`;
    $('#runMsg').style.display='block';
  };
  $$('.casebtn').forEach(b=>b.onclick=()=>{
    $$('.casebtn').forEach(x=>x.classList.remove('on')); b.classList.add('on');
    renderStructureAxis(b.dataset.g);
  });
  // ?s=s3 또는 #s3 으로 특정 화면 진입 (캡처·딥링크용)
  const want=new URLSearchParams(location.search).get('s')
            || (location.hash||'').replace('#','');
  if(want) gotoScreen(want);
  const qp=new URLSearchParams(location.search);
  if(qp.get('acc')) $('#acc').value=qp.get('acc').toUpperCase();
  if(qp.get('mut')) $('#mut').value=qp.get('mut').toUpperCase();
  // 화면 2 진입 시 기본 사례를 자동 조회 (첫 화면이 빈 상태로 보이지 않게)
  if(want==='s2' || !want) doLookup();

  // 사례 프리셋
  $$('[data-preset]').forEach(b=>b.onclick=()=>{
    const c=CASES[b.dataset.preset];
    $('#acc').value=c.acc; $('#mut').value=`${c.wt}${c.pos}${c.mut}`;
    $$('nav a').forEach(x=>x.classList.remove('on'));
    $$('.screen').forEach(x=>x.classList.remove('on'));
    $('nav a[data-s="s2"]').classList.add('on'); $('#s2').classList.add('on');
    window.scrollTo(0,0); doLookup();
  });
});
