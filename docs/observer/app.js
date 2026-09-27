(() => {
  'use strict';

  const DATA = {
    index: './data/observer-index.json',
    earth: './data/earth-ledger.json',
    space: './data/space-ledger.json',
    bridge: './data/bridge-ledger.json',
    plates: './data/plate_manifest.csv'
  };

  const fallbackEarth = [
    {site:'Pantex Plant',lat:35.321,lng:-101.563,date:'2015-09-01',row_kind:'EVENT',source_grade:'P0_PRIMARY_OFFICIAL_INCIDENT_RECORD',classification_state:'UNRESOLVED_SOURCE_IDENTITY',source_family:'JANUS_COSMOS'},
    {site:'Pilgrim Nuclear Power Station',lat:41.944,lng:-70.579,date:'2015-10-01',row_kind:'EVENT',source_grade:'P0_NARA_SERIES_BOUND',classification_state:'UNRESOLVED_OBSERVATION',source_family:'JANUS_COSMOS'},
    {site:'Los Alamos',lat:35.881,lng:-106.299,date:'1949',row_kind:'HISTORICAL_ARCHIVAL_CASE_FAMILY',source_grade:'P0_OFFICIAL_ARCHIVAL_RELEASE',classification_state:'INSTITUTIONAL_INVESTIGATION',source_family:'JANUS_COSMOS'}
  ];

  const fallbackPlateCSV = `plate_id,ra_deg,dec_deg
XE002,2.348777,84.754856
XE003,44.912000,84.549934
XE004,83.899465,84.044706
XE005,120.975302,83.586154
XE006,154.842758,83.263035
XE007,188.661150,83.255440
XE008,221.640565,83.384269`;

  const continents = [
    [[72,-165],[62,-150],[53,-135],[48,-125],[35,-120],[23,-108],[19,-98],[25,-82],[42,-70],[52,-58],[61,-74],[70,-105],[72,-165]],
    [[12,-80],[4,-76],[-12,-72],[-30,-70],[-55,-68],[-43,-55],[-22,-45],[0,-50],[10,-65],[12,-80]],
    [[70,-10],[61,17],[49,30],[38,24],[35,8],[35,-10],[19,-17],[2,-8],[-18,12],[-35,20],[-34,35],[-12,45],[11,50],[31,35],[43,25],[52,10],[70,-10]],
    [[70,28],[72,80],[62,128],[52,155],[36,142],[22,112],[8,80],[24,57],[39,45],[55,48],[70,28]],
    [[-10,112],[-21,114],[-35,116],[-42,145],[-25,154],[-11,145],[-10,112]],
    [[83,-45],[75,-18],[66,-35],[62,-52],[70,-65],[83,-45]]
  ];

  const state = {index:null,earth:null,space:null,bridge:null,earthPoints:[],platePoints:[],skyTargets:[],spacePoints:[],earthSphere:null,spaceSphere:null,mode:'earth'};
  const $ = id => document.getElementById(id);
  const rad = d => d * Math.PI / 180;
  const clamp = (v,a,b) => Math.max(a,Math.min(b,v));
  const safe = (v,f='—') => v === undefined || v === null || v === '' ? f : String(v);

  function tickClock(){ $('clock').textContent = new Date().toISOString().replace('T',' ').slice(0,19)+' UTC'; }
  tickClock(); setInterval(tickClock,1000);

  function gradeShort(v){
    if(!v) return 'OPEN';
    if(String(v).startsWith('P0_')) return 'P0';
    if(String(v).includes('PRIMARY')) return 'PRIMARY';
    if(String(v).includes('META')) return 'META';
    return String(v).slice(0,12);
  }

  function earthColor(p){
    const s=(p.classification_state||'').toUpperCase();
    if(p.source_family==='JANUS_META') return '#48e59a';
    if(s.includes('POLICY')) return '#ffbd5a';
    if(s.includes('UAS')) return '#9a8cff';
    if(s.includes('UNRESOLVED')) return '#53d7ff';
    return '#65b7ff';
  }

  function decorateEarth(p){ return {...p,color:earthColor(p),size:(p.row_kind==='EVENT'?6.5:(p.source_family==='JANUS_META'?4.3:4.8))}; }

  function skyPoint(p, layer){
    const ra=Number(p.ra), dec=Number(p.dec);
    return {...p,ra,dec,lat:dec,lng:(((-ra + 180) % 360) + 360) % 360 - 180,layer,color:layer==='janus' ? '#ffbd5a' : '#5a9dff',size:layer==='janus' ? 5.6 : 1.65};
  }

  async function fetchJSON(url){ const r=await fetch(url,{cache:'no-store'}); if(!r.ok) throw new Error(`${r.status} ${url}`); return r.json(); }
  async function fetchText(url){ const r=await fetch(url,{cache:'no-store'}); if(!r.ok) throw new Error(`${r.status} ${url}`); return r.text(); }

  function parsePlateCSV(csv){
    const lines=csv.trim().split(/\r?\n/); if(lines.length<2) return [];
    const h=lines.shift().split(',').map(x=>x.trim()); const i=Object.fromEntries(h.map((x,n)=>[x,n]));
    return lines.map(line=>{ const c=line.split(','),ra=Number(c[i.ra_deg]),dec=Number(c[i.dec_deg]); if(!c[i.plate_id]||!Number.isFinite(ra)||!Number.isFinite(dec)) return null; return skyPoint({plate_id:c[i.plate_id],label:c[i.plate_id],ra,dec,row_kind:'POSS1_PLATE_CENTER',classification_state:'SURVEY_POSITION_NOT_ANOMALY',source_family:'POSS1_PUBLIC_MANIFEST'},'plate'); }).filter(Boolean);
  }

  function seededStars(count){ let seed=0x7f4a7c15; const rnd=()=>{seed=(1664525*seed+1013904223)>>>0;return seed/4294967296;}; return Array.from({length:count},()=>({x:rnd(),y:rnd(),a:.12+rnd()*.65,s:.35+rnd()*1.25})); }

  class Sphere {
    constructor(node,mode,onHover,onSelect){
      this.node=node; this.mode=mode; this.points=[]; this.screenPoints=[]; this.yaw=mode==='earth'?rad(-38):rad(15); this.pitch=mode==='earth'?rad(-8):rad(12); this.zoom=1; this.dragging=false; this.lastX=0; this.lastY=0; this.moved=0; this.lastInteraction=0; this.hovered=null; this.selected=null; this.onHover=onHover; this.onSelect=onSelect; this.bgStars=seededStars(620);
      this.canvas=document.createElement('canvas'); Object.assign(this.canvas.style,{position:'absolute',inset:'0',width:'100%',height:'100%',display:'block',touchAction:'none'}); node.replaceChildren(this.canvas); this.ctx=this.canvas.getContext('2d',{alpha:true}); this.bind(); this.resize(); this.ro=new ResizeObserver(()=>this.resize()); this.ro.observe(node); this.loop=this.loop.bind(this); requestAnimationFrame(this.loop);
    }
    setPoints(points){this.points=points||[];return this;}
    resize(){ const r=this.node.getBoundingClientRect(),dpr=Math.min(devicePixelRatio||1,2); this.w=Math.max(1,r.width);this.h=Math.max(1,r.height);this.canvas.width=Math.floor(this.w*dpr);this.canvas.height=Math.floor(this.h*dpr);this.ctx.setTransform(dpr,0,0,dpr,0,0);this.radius=Math.min(this.w,this.h)*(this.mode==='earth'?.37:.405)*this.zoom;this.cx=this.w*.5;this.cy=this.h*.54; }
    bind(){
      this.canvas.addEventListener('pointerdown',e=>{this.dragging=true;this.lastX=e.clientX;this.lastY=e.clientY;this.moved=0;this.lastInteraction=performance.now();this.canvas.setPointerCapture?.(e.pointerId);});
      this.canvas.addEventListener('pointermove',e=>{ if(this.dragging){const dx=e.clientX-this.lastX,dy=e.clientY-this.lastY;this.lastX=e.clientX;this.lastY=e.clientY;this.moved+=Math.hypot(dx,dy);this.yaw+=dx*.007;this.pitch=clamp(this.pitch+dy*.006,-1.28,1.28);this.lastInteraction=performance.now();} else this.pick(e.offsetX,e.offsetY,false); });
      this.canvas.addEventListener('pointerup',e=>{this.dragging=false;this.lastInteraction=performance.now();if(this.moved<8)this.pick(e.offsetX,e.offsetY,true);});
      this.canvas.addEventListener('pointercancel',()=>this.dragging=false);
      this.canvas.addEventListener('wheel',e=>{e.preventDefault();this.zoom=clamp(this.zoom*(e.deltaY>0?.92:1.08),.7,1.48);this.lastInteraction=performance.now();this.resize();},{passive:false});
    }
    xyz(lat,lng){ const la=rad(lat),lo=rad(lng);let x=Math.cos(la)*Math.sin(lo),y=Math.sin(la),z=Math.cos(la)*Math.cos(lo); const cy=Math.cos(this.yaw),sy=Math.sin(this.yaw);[x,z]=[x*cy+z*sy,-x*sy+z*cy]; const cp=Math.cos(this.pitch),sp=Math.sin(this.pitch);[y,z]=[y*cp-z*sp,y*sp+z*cp];return{x,y,z}; }
    project(lat,lng,alt=0){const p=this.xyz(lat,lng),s=1+alt;return{x:this.cx+p.x*this.radius*s,y:this.cy-p.y*this.radius*s,z:p.z};}
    pick(x,y,select){ let best=null,bd=18; for(const s of this.screenPoints){const d=Math.hypot(x-s.x,y-s.y);if(d<bd){best=s.point;bd=d;}} if(select&&best){this.selected=best;this.onSelect?.(best);} else if(!select&&best!==this.hovered){this.hovered=best;if(best)this.onHover?.(best);} this.canvas.style.cursor=best?'pointer':(this.dragging?'grabbing':'grab'); }
    gridLat(lat,color,alpha){ const c=this.ctx;c.beginPath();let open=false;for(let lng=-180;lng<=180;lng+=3){const p=this.project(lat,lng);if(p.z>0){if(!open){c.moveTo(p.x,p.y);open=true;}else c.lineTo(p.x,p.y);}else open=false;}c.strokeStyle=color;c.globalAlpha=alpha;c.lineWidth=1;c.stroke();c.globalAlpha=1; }
    gridLng(lng,color,alpha){ const c=this.ctx;c.beginPath();let open=false;for(let lat=-88;lat<=88;lat+=2){const p=this.project(lat,lng);if(p.z>0){if(!open){c.moveTo(p.x,p.y);open=true;}else c.lineTo(p.x,p.y);}else open=false;}c.strokeStyle=color;c.globalAlpha=alpha;c.lineWidth=1;c.stroke();c.globalAlpha=1; }
    drawContinents(){ const c=this.ctx;c.save();c.strokeStyle='#60dfff';c.lineWidth=1.25;c.shadowColor='#1ca9ff';c.shadowBlur=7;c.globalAlpha=.58;for(const poly of continents){c.beginPath();let open=false;for(const [lat,lng] of poly){const p=this.project(lat,lng,.004);if(p.z>.02){if(!open){c.moveTo(p.x,p.y);open=true;}else c.lineTo(p.x,p.y);}else open=false;}c.stroke();}c.restore(); }
    drawBackgroundStars(){ const c=this.ctx;c.save();for(const s of this.bgStars){c.globalAlpha=s.a;c.fillStyle='#dcecff';c.beginPath();c.arc(s.x*this.w,s.y*this.h,s.s,0,Math.PI*2);c.fill();}c.restore(); }
    drawSpaceLabels(){ const c=this.ctx;c.save();c.fillStyle='rgba(135,179,255,.72)';c.font='9px system-ui,sans-serif';c.textAlign='center';[['0h',0],['6h',90],['12h',180],['18h',270]].forEach(([label,ra])=>{const lng=(((-ra+180)%360)+360)%360-180,p=this.project(0,lng,.06);if(p.z>.05)c.fillText(label,p.x,p.y-8);});c.textAlign='left';c.fillStyle='rgba(255,189,90,.8)';c.fillText('RA ←  ·  Dec ↕  ·  ANGULAR SHELL (NOT DISTANCE)',18,this.h-50);c.restore(); }
    draw(){
      if(!this.ctx||!this.w)return;const c=this.ctx;c.clearRect(0,0,this.w,this.h);if(this.mode==='space')this.drawBackgroundStars();const R=this.radius;const glow=c.createRadialGradient(this.cx-R*.28,this.cy-R*.28,R*.05,this.cx,this.cy,R*1.15);
      if(this.mode==='earth'){glow.addColorStop(0,'#1765a5');glow.addColorStop(.38,'#092b54');glow.addColorStop(.82,'#041126');glow.addColorStop(1,'#010611');}else{glow.addColorStop(0,'rgba(25,45,105,.18)');glow.addColorStop(.62,'rgba(8,15,48,.24)');glow.addColorStop(.9,'rgba(3,7,24,.55)');glow.addColorStop(1,'rgba(1,3,12,.82)');}
      c.save();c.shadowColor=this.mode==='earth'?'#159cff':'#725cff';c.shadowBlur=35;c.fillStyle=glow;c.beginPath();c.arc(this.cx,this.cy,R,0,Math.PI*2);c.fill();c.restore();c.save();c.beginPath();c.arc(this.cx,this.cy,R-1,0,Math.PI*2);c.clip();const grid=this.mode==='earth'?'#45bfff':'#6e82ff';for(let lat=-60;lat<=60;lat+=30)this.gridLat(lat,grid,lat===0?.32:.15);for(let lng=-150;lng<=180;lng+=30)this.gridLng(lng,grid,.13);if(this.mode==='earth')this.drawContinents();c.restore();c.save();c.strokeStyle=this.mode==='earth'?'rgba(83,215,255,.72)':'rgba(119,104,255,.72)';c.lineWidth=1.5;c.shadowColor=this.mode==='earth'?'#42d0ff':'#715cff';c.shadowBlur=18;c.beginPath();c.arc(this.cx,this.cy,R,0,Math.PI*2);c.stroke();c.restore();
      this.screenPoints=[];const visible=this.points.map(point=>({point,p:this.project(point.lat,point.lng,.024)})).filter(o=>o.p.z>0).sort((a,b)=>a.p.z-b.p.z);
      for(const {point,p} of visible){const depth=.55+p.z*.72,base=point.size||3,size=base*depth;c.save();c.globalAlpha=point.layer==='plate'?(.30+p.z*.30):(.62+p.z*.35);c.fillStyle=point.color||'#53d7ff';c.shadowColor=point.color||'#53d7ff';c.shadowBlur=point.layer==='plate'?5:15;c.beginPath();c.arc(p.x,p.y,size,0,Math.PI*2);c.fill();if(point.layer!=='plate'){c.globalAlpha=.7;c.strokeStyle='#fff';c.lineWidth=.7;c.beginPath();c.arc(p.x,p.y,size+3,0,Math.PI*2);c.stroke();}c.restore();this.screenPoints.push({point,x:p.x,y:p.y});if(this.mode==='earth'&&(point.row_kind==='EVENT'||point===this.selected)){c.save();c.fillStyle='#e8f8ff';c.font='600 10px system-ui,sans-serif';c.shadowColor='#00101f';c.shadowBlur=5;c.fillText(safe(point.site),p.x+9,p.y-6);c.fillStyle='#79bde4';c.font='8px system-ui,sans-serif';c.fillText(safe(point.date||point.row_kind),p.x+9,p.y+5);c.restore();}}
      if(this.mode==='space')this.drawSpaceLabels();if(this.hovered){const s=this.screenPoints.find(x=>x.point===this.hovered);if(s){const label=this.mode==='earth'?safe(this.hovered.site):`${safe(this.hovered.label||this.hovered.plate_id)} · RA ${Number(this.hovered.ra).toFixed(2)}° · Dec ${Number(this.hovered.dec).toFixed(2)}°`;c.save();c.font='600 10px system-ui,sans-serif';const tw=Math.min(c.measureText(label).width,330);c.fillStyle='rgba(3,10,24,.9)';c.fillRect(s.x+10,s.y-26,tw+16,22);c.strokeStyle='rgba(83,215,255,.35)';c.strokeRect(s.x+10,s.y-26,tw+16,22);c.fillStyle='#dff6ff';c.fillText(label,s.x+18,s.y-11,330);c.restore();}}
    }
    loop(t){if(!this.dragging&&t-this.lastInteraction>1800)this.yaw+=this.mode==='earth'?.00055:-.00034;this.draw();requestAnimationFrame(this.loop);}
    focus(point){if(!point)return;this.selected=point;this.yaw=-rad(point.lng);this.pitch=rad(point.lat)*.42;this.lastInteraction=performance.now();}
  }

  function showDetail(p,domain){
    if(!p)return;
    if(domain==='space'){$('detail-grade').textContent=p.layer==='janus'?'JANUS':'POSS-I';const isJanus=p.layer==='janus';$('detail').className='detail-content';$('detail').innerHTML=`<span class="sub">SPACE · ${isJanus?'JANUS SOURCE TARGET':'SURVEY POSITION'}</span><h3>${safe(p.label||p.plate_id)}</h3><div class="detail-table"><div><span>RIGHT ASC.</span><b>${Number(p.ra).toFixed(6)}°</b></div><div><span>DECLINATION</span><b>${Number(p.dec).toFixed(6)}°</b></div><div><span>LAYER</span><b>${isJanus?'explicit RA/Dec discovered by autosync':'POSS-I public plate center'}</b></div><div><span>SOURCE</span><b>${safe(p.source_path)}</b></div><div><span>CONTRACT</span><b>Angular location only; no physical distance or anomaly identity inferred.</b></div></div>${p.source_url?`<a class="source-link" target="_blank" rel="noreferrer" href="${p.source_url}">OPEN SOURCE ↗</a>`:''}`;return;}
    $('detail-grade').textContent=gradeShort(p.source_grade);$('detail').className='detail-content';$('detail').innerHTML=`<span class="sub">EARTH · ${safe(p.row_kind)}</span><h3>${safe(p.site)}</h3><div class="detail-table"><div><span>DATE</span><b>${safe(p.date)}</b></div><div><span>REGION</span><b>${safe(p.region)}</b></div><div><span>STATE</span><b>${safe(p.classification_state)}</b></div><div><span>COORDS</span><b>${Number(p.lat).toFixed(5)}, ${Number(p.lng).toFixed(5)}</b></div><div><span>DISCOVERY</span><b>${safe(p.coordinate_method)}</b></div><div><span>BRIDGE</span><b>${safe(p.palomar_match_state)}</b></div><div><span>SOURCE</span><b>${safe(p.source_family)} · ${safe(p.source_path)}</b></div></div>${p.source_url?`<a class="source-link" target="_blank" rel="noreferrer" href="${p.source_url}">OPEN SOURCE ↗</a>`:''}`;
  }

  function buildEventList(points,unmapped){
    const list=$('earth-list');list.classList.remove('skeleton-list');list.innerHTML='';const unique=[];const seen=new Set();for(const p of points){const k=`${p.site}|${Number(p.lat).toFixed(5)}|${Number(p.lng).toFixed(5)}|${p.date||''}`;if(seen.has(k))continue;seen.add(k);unique.push(p);}unique.slice(0,80).forEach(p=>{const el=document.createElement('div');el.className='event-card';el.innerHTML=`<div class="row"><b>${safe(p.site)}</b><small class="badge">${gradeShort(p.source_grade)}</small></div><div class="row kind"><span>${safe(p.date||p.row_kind)}</span><small>${safe(p.source_family)}</small></div>`;el.onclick=()=>{showDetail(p,'earth');state.earthSphere?.focus(p);};list.appendChild(el);});if(unique.length>80){const m=document.createElement('div');m.className='microcopy';m.textContent=`${unique.length-80} more mapped rows are on the globe; sidebar is capped for readability.`;list.appendChild(m);}if(unmapped){const m=document.createElement('div');m.className='microcopy';m.textContent=`${unmapped} source rows have a site/name but no legal coordinate yet. They are preserved as UNRESOLVED, not dropped.`;list.appendChild(m);}
  }

  function initSpheres(){state.earthSphere=new Sphere($('earth-globe'),'earth',p=>showDetail(p,'earth'),p=>showDetail(p,'earth')).setPoints(state.earthPoints);state.spaceSphere=new Sphere($('space-globe'),'space',p=>showDetail(p,'space'),p=>showDetail(p,'space')).setPoints(state.spacePoints);}

  function setMode(mode){
    state.mode=mode;document.querySelectorAll('.mode').forEach(b=>b.classList.toggle('active',b.dataset.mode===mode));$('earth-globe').classList.toggle('active',mode==='earth');$('space-globe').classList.toggle('active',mode==='space');$('bridge-view').classList.toggle('active',mode==='bridge');
    if(mode==='earth'){$('stage-eyebrow').textContent='AUTOSYNC TERRESTRIAL WITNESS LANE';$('stage-title').textContent='JANUS Earth Evidence Globe';$('stage-big').textContent=state.earthPoints.length;$('stage-small').textContent='source-coordinate rows discovered automatically';$('hud-source').textContent='AUTOSYNC INDEX';$('hud-note').textContent='mapped + unresolved counts preserve source completeness';}
    else if(mode==='space'){$('stage-eyebrow').textContent='CELESTIAL ANGULAR WITNESS LANE';$('stage-title').textContent='JANUS Celestial Angular Shell';$('stage-big').textContent=state.spacePoints.length;$('stage-small').textContent=`${state.platePoints.length} POSS-I centers + ${state.skyTargets.length} JANUS sky targets`;$('hud-source').textContent='RA / DEC SHELL';$('hud-note').textContent='angular coordinates only — no fake 3D distance';}
    else{$('stage-eyebrow').textContent='CONTROLLED CROSS-DOMAIN LANE';$('stage-title').textContent='Earth ↔ Space Bridge';$('stage-big').textContent=state.index?.bridge?.palomar_day_rows_materialized||'0/312';$('stage-small').textContent='Palomar day rows materialized';$('hud-source').textContent='JOIN CONTRACT';$('hud-note').textContent='exact-day join remains blocked until denominator rights exist';}
  }

  document.querySelectorAll('.mode').forEach(btn=>btn.addEventListener('click',()=>setMode(btn.dataset.mode)));

  async function load(){
    $('load-state').textContent='AUTOSYNC';const tasks=await Promise.allSettled([fetchJSON(DATA.index),fetchJSON(DATA.earth),fetchJSON(DATA.space),fetchJSON(DATA.bridge),fetchText(DATA.plates)]);state.index=tasks[0].status==='fulfilled'?tasks[0].value:null;state.earth=tasks[1].status==='fulfilled'?tasks[1].value:null;state.space=tasks[2].status==='fulfilled'?tasks[2].value:null;state.bridge=tasks[3].status==='fulfilled'?tasks[3].value:null;
    if(state.index){state.earthPoints=(state.index.earth_points||[]).map(decorateEarth);state.platePoints=(state.index.poss1_plate_centers||[]).map(p=>skyPoint(p,'plate'));state.skyTargets=(state.index.janus_sky_targets||[]).map(p=>skyPoint(p,'janus'));}else{state.earthPoints=fallbackEarth.map(decorateEarth);const csv=tasks[4].status==='fulfilled'?tasks[4].value:fallbackPlateCSV;state.platePoints=parsePlateCSV(csv);state.skyTargets=[];}
    if(!state.earthPoints.length)state.earthPoints=fallbackEarth.map(decorateEarth);if(!state.platePoints.length)state.platePoints=parsePlateCSV(fallbackPlateCSV);state.spacePoints=[...state.platePoints,...state.skyTargets];const unmapped=state.index?.counts?.earth_unmapped_source_rows||0;$('earth-count').textContent=`${state.earthPoints.length} MAP · ${unmapped} OPEN`;buildEventList(state.earthPoints,unmapped);
    const cohort=state.space?.rows?.find(r=>r.row_kind==='ANALYSIS_COHORT'),dataset=state.space?.rows?.find(r=>r.row_kind==='PUBLIC_DATASET_RELEASE');$('m-plates').textContent=safe(dataset?.plate_count||cohort?.plates_in_window||state.platePoints.length);$('m-tiles').textContent=safe(dataset?.tile_rows);const cc=dataset?.catalogue_rows||cohort?.candidate_rows_on_in_window_plates;$('m-candidates').textContent=cc?Number(cc).toLocaleString():'—';$('m-days').textContent=safe(cohort?.unique_observation_dates);$('space-count').textContent=`${state.platePoints.length} PLATES · ${state.skyTargets.length} JANUS`;
    const gate=state.index?.bridge?.exact_day_join_gate;$('bridge-status').textContent=gate==='OPEN'?'L2 EXACT-DAY JOIN: OPEN':'L2 EXACT-DAY JOIN: BLOCKED';const ms=state.index?.bridge?.match_state_counts||{};$('bridge-reason').textContent=`Palomar 312-day denominator: ${state.index?.bridge?.palomar_day_rows_materialized||0}/312. Current Earth match states: ${Object.entries(ms).map(([k,v])=>`${k}=${v}`).join(' · ')||'pending index'}.`;initSpheres();setMode('earth');$('load-state').textContent=state.index?'LIVE AUTOSYNC':'RESILIENT FALLBACK';$('load-state').style.color=state.index?'#48e59a':'#ffbd5a';
  }

  load().catch(err=>{console.error(err);$('load-state').textContent='DEGRADED';$('load-state').style.color='#ff6b6b';state.earthPoints=fallbackEarth.map(decorateEarth);state.platePoints=parsePlateCSV(fallbackPlateCSV);state.spacePoints=state.platePoints;buildEventList(state.earthPoints,0);initSpheres();setMode('earth');});
})();
