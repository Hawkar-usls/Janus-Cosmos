/* JANUS COSMOS FERMI-V7: read-only, source-bound event inspector.
 * Scatter coordinates are observation FEATURES, never positions or propagation distances.
 * This UI does not perform any LIV inference or modify source evidence.
 */
const $=id=>document.getElementById(id);
const colors={"090510":"#6cdaef","090902B":"#efc07b","130427A":"#baa5ff"};
const state={dataset:null,events:[],filtered:[],selected:null,hitPoints:[],source:null};
const canvas=$("plot"),ctx=canvas.getContext("2d"),limits={minTime:0,maxTime:1,minEnergy:0,maxEnergy:2};
const fmt=(v,dec=3)=>Number(v).toLocaleString("ru-RU",{maximumFractionDigits:dec});
function string(v){return String(v==null?"—":v);}
function setText(id,value){$(id).textContent=string(value);}
function finite(v){return typeof v==="number"&&Number.isFinite(v);}
function validate(d){
  if(!d||d.schema!=="janus.cosmos.fermi.v7.event-evidence-browser.v1"||
     !Array.isArray(d.events)||!Array.isArray(d.sources)||d.sources.length!==3||
     d.events.length!==175||d.counts?.archive_query_roi_events!==6520||
     d.gates?.multi_burst_LIV_inference!=="BLOCKED"||
     d.gates?.independent_external_replication!=="NOT_DONE"||
     d.gates?.observational_LIV_p_value!==null||d.gates?.lattice_spacing_fit!==null)
    throw Error("Несовпадение схемы, числа событий или научных ограничений V7.");
  const keys=new Set();
  for(const e of d.events){
    const id=[e.grb,e.run_id,e.event_id].join(":");
    if(keys.has(id)||!finite(e.met_s)||!finite(e.after_trigger_s)||
      !finite(e.energy_MeV)||e.energy_MeV<1000||!finite(e.angular_sep_deg)||!colors[e.grb])
      throw Error("Повреждён идентификатор, энергия, время или числовое поле V7.");
    keys.add(id);
  }
  return d;
}
function axisTime(t){return Math.sign(t)*Math.log10(1+Math.abs(t));}
function energyGeV(e){return e.energy_MeV/1000;}
function axisEnergy(e){return Math.log10(energyGeV(e));}
function computeBounds(){
  const t=state.events.map(e=>axisTime(e.after_trigger_s)),y=state.events.map(axisEnergy);
  const a=Math.min(...t),b=Math.max(...t),c=Math.min(...y),d=Math.max(...y);
  limits.minTime=a-.08;limits.maxTime=b+.08;
  limits.minEnergy=Math.min(0,c-.12);limits.maxEnergy=Math.max(2.1,d+.12);
}
function chartBounds(){
  const r=canvas.getBoundingClientRect(),ratio=Math.min(window.devicePixelRatio||1,2);
  const width=Math.max(1,Math.round(r.width*ratio)),height=Math.max(1,Math.round(r.height*ratio));
  if(canvas.width!==width||canvas.height!==height){canvas.width=width;canvas.height=height;}
  const scale=ratio;
  return {w:r.width,h:r.height,scale,x0:65,x1:r.width-18,y0:25,y1:r.height-53};
}
function scaleX(value,b){return b.x0+(axisTime(value)-limits.minTime)/(limits.maxTime-limits.minTime)*(b.x1-b.x0);}
function scaleY(value,b){return b.y1-(Math.log10(value)-limits.minEnergy)/(limits.maxEnergy-limits.minEnergy)*(b.y1-b.y0);}
function background(b){
  ctx.fillStyle="#071528";ctx.fillRect(0,0,b.w,b.h);
  const grad=ctx.createRadialGradient(b.w*.48,b.h*.38,0,b.w*.48,b.h*.38,b.w*.68);
  grad.addColorStop(0,"rgba(46,99,149,.19)");grad.addColorStop(1,"rgba(5,11,28,0)");
  ctx.fillStyle=grad;ctx.fillRect(0,0,b.w,b.h);
  ctx.font="10px ui-monospace,Consolas,monospace";
  ctx.textBaseline="middle";
  for(const val of [1,2,3,5,10,20,50,100,200,300]){
    if(axisEnergy({energy_MeV:val*1000})<limits.minEnergy||
       Math.log10(val)>limits.maxEnergy)continue;
    const y=scaleY(val,b);
    ctx.strokeStyle="rgba(126,176,231,.12)";ctx.lineWidth=1;
    ctx.beginPath();ctx.moveTo(b.x0,y);ctx.lineTo(b.x1,y);ctx.stroke();
    ctx.fillStyle="#8aa7c8";ctx.textAlign="right";ctx.fillText(string(val),b.x0-9,y);
  }
  for(const sec of [-100,-10,-1,0,1,10,100,1000,10000]){
    const x=scaleX(sec,b);
    if(x<b.x0||x>b.x1)continue;
    ctx.strokeStyle="rgba(130,183,233,.12)";ctx.beginPath();
    ctx.moveTo(x,b.y0);ctx.lineTo(x,b.y1);ctx.stroke();
    ctx.fillStyle="#86a3c0";ctx.textAlign="center";ctx.fillText(string(sec),x,b.y1+22);
  }
  ctx.strokeStyle="rgba(147,195,240,.43)";ctx.beginPath();
  ctx.moveTo(b.x0,b.y0);ctx.lineTo(b.x0,b.y1);ctx.lineTo(b.x1,b.y1);ctx.stroke();
}
function renderPlot(){
  if(!ctx||!state.dataset)return;
  const b=chartBounds();ctx.setTransform(b.scale,0,0,b.scale,0,0);ctx.clearRect(0,0,b.w,b.h);
  background(b);state.hitPoints=[];
  for(const e of state.filtered){
    const x=scaleX(e.after_trigger_s,b),y=scaleY(energyGeV(e),b);
    const selected=state.selected&&key(e)===key(state.selected),radius=selected?6:3.8;
    if(x<b.x0||x>b.x1||y<b.y0||y>b.y1)continue;
    ctx.beginPath();ctx.arc(x,y,radius+2,0,2*Math.PI);
    ctx.fillStyle=selected?"rgba(255,255,255,.38)":"rgba(93,165,226,.07)";ctx.fill();
    ctx.beginPath();ctx.arc(x,y,radius,0,Math.PI*2);
    ctx.fillStyle=colors[e.grb];ctx.globalAlpha=selected?1:.8;ctx.fill();
    ctx.globalAlpha=1;
    if(selected){
      ctx.strokeStyle="#eaf9ff";ctx.lineWidth=1.3;ctx.stroke();
    }
    state.hitPoints.push({event:e,x,y});
  }
}
function key(e){return e.grb+":"+e.run_id+":"+e.event_id;}
function selectEvent(e,fromTable=false){
  state.selected=e;
  setText("selected-identity",key(e));
  $("event-detail").textContent=[
    "GRB: "+e.grb+" · RUN_ID "+e.run_id+" · EVENT_ID "+e.event_id,
    "LAT TIME (MET): "+fmt(e.met_s,7)+" s",
    "От оценочного GBM-триггера: "+fmt(e.after_trigger_s,6)+" s",
    "Реконструированная энергия: "+fmt(energyGeV(e),6)+" GeV",
    "Угловое расстояние от центра: "+fmt(e.angular_sep_deg,6)+"°",
    "THETA: "+fmt(e.theta_deg,6)+"° · ZENITH: "+fmt(e.zenith_deg,6)+"°",
    "EVENT_CLASS mask: "+e.event_class_mask+" · EVENT_TYPE mask: "+e.event_type_mask,
    "Источник: LAT Extended Pass8 · V4 diagnostic export, NOT event association",
    "Интерпретация: arrival event; source emission time UNKNOWN."
  ].join("\n");
  showSource(e.grb);
  renderPlot();
  if(!fromTable)renderTable();
}
function showSource(grb){
  const source=state.dataset?.sources.find(x=>x.grb===grb);
  if(!source)return;
  setText("source-title","GRB "+grb+" · "+source.class);
  setText("source-criteria",source.counts.diagnostic_selected+" отобранных из "+
    fmt(source.counts.events_in_query_roi,0)+" записей в области поиска. "+
    "GBM trigger MET = "+source.trigger_met_estimated+" (оценка).");
  setText("source-hash",source.event_file.sha256);
  const link=$("source-fits");
  // Only NASA FSSC URLs; do not allow source JSON to inject arbitrary navigations.
  if(/^https:\/\/fermi\.gsfc\.nasa\.gov\/FTP\/fermi\/data\/lat\/queries\/L[A-Z0-9]+_EV00\.fits$/.test(source.event_file.url)){
    link.href=source.event_file.url;link.textContent="Архивный EV00 FITS · может устареть";
  }else{
    link.href="https://fermi.gsfc.nasa.gov/ssc/data/access/";
    link.textContent="NASA FSSC · искать оригинал";
  }
}
function eventMatches(e){
  const grb=$("burst").value;
  if(grb!=="all"&&e.grb!==grb)return false;
  if(energyGeV(e)<Number($("min-energy").value))return false;
  const until=$("time-window").value;
  if(until!=="all"&&!(e.after_trigger_s>=0&&e.after_trigger_s<=Number(until)))return false;
  const q=$("event-search").value.trim().toLowerCase();
  return !q||[e.grb,e.run_id,e.event_id,key(e)].some(v=>string(v).toLowerCase().includes(q));
}
function filtered(){
  state.filtered=state.events.filter(eventMatches);
  setText("filtered-count",state.filtered.length.toLocaleString("ru-RU"));
  setText("plot-meta",state.filtered.length+" / 175");
  setText("ledger-count",state.filtered.length+" записей");
  setText("min-energy-value",$("min-energy").value+" ГэВ");
  if(state.selected&&!state.filtered.some(e=>key(e)===key(state.selected)))state.selected=null;
  if(!state.selected){setText("selected-identity","NOT SELECTED");
    setText("event-detail",state.filtered.length?"Нажми на точку или строку, чтобы открыть запись.":
      "Совпадений по выбранным фильтрам нет.");}
  renderPlot();renderTable();
  const s=state.filtered[0]||state.events[0];
  if(s)showSource(s.grb);
}
function renderTable(){
  const body=$("ledger");body.replaceChildren();
  if(!state.filtered.length){
    const tr=document.createElement("tr"),td=document.createElement("td");td.colSpan=9;
    td.textContent="Нет событий с такими фильтрами.";tr.append(td);body.append(tr);return;
  }
  for(const e of state.filtered){
    const tr=document.createElement("tr");
    tr.tabIndex=0;tr.setAttribute("aria-label","Открыть событие "+key(e));
    if(state.selected&&key(e)===key(state.selected))tr.classList.add("selected");
    const values=[e.grb,e.run_id,e.event_id,fmt(e.met_s,6),fmt(e.after_trigger_s,6),
      fmt(energyGeV(e),5),fmt(e.angular_sep_deg,5),fmt(e.theta_deg,4),fmt(e.zenith_deg,4)];
    for(const v of values){const td=document.createElement("td");td.textContent=string(v);tr.append(td);}
    tr.addEventListener("click",()=>selectEvent(e,true));
    tr.addEventListener("keydown",ev=>{if(ev.key==="Enter"||ev.key===" "){ev.preventDefault();selectEvent(e,true);}});
    body.append(tr);
  }
}
canvas.addEventListener("click",event=>{
  const r=canvas.getBoundingClientRect(),x=event.clientX-r.left,y=event.clientY-r.top;
  let best=null,dist=15*15;
  for(const p of state.hitPoints){
    const d=(p.x-x)**2+(p.y-y)**2;
    if(d<dist){dist=d;best=p.event;}
  }
  if(best)selectEvent(best);
});
canvas.addEventListener("keydown",e=>{
  if(e.key==="Enter"&&state.filtered.length){selectEvent(state.filtered[0]);e.preventDefault();}
});
window.addEventListener("resize",renderPlot);
for(const id of ["burst","min-energy","time-window","event-search"]){
  $(id).addEventListener(id==="event-search"||id==="min-energy"?"input":"change",filtered);
}
$("reset").addEventListener("click",()=>{
  $("burst").value="all";$("min-energy").value="1";$("time-window").value="all";$("event-search").value="";
  filtered();
});
window.__JANUS_FERMI_V7_DIAG=()=>({
  loaded:Boolean(state.dataset),sourceEvents:state.events.length,
  filteredEvents:state.filtered.length,burstCount:state.dataset?.counts?.grbs||0,
  rawStatus:state.dataset?.gates?.source_byte_recheck_this_build||"",
  scientificStatus:state.dataset?.gates?.multi_burst_LIV_inference||"",
  selectedEvent:state.selected?key(state.selected):null
});
async function start(){
  try{
    const response=await fetch("./data/v7-events.json",{credentials:"omit",cache:"no-store"});
    if(!response.ok)throw Error("HTTP "+response.status);
    const data=validate(await response.json());
    state.dataset=data;state.events=data.events;
    computeBounds();
    setText("metric-bursts",data.counts.grbs);
    setText("metric-roi",fmt(data.counts.archive_query_roi_events,0));
    setText("metric-selected",fmt(data.counts.diagnostic_selected,0));
    setText("gate-manifest","JSON PASS");
    setText("gate-events","JSON PASS");
    setText("gate-raw",data.gates.source_byte_recheck_this_build==="RAW_FITS_NOT_RECHECKED_THIS_BUILD"?
      "НЕ ПОВТОРЯЛАСЬ":"SHA + REPLAY PASS");
    setText("gate-note",data.gates.source_byte_recheck_this_build+
      " · сводка построена по исторически сохранённым событиям V4. "+
      "Второй программный контроль и тесты в CI не являются внешней лабораторной проверкой.");
    for(const s of data.sources){
      const option=document.createElement("option");option.value=s.grb;
      option.textContent="GRB "+s.grb+" · "+s.counts.diagnostic_selected+" кандидатов";
      $("burst").append(option);
    }
    filtered();
  }catch(error){
    const msg="FERMI-V7 БЛОКИРОВАН: "+error.message;
    setText("plot-meta","DATA BLOCKED");setText("event-detail",msg);
    setText("gate-note",msg);setText("ledger-count",msg);
    $("ledger").replaceChildren();
  }
}
start();
