// JANUS Cosmos Volume — deterministic model and source-bound scientific coordinates.
// A galaxy illustration is NOT a 3D HI4PI reconstruction. Astronomical event coordinates
// are displayed ONLY in a labeled feature-space; sky angular samples use a unit shell.
export const FERMI_SOURCE = "data/cosmos/COSMOS-ALI-V3-FERMI-PUBLISHED-HE-PHOTONS-2026-10-18.json";
export const FERMI_CATALOG = "https://arxiv.org/abs/1906.11403";
export const HI_SCHEMA = "janus.cosmos.hi.angular_samples.v1";
export const scientificLimits = Object.freeze({
  galaxy: "ILLUSTRATIVE_DENSITY_MODEL_NOT_OBSERVATION",
  fermi: "OBSERVATIONS_IN_FEATURE_SPACE_NOT_ASTROMETRIC_POSITIONS",
  hi: "OBSERVED_ANGULAR_SHELL_NO_DISTANCE_OR_TRUE_VOLUME",
  interpretation: "NO_FASTER_THAN_LIGHT_OR_RETROCAUSAL_INFERENCE"
});
export function rng(seed=73019) {
  let s=seed>>>0;
  return function () { s+=0x6D2B79F5; let t=s; t=Math.imul(t^(t>>>15),t|1);
    t^=t+Math.imul(t^(t>>>7),t|61); return ((t^(t>>>14))>>>0)/4294967296; };
}
function gaussian(random) {
  const u=Math.max(1e-10,random()),v=random();
  return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v);
}
function point(x,y,z,r,g,b,size,alpha=1,extra={}) {
  if(![x,y,z,r,g,b,size,alpha].every(Number.isFinite))throw Error("Non-finite scene value");
  return {x,y,z,r,g,b,size,alpha,...extra};
}
// Stable visualization-only galaxy particles; never present them as observational samples.
export function makeIllustration({gasCount=18000,starCount=1700,seed=20261114}={}) {
  if(!Number.isInteger(gasCount)||gasCount<0||gasCount>60000)throw Error("Invalid gasCount");
  if(!Number.isInteger(starCount)||starCount<0||starCount>10000)throw Error("Invalid starCount");
  const rand=rng(seed),particles=[];
  for(let i=0;i<gasCount;i++){
    const arm=Math.floor(rand()*4),rad=.35+4.45*Math.pow(rand(),.73);
    const phase=arm*Math.PI/2+rad*1.47+gaussian(rand)*(.16+.035*rad);
    const x=Math.cos(phase)*rad+gaussian(rand)*.085,z=Math.sin(phase)*rad+gaussian(rand)*.085;
    const y=gaussian(rand)*(.085+.047*rad);
    const core=Math.exp(-rad*.48),p=rand();
    const warm=Math.max(0,Math.min(1,.18+.40*Math.sin(rad*1.38+arm)+p*.22));
    particles.push(point(x,y,z,.19+.78*warm+.55*core,.33+.35*(1-warm)+.36*core,
      .68+.27*(1-warm),1.5+rand()*3.5,.10+rand()*.26,{kind:"model_gas"}));
  }
  for(let i=0;i<starCount;i++){
    const rad=Math.pow(rand(),.64)*5.1,theta=rand()*Math.PI*2;
    const y=gaussian(rand)*(.12+.05*rad);
    const hue=rand(),c=hue>.75?[1,.75,.48]:hue>.3?[.63,.82,1]:[.9,.95,1];
    particles.push(point(Math.cos(theta)*rad,y,Math.sin(theta)*rad,...c,
      1.1+rand()*2.6,.7+rand()*.3,{kind:"model_star"}));
  }
  return particles;
}
export function parseFermiCatalog(data) {
  if(!data || data.schema!=="janus.cosmos.published-fermi-photon-summary.v3"||
      !Array.isArray(data.bursts)||data.bursts.length===0)throw Error("Unexpected Fermi catalog schema");
  const records=[];
  for(const burst of data.bursts){
    if(typeof burst.grb!=="string" || !Array.isArray(burst.events))throw Error("Invalid burst");
    const z=Number(burst.redshift);
    if(!Number.isFinite(z)||z<0||z>30)throw Error("Invalid redshift");
    burst.events.forEach((event,index)=>{
      const energy=Number(event.energy_GeV),time=Number(event.time_s);
      if(!Number.isFinite(energy)||energy<=0||!Number.isFinite(time)||time<=0)
        throw Error("Invalid photon event");
      records.push({
        id:burst.grb+":"+(index+1),burst:burst.grb,class:burst.class,
        z,energy,time,observed:true,
        // Feature dimensions are NOT physical spatial coordinates:
        x:(z/4.5-.5)*8.2,
        y:((Math.log10(time)+.2)/4.8-.5)*6.0,
        zPosition:(energy/100-.5)*5.5
      });
    });
  }
  return records;
}
export function fermiPoints(events) {
  return events.map(e=>point(e.x,e.y,e.zPosition,.22+.74*(e.energy/100),
    .72-.30*(e.energy/100),1,.75+e.energy/22,.96,
    {kind:"fermi",event:e}));
}
// Explicitly angular, not physical 3D; NHI units and source record provenance mandatory.
// Samples JSON: {"schema":"janus.cosmos.hi.angular_samples.v1","source_url":"https://...",
// "nhi_unit":"cm-2","samples":[{"l_deg":0,"b_deg":0,"nhi_cm2":1e21},...]}.
export function parseHI(data) {
  if(!data||data.schema!==HI_SCHEMA||data.nhi_unit!=="cm-2"||
    typeof data.source_url!=="string"||!/^https:\/\//.test(data.source_url)||
    !Array.isArray(data.samples)||data.samples.length<1||data.samples.length>70000)
      throw Error("HI: source URL, schema, unit or sample count invalid");
  const rows=[];
  for(let i=0;i<data.samples.length;i++){
    const s=data.samples[i],l=Number(s.l_deg),b=Number(s.b_deg),nh=Number(s.nhi_cm2);
    if(![l,b,nh].every(Number.isFinite)||l<0||l>=360||b< -90||b>90||nh<=0)
      throw Error("HI invalid row "+i);
    const lon=l*Math.PI/180,lat=b*Math.PI/180,rad=4.25;
    const intensity=Math.max(0,Math.min(1,(Math.log10(nh)-19.2)/3.0));
    rows.push(point(rad*Math.cos(lat)*Math.cos(lon),rad*Math.sin(lat),
      rad*Math.cos(lat)*Math.sin(lon),.10+.88*intensity,.22+.65*intensity,
      1-.46*intensity,1.3+intensity*5.3,.20+.65*intensity,
      {kind:"hi_angular",l,b,nhi:nh,source:data.source_url}));
  }
  return rows;
}
export function bounds(records) {
  return {count:records.length,finite:records.every(p=>[p.x,p.y,p.z].every(Number.isFinite))};
}
