import {makeIllustration,parseFermiCatalog,fermiPoints,parseHI,scientificLimits} from "./model.mjs";
const $=id=>document.getElementById(id),canvas=$("cosmos");
let mode="galaxy",points=makeIllustration(),events=[],hiPoints=[],yaw=.4,pitch=.55,dist=12.2,slice=1,spinning=true,drag=null,last=0;
let gl,program,buffer; const MAX_DPR=1.5;
function status(v,bad=false){$("status").textContent=v;$("status").dataset.error=String(bad);}
const VERT=[
"#version 300 es","precision highp float;",
"layout(location=0) in vec3 aPosition;","layout(location=1) in vec3 aColor;",
"layout(location=2) in float aSize;","layout(location=3) in float aAlpha;",
"uniform float uYaw,uPitch,uDist,uAspect,uScale,uRatio,uSlice;",
"out vec3 vColor;out float vAlpha;",
"void main(){",
"float cy=cos(uYaw),sy=sin(uYaw),cp=cos(uPitch),sp=sin(uPitch);",
"vec3 q=vec3(cy*aPosition.x+sy*aPosition.z,aPosition.y,-sy*aPosition.x+cy*aPosition.z);",
"q=vec3(q.x,cp*q.y-sp*q.z,sp*q.y+cp*q.z);",
"float depth=max(.1,uDist-q.z);",
"bool cut=uSlice<.999 && aPosition.y>uSlice*6.-3.;",
"gl_Position=vec4(q.x*uScale/uAspect,q.y*uScale,depth*.48,depth);",
"gl_PointSize=cut?0.:clamp(aSize*uRatio*31./depth,1.,72.);",
"vColor=aColor;vAlpha=cut?0.:aAlpha;}"
].join("\n");
const FRAG=[
"#version 300 es","precision highp float;","in vec3 vColor;in float vAlpha;",
"uniform float uOpacity;out vec4 color;","void main(){",
"vec2 p=gl_PointCoord*2.-1.;float r=dot(p,p);if(r>1.)discard;",
"float glow=exp(-r*3.6);color=vec4(vColor*(.65+.5*glow),uOpacity*vAlpha*glow);}"
].join("\n");
function initGL(){
  try{
    gl=canvas.getContext("webgl2",{alpha:false,antialias:false});
    if(!gl)throw Error("WebGL2 unavailable");
    function compile(type,src){
      const shader=gl.createShader(type);gl.shaderSource(shader,src);gl.compileShader(shader);
      if(!gl.getShaderParameter(shader,gl.COMPILE_STATUS))throw Error(gl.getShaderInfoLog(shader));
      return shader;
    }
    program=gl.createProgram();const v=compile(gl.VERTEX_SHADER,VERT),f=compile(gl.FRAGMENT_SHADER,FRAG);
    gl.attachShader(program,v);gl.attachShader(program,f);gl.linkProgram(program);
    if(!gl.getProgramParameter(program,gl.LINK_STATUS))throw Error(gl.getProgramInfoLog(program));
    gl.deleteShader(v);gl.deleteShader(f);gl.useProgram(program);
    buffer=gl.createBuffer();gl.bindBuffer(gl.ARRAY_BUFFER,buffer);
    for(const [idx,n,offset] of [[0,3,0],[1,3,12],[2,1,24],[3,1,28]]){
      gl.enableVertexAttribArray(idx);gl.vertexAttribPointer(idx,n,gl.FLOAT,false,32,offset);
    }
    gl.disable(gl.DEPTH_TEST);gl.enable(gl.BLEND);gl.blendFunc(gl.SRC_ALPHA,gl.ONE);
    $("engine").textContent="WebGL2 · particle volume";
  }catch(e){
    gl=null;$("engine").textContent="Canvas 2D · fallback";
    status("GPU недоступен, включён 2D режим: "+e.message,true);
  }
}
function axes(pts){
  for(let i=-40;i<=40;i++)for(let axis=0;axis<3;axis++){
    const v=[0,0,0];v[axis]=i/10;
    pts.push({x:v[0],y:v[1],z:v[2],r:axis===0?.40:.15,g:axis===1?.70:.45,
      b:axis===2?.97:.70,size:1.1,alpha:.22,kind:"axis"});
  }
  return pts;
}
function upload(){
  $("points").textContent=points.length.toLocaleString("ru-RU");
  if(!gl)return;
  const raw=new Float32Array(points.length*8);
  points.forEach((p,i)=>{
    let j=i*8;raw[j++]=p.x;raw[j++]=p.y;raw[j++]=p.z;
    raw[j++]=p.r;raw[j++]=p.g;raw[j++]=p.b;raw[j++]=p.size;raw[j]=p.alpha;
  });
  gl.bindBuffer(gl.ARRAY_BUFFER,buffer);gl.bufferData(gl.ARRAY_BUFFER,raw,gl.STATIC_DRAW);
}
function choose(next){
  if(next==="hi"&&!hiPoints.length){$("hi-upload").click();return;}
  mode=next;
  document.querySelectorAll("[data-mode]").forEach(el=>{
    const yes=el.dataset.mode===mode;el.classList.toggle("active",yes);el.setAttribute("aria-pressed",String(yes));
  });
  if(mode==="galaxy"){
    points=makeIllustration();dist=12.2;pitch=.5;
    $("heading").textContent="Галактика · объемная модель";
    $("subheading").textContent="Водородные облака и спиральные рукава: визуальная схема, не измеренное положение газа.";
    $("axes").textContent="Модельные координаты · без физических единиц";
    $("scientific-note").textContent=scientificLimits.galaxy;
    $("source").href="https://www.astro.uni-bonn.de/hisurvey/hi4pi/";
    $("source").textContent="Научный контекст HI4PI (данные не встроены)";
    status("Галактический газ построен генератором, НЕ по наблюдениям.");
  }else if(mode==="fermi"){
    points=axes(fermiPoints(events));dist=11.8;pitch=.3;
    $("heading").textContent="Fermi LAT · высокоэнергетические фотоны";
    $("subheading").textContent="Оригинальные записи Fermi: нажмите на точку для энергии, времени и GRB.";
    $("axes").textContent="X: красное смещение z · Y: log10(сек после GBM) · Z: энергия GeV";
    $("scientific-note").textContent=scientificLimits.fermi+
      " · время триггера не равно моменту излучения источника.";
    $("source").href="https://github.com/Hawkar-usls/Janus-Cosmos/blob/6b902d4cd71971aacb4b504d680692fa798600a1/data/cosmos/COSMOS-ALI-V3-FERMI-PUBLISHED-HE-PHOTONS-2026-10-18.json";
    $("source").textContent="Исходный JSON Fermi · закреплённая версия";
    status(events.length?"Fermi: "+events.length+" событий из опубликованного каталога.":"Файл Fermi пока не загружен.",!events.length);
  }else{
    points=hiPoints;dist=12.4;pitch=.3;
    $("heading").textContent="HI · наблюдённый угловой слой";
    $("subheading").textContent="Направления на небе, а не реконструированные расстояния до облаков.";
    $("axes").textContent="l, b в галактических градусах · NHI в cm⁻²";
    $("scientific-note").textContent=scientificLimits.hi;
    status("Импортировано "+hiPoints.length+" направлений HI.");
  }
  $("detail").textContent=mode==="fermi"?"Выберите фотон для чтения наблюдения.":
    "Вращайте поле мышью или пальцем, увеличивайте колесом.";
  $("view").textContent=mode==="galaxy"?"ILLUSTRATIVE / NOT MEASURED":mode==="fermi"?
    "REAL DATA / FEATURE SPACE":"REAL DATA / ANGULAR SKY";
  upload();
}
function project(p,w,h){
  const cy=Math.cos(yaw),sy=Math.sin(yaw),cp=Math.cos(pitch),sp=Math.sin(pitch);
  const x=cy*p.x+sy*p.z,z0=-sy*p.x+cy*p.z;
  const y=cp*p.y-sp*z0,z=sp*p.y+cp*z0,depth=Math.max(.1,dist-z);
  return {x:w/2+x*2.35/depth*h/2,y:h/2-y*2.35/depth*h/2};
}
function resize(){
  const r=canvas.getBoundingClientRect(),d=Math.min(MAX_DPR,window.devicePixelRatio||1);
  const w=Math.max(1,Math.round(r.width*d)),h=Math.max(1,Math.round(r.height*d));
  if(canvas.width!==w||canvas.height!==h){canvas.width=w;canvas.height=h;}
  if(gl)gl.viewport(0,0,w,h);
}
function frame(t){
  if(spinning&&!window.matchMedia("(prefers-reduced-motion: reduce)").matches)
    yaw+=Math.min(50,t-(last||t))*.00011;
  last=t;resize();
  if(gl){
    gl.clearColor(.002,.005,.018,1);gl.clear(gl.COLOR_BUFFER_BIT);gl.useProgram(program);
    function u(k,v){gl.uniform1f(gl.getUniformLocation(program,k),v);}
    u("uYaw",yaw);u("uPitch",pitch);u("uDist",dist);u("uAspect",canvas.width/canvas.height);
    u("uScale",2.35);u("uRatio",Math.min(MAX_DPR,window.devicePixelRatio||1));
    u("uSlice",slice);u("uOpacity",Number($("brightness").value)/100);
    gl.drawArrays(gl.POINTS,0,points.length);
  }else{
    const ctx=canvas.getContext("2d");if(ctx){
      ctx.fillStyle="#020613";ctx.fillRect(0,0,canvas.width,canvas.height);
      const step=Math.max(1,Math.ceil(points.length/3500));
      for(let i=0;i<points.length;i+=step){
        const p=points[i];if(slice<.999&&p.y>slice*6-3)continue;
        const pos=project(p,canvas.width,canvas.height);
        ctx.fillStyle="rgba("+Math.round(p.r*255)+","+Math.round(p.g*255)+","+Math.round(p.b*255)+","+Math.min(.9,p.alpha)+")";
        ctx.beginPath();ctx.arc(pos.x,pos.y,Math.max(.6,p.size*.7),0,Math.PI*2);ctx.fill();
      }
    }
  }
  requestAnimationFrame(frame);
}
function pick(x,y){
  if(mode==="galaxy"){status("Эти частицы — графическая модель, а не астрономические измерения.");return;}
  const r=canvas.getBoundingClientRect(),scale=canvas.width/r.width;
  let best=null,distance=23*23*scale*scale;
  for(const p of points){
    if(mode==="fermi"?p.kind!=="fermi":p.kind!=="hi_angular")continue;
    if(slice<.999&&p.y>slice*6-3)continue;
    const q=project(p,canvas.width,canvas.height),v=(q.x-x*scale)**2+(q.y-y*scale)**2;
    if(v<distance){distance=v;best=p;}
  }
  if(!best){$("detail").textContent="Нажмите ближе к светящейся точке.";return;}
  if(best.kind==="fermi"){
    const e=best.event;
    $("detail").textContent="GRB "+e.burst+" · "+e.class+"\nЭнергия: "+e.energy+
      " GeV\nОт триггера GBM: "+e.time+" s\nКрасное смещение: z="+e.z+
      "\nПроисхождение: опубликованная таблица Fermi-LAT (Ajello et al.)"+
      "\nЗапрещённая интерпретация: время после триггера ≠ измерение скорости света.";
  }else{
    $("detail").textContent="HI · направление l="+best.l+"°, b="+best.b+
      "°\nNHI="+best.nhi.toExponential(4)+" cm⁻²\nИсточник: "+best.source+
      "\nРасстояние до газа неизвестно.";
  }
}
canvas.addEventListener("pointerdown",e=>{
  drag={x:e.clientX,y:e.clientY,yaw,pitch,moved:false};canvas.setPointerCapture(e.pointerId);
});
canvas.addEventListener("pointermove",e=>{
  if(!drag)return;
  const dx=e.clientX-drag.x,dy=e.clientY-drag.y;
  if(Math.abs(dx)+Math.abs(dy)>5)drag.moved=true;
  yaw=drag.yaw+dx*.006;pitch=Math.min(1.4,Math.max(-1.4,drag.pitch+dy*.006));
});
canvas.addEventListener("pointerup",e=>{
  if(!drag)return;
  if(!drag.moved){const r=canvas.getBoundingClientRect();pick(e.clientX-r.left,e.clientY-r.top);}
  drag=null;
});
canvas.addEventListener("wheel",e=>{
  e.preventDefault();dist=Math.max(7,Math.min(30,dist*Math.exp(e.deltaY*.001)));
},{passive:false});
document.querySelectorAll("[data-mode]").forEach(b=>b.addEventListener("click",()=>choose(b.dataset.mode)));
$("brightness").addEventListener("input",e=>$("brightness-value").textContent=e.target.value+"%");
$("slice").addEventListener("input",e=>{
  slice=Number(e.target.value)/100;$("slice-value").textContent=e.target.value+"%";
});
$("spin").addEventListener("click",()=>{
  spinning=!spinning;$("spin").textContent=spinning?"Остановить вращение":"Включить вращение";
});
$("reset").addEventListener("click",()=>{
  yaw=.4;pitch=.5;dist=12.2;slice=1;$("slice").value=100;$("slice-value").textContent="100%";
});
$("hi-upload").addEventListener("change",async e=>{
  const file=e.target.files?.[0];if(!file)return;
  try{
    if(file.size>14e6)throw Error("Слишком большой файл (максимум 14 MB)");
    const payload=JSON.parse(await file.text());
    hiPoints=parseHI(payload);choose("hi");
    $("source").href=payload.source_url;$("source").textContent="Оригинал импортированного слоя";
  }catch(error){status("HI JSON не принят: "+error.message,true);}
  e.target.value="";
});
async function loadFermi(){
  let error="";
  for(const path of ["../../data/cosmos/COSMOS-ALI-V3-FERMI-PUBLISHED-HE-PHOTONS-2026-10-18.json","../data/fermi.json"]){
    try{
      const r=await fetch(path,{cache:"no-store",credentials:"omit"});
      if(!r.ok)throw Error("HTTP "+r.status);
      const data=await r.json();
      events=parseFermiCatalog(data);
      $("fermi-count").textContent=events.length+" фотонов · "+data.bursts.length+" GRB";
      return;
    }catch(e){error=e.message;}
  }
  $("fermi-count").textContent="НЕТ ДАННЫХ";
  status("Оригинальный Fermi JSON не загружен: "+error,true);
}
canvas.addEventListener("webglcontextlost",e=>{
  e.preventDefault();status("Потерян графический контекст — обновите страницу.",true);
});
initGL();upload();requestAnimationFrame(frame);await loadFermi();
