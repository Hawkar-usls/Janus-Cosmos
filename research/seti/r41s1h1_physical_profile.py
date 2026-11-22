"""Preregistered synthetic physics stress test. No archive labels or production decisions."""
import hashlib,json,pathlib,time
import numpy as np
from scipy.special import erf
from scipy.ndimage import gaussian_filter
from sklearn.metrics import roc_auc_score
P=pathlib.Path("research/seti/r41s1h1_physical_protocol.json")
S=json.loads(P.read_text())
O=pathlib.Path("runtime/r41s1h1");O.mkdir(parents=True,exist_ok=True)
def save(k,v):(O/k).write_text(json.dumps(v,indent=2,allow_nan=False)+"\n")
def sha(x):return hashlib.sha256(x).hexdigest()
def integrated(x,mu,sigma,df):
    z=np.sqrt(2)*sigma
    return .5*(erf((x+df/2-mu)/z)-erf((x-df/2-mu)/z))
def templates(df,dt):
    t=(np.arange(16)-7.5)*dt;x=np.arange(-24,25)*df
    bank=[];params=[]
    for width in (.25,.6,1.2,2.5,5.,10.):
        for off in (-6,-3,0,3,6):
            for drift in (-.10,0,.10):
                for acc in (-.001,0,.001):
                    mid=off*df+drift*t+.5*acc*t*t
                    mu=mid[:,None,None]+(drift+acc*t)[:,None,None]*np.linspace(-dt/2,dt/2,5)[None,:,None]
                    v=integrated(x[None,None,:],mu,width*df,df).mean(1)
                    v/=max(v.max(),1e-12)
                    bank.append(v);params.append((width,off,drift,acc))
    return np.asarray(bank),np.asarray(params)
def normalize(y):
    side=np.r_[0:7,42:49];bg=np.median(y[:,side],axis=1,keepdims=True)
    v=y-bg;mad=1.4826*np.median(abs(v[:,side]-np.median(v[:,side])))
    assert np.isfinite(mad) and mad>1e-10
    return np.clip(v/mad,-12,12)
def role_score(y,bank,params):
    y=normalize(y);ret=[]
    for fold in (0,1):
        tr=(np.arange(16)//4)%2==fold;te=~tr
        a=bank[:,tr,:];b=y[tr];norm=(a*a).sum((1,2))
        amp=np.maximum(0,np.einsum("btf,tf->b",a,b)/np.maximum(norm,1e-12))
        gain=(2*amp*np.einsum("btf,tf->b",a,b)-amp**2*norm)/b.size
        # predeclared families; choose only on training blocks.
        masks=[(params[:,0]<=1.2)&(abs(params[:,1])<=3),params[:,0]>=5,abs(params[:,1])>=6,(params[:,0]==.6)&(params[:,1]==0)&(params[:,2]==0)&(params[:,3]==0)]
        vals=[];widths=[]
        for m in masks:
            j=np.flatnonzero(m)[np.argmax(gain[m])]
            pred=amp[j]*bank[j,te];w=y[te]
            # held-out reduction in Huber loss; no test-block fitting.
            def huber(v):
                v=abs(v);return np.where(v<=2,.5*v*v,2*v-2)
            vals.append(float(np.mean(huber(w)-huber(w-pred))))
            widths.append(float(params[j,0]))
        ret.append((vals,widths))
    vals=np.asarray([r[0] for r in ret]);n,b,r,h=vals.mean(0)
    return {"narrow":float(n),"broad":float(b),"rival":float(r),"fixed":float(h),"width_bins":float(np.mean([z[1][0] for z in ret])),"narrow_conservative":float(vals[:,0].min())}
def score(cube,bank,params,df):
    on,off,repeat,roff=[role_score(y,bank,params) for y in cube]
    # Conservative paired evidence; scores are diagnostics, not a veto.
    def excess(q):return q["narrow"]-max(0,q["broad"],q["rival"])
    s=min(excess(on),excess(repeat))-max(0,off["narrow"],roff["narrow"])
    fixed=min(on["fixed"],repeat["fixed"])-max(0,off["fixed"],roff["fixed"])
    return {"physical":float(s),"fixed":float(fixed),"width_hz_on":on["width_bins"]*df,"width_hz_repeat":repeat["width_bins"]*df,"on_off":on["narrow"]-off["narrow"],"repeat_on_off":repeat["narrow"]-roff["narrow"]}
def generate(rng,kind,df,dt,held):
    t=(np.arange(16)-7.5)*dt;x=np.arange(-24,25)*df
    # Correlated receiver noise and bandpass slopes; held set includes unseen heavier correlation.
    z=rng.normal(size=(4,16,49));corr=gaussian_filter(rng.normal(size=z.shape),sigma=(0,.7 if held else .35,1.4 if held else .7))
    z=z+.9*corr/max(corr.std(),1e-9)+rng.normal(0,.1,(4,16,1))*x[None,None,:]/(24*df)
    amp=rng.uniform(2,7);wid=rng.uniform(.35,1.1)*df
    drift=rng.uniform(-.1,.1);acc=rng.uniform(-.001,.001);offset=rng.uniform(-1,1)*df
    roles=(0,2);intermittent=False
    if kind=="noise":return z
    if kind=="target_weak":amp=rng.uniform(.5,1.5)
    if kind=="target_intermittent":intermittent=True
    if kind=="shared_narrow":roles=(0,1,2,3)
    if kind=="broad_rfi":wid=rng.uniform(5,10)*df;roles=(0,2) if rng.random()<.5 else (0,1,2,3)
    if kind=="stationary_off":roles=(0,1,2,3);drift=acc=0;offset=0
    if kind=="shifted_rival":offset=rng.choice([-1,1])*rng.uniform(6,10)*df
    for k in roles:
        mu=offset+drift*t+.5*acc*t*t+rng.normal(0,.12*df)
        if held and rng.random()<.4:
            v=1/(1+((x[None,:]-mu[:,None])/wid)**2)
        else:
            v=integrated(x[None,:],mu[:,None],wid,df)
        v/=max(v.max(),1e-10)
        aa=amp*rng.uniform(.8,1.2,16)
        if intermittent:aa[rng.random(16)<.35]=0
        z[k]+=aa[:,None]*v
    return z
def evaluate(rows):
    out={}
    for seed in sorted({r["seed"] for r in rows}):
        rr=[r for r in rows if r["seed"]==seed]
        y=[r["positive"] for r in rr]
        out[str(seed)]={k:float(roc_auc_score(y,[r[k] for r in rr])) for k in ("physical","fixed")}
        out[str(seed)]["by_class"]={c:{k:float(np.median([r[k] for r in rr if r["class"]==c])) for k in ("physical","fixed","width_hz_on","width_hz_repeat")} for c in S["classes"]}
    return out
def main():
    t0=time.perf_counter();rows=[]
    contexts=[(1.,10.),(2.794,18.25),(5.,8.)];banks={v:templates(*v) for v in contexts}
    for split,key in (("development","development_seeds"),("heldout","heldout_seeds")):
        for seed in S["split"][key]:
            rng=np.random.default_rng(seed)
            for kind in S["classes"]:
                for i in range(S["split"]["examples_per_class_per_seed"]):
                    df,dt=contexts[int(rng.integers(len(contexts)))];bank,params=banks[(df,dt)]
                    cube=generate(rng,kind,df,dt,split=="heldout");q=score(cube,bank,params,df)
                    rows.append({"seed":seed,"split":split,"class":kind,"index":i,"positive":kind.startswith("target_"),"channel_hz":df,"integration_s":dt,**q})
            print("R41S1H1_SEED",seed,flush=True)
    assert len(rows)==1920
    held=[r for r in rows if r["split"]=="heldout"]
    y=[r["positive"] for r in held]
    auc={k:float(roc_auc_score(y,[r[k] for r in held])) for k in ("physical","fixed")}
    d={"status":"SYNTHETIC_DIAGNOSTIC_ONLY__HOLD_NATIVE_VALIDATION","protocol_sha256":sha(P.read_bytes()),"code_sha256":sha(pathlib.Path(__file__).read_bytes()),"heldout_examples":len(held),"development_examples":len(rows)-len(held),"heldout_auc":auc,"per_seed":evaluate(held),"development":evaluate([r for r in rows if r["split"]=="development"]),"decision_changes":0,"historical_tuning":False,"production_adoption":False,"external_validation":False,"elapsed_seconds":time.perf_counter()-t0}
    save("rows.json",rows);save("result.json",d);print("R41S1H1_SUMMARY",json.dumps({k:v for k,v in d.items() if k not in ("per_seed","development")}),flush=True)
if __name__=="__main__":main()
