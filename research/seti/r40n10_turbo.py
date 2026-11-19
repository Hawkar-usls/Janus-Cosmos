"""Exact N9 scoring accelerator and signed-frequency repeat witness."""
import time
import numpy as np
from numba import njit
import r40n9_halo as ref9
@njit(cache=True,fastmath=False)
def legal_mask(paths,p):
    N=len(p);M=len(paths);mask=np.empty((N,M),np.bool_)
    for i in range(N):
        for k in range(M):
            num=0
            for j in range(16):
                if abs(paths[k,j]-p[i,j])>2:
                    num+=1
                    if num>=5:break
            mask[i,k]=num>=5
    return mask
@njit(cache=True,fastmath=False)
def bank_scores(U,paths):
    N=len(U);M=len(paths);out=np.empty((N,M),np.float64)
    for i in range(N):
        for k in range(M):
            if N==1: # NumPy's contiguous 16-element pairwise sum
                a0=U[i,0,paths[k,0]]+U[i,8,paths[k,8]]
                a1=U[i,1,paths[k,1]]+U[i,9,paths[k,9]]
                a2=U[i,2,paths[k,2]]+U[i,10,paths[k,10]]
                a3=U[i,3,paths[k,3]]+U[i,11,paths[k,11]]
                a4=U[i,4,paths[k,4]]+U[i,12,paths[k,12]]
                a5=U[i,5,paths[k,5]]+U[i,13,paths[k,13]]
                a6=U[i,6,paths[k,6]]+U[i,14,paths[k,14]]
                a7=U[i,7,paths[k,7]]+U[i,15,paths[k,15]]
                s=((a0+a1)+(a2+a3))+((a4+a5)+(a6+a7))
            else: # NumPy strided reduction in original advanced indexing layout
                s=0.
                for j in range(16):s+=U[i,j,paths[k,j]]
            out[i,k]=s/16.
    return out
ORACLE=ref9.context_scores
class Turbo:
    def __init__(self,paths,verify=6):
        self.paths=paths;self.offsets={o:np.ascontiguousarray(paths-o) for o in (0,16,32)}
        self.p=None;self.mask=None;self.verify=verify;self.calls=0;self.checked=0;self.fast_s=0.;self.reference_s=0.;self.paired_fast_s=0.;self.candidates=0
    def __call__(self,Z,p,offset,paths):
        assert paths is self.paths or np.array_equal(paths,self.paths)
        ts=time.perf_counter()
        pp=self.offsets[offset];ap=p+offset
        if self.p is None or not np.array_equal(ap,self.p):
            self.p=ap.copy();self.mask=legal_mask(self.paths,ap)
        U=ref9.grid(Z);sc=bank_scores(U,pp);mask=self.mask
        a=np.where(mask,sc,-np.inf);j=a.argmax(1);best=a[np.arange(len(Z)),j]
        assert np.isfinite(best).all() and np.allclose(ref9.contrast(Z,pp[j]).mean(1),best,rtol=1e-12,atol=1e-12)
        out=(sc,mask,best,self.paths[j]-ref9.HALO)
        dt=time.perf_counter()-ts;self.fast_s+=dt;self.calls+=1;self.candidates+=len(Z)*len(paths)
        if self.checked<self.verify:
            ts=time.perf_counter();expected=ORACLE(Z,p,offset,paths);self.reference_s+=time.perf_counter()-ts;self.paired_fast_s+=dt
            for a,b in zip(out,expected):assert np.array_equal(a,b),"R40N10 full score bank parity"
            self.checked+=1
        return out
    def report(self):
        return {"calls":self.calls,"oracle_calls":self.checked,"candidate_scores":self.candidates,"turbo_seconds":self.fast_s,"paired_turbo_seconds":self.paired_fast_s,"paired_reference_seconds":self.reference_s,"speedup":self.reference_s/self.paired_fast_s if self.paired_fast_s else None}
def selftest(norm,paths):
    rng=np.random.default_rng(1010)
    # Compile and warm before timing.
    bank_scores(np.zeros((2,16,176)),paths-32);legal_mask(paths,np.full((2,16),100))
    turbo=Turbo(paths,100)
    for N in (1,2,7,8):
        A=rng.normal(size=(N,16,208));B=rng.normal(size=(N,16,208));p=rng.integers(0,128,(N,16))
        for sh in ref9.SHIFTS:
            off=16+sh;pa=p+40-off
            for raw in (A,B):
                Z=norm(raw[:,:,off:off+176]);turbo(Z,pa,off,paths)
        # deliberately clipped values and exact ties
        Z=np.zeros((N,16,176));Z[:,:,70:72]=10.;turbo(Z,p+24,16,paths)
    return turbo.report()
def repeat_witness(rows,meta,pre,idx,p,ref,kind):
    # Distinct physical scans; same receiver-native frequency only, not source-frame extrapolation.
    N=len(idx);F=np.zeros((N,12));stress=np.zeros((N,12));valid=np.zeros(N,bool);stress_valid=np.zeros(N,bool)
    receipt=np.full((N,2),"",dtype="<U64")
    offsets=np.arange(-8,9);r=np.arange(16)[:,None]
    def features(A,B,path):
        scores=[];pos=[];best=[];sel=[];frac=[]
        for raw in (A,B):
            Z=ref.norm(raw[None,:,:]);P=path[:,None]+offsets[None,:]
            sides=Z[0,r[:,:,None],P[:,:,None]+np.array([-5,-4,-3,3,4,5])];V=np.clip(Z[0,r,P]-np.median(sides,2),-4,4)
            s=V.mean(0);j=int(s.argmax());scores.append(s);pos.append(j);best.append(float(s[j]));sel.append(float(s[8]));frac.append(float((V[:,8]>0).mean()))
        return np.asarray([sel[0],sel[1],best[0],best[1],frac[0],frac[1],sel[0]-sel[1],best[0]-best[1],abs(pos[0]-pos[1])/16.,min(best),float(np.minimum(scores[0],scores[1]).max()),float((scores[0]>0).mean()-(scores[1]>0).mean())])
    for j,i in enumerate(idx):
        slot,wi,f0,df=meta[i];oth=2-slot;ctx=[]
        for s in (oth,oth+1):
            a,X=pre[s][wi];vdf=a["df_hz"]
            if not np.isclose(df,vdf,rtol=1e-12,atol=1e-12):break
            off=(f0-a["channel0_mhz"])*1e6/vdf;k=int(round(off))
            if abs(off-k)>=.001 or k<40 or k+168>X.shape[1]:break
            ctx.append(X[:,k-40:k+168].astype(float).copy());receipt[j,len(ctx)-1]=a["slice_sha256"]
        if len(ctx)!=2:continue
        A,B=ctx;path=p[j]+40
        F[j]=features(A,B,path);valid[j]=True
        # Additional synthetic repeated-stationary stress only, no training: inject in other ON;
        # shared original challenge additionally injects into other OFF. Fixed native bin and original active rows.
        row=rows[i]
        if row[1]==1 and (row[7]==row[7][0]).all():
            truth=row[7]+40;active=row[8]
            def inject(X):
                xx=X[:,40:168];scale=1.4826*np.median(abs(xx-np.median(xx,axis=1)[:,None]),axis=1)
                assert np.all(scale>0)
                X=X.copy();X[np.arange(16)[active],truth[active]]+=row[4]*scale[active];return X
            stress[j]=features(inject(A),inject(B) if kind=="shared" else B,path);stress_valid[j]=True
    assert np.isfinite(F).all() and np.isfinite(stress).all()
    return {"repeat":F,"repeat_supported":valid,"repeat_stress":stress,"repeat_stress_supported":stress_valid,"repeat_source_sha":receipt}
