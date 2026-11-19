"""Native spectral context and ArtMoney-inspired physical candidate intersection."""
import numpy as np
HALO=40
SHIFTS=(-16,0,16)
DOMAIN=(37,171)
def contrast(Z,p):
    i=np.arange(len(Z))[:,None];j=np.arange(16)[None,:]
    assert p.min()>=5 and p.max()<Z.shape[2]-5
    side=Z[i[:,:,None],j[:,:,None],p[:,:,None]+np.array([-5,-4,-3,3,4,5])]
    return np.clip(Z[i,j,p]-np.median(side,axis=2),-4,4)
def grid(Z):
    i=np.arange(5,Z.shape[2]-5);out=np.zeros_like(Z)
    out[:,:,i]=np.clip(Z[:,:,i]-np.median(Z[:,:,i[:,None]+np.array([-5,-4,-3,3,4,5])],axis=3),-4,4)
    return out
def candidates(bank):
    paths=[]
    for delta in bank:
        lo=DOMAIN[0]-int(delta.min());hi=DOMAIN[1]-int(delta.max())
        if hi>lo:paths.append(np.arange(lo,hi)[:,None]+delta)
    return np.concatenate(paths,axis=0)
def context_scores(Z,p,offset,paths):
    U=grid(Z);N=len(Z);ppaths=paths-offset
    assert ppaths.min()>=5 and ppaths.max()<Z.shape[2]-5
    rows=np.arange(16);allsc=np.empty((N,len(paths)))
    for a in range(0,len(paths),1024):
        pp=ppaths[a:a+1024];allsc[:,a:a+len(pp)]=U[:,rows[None,:],pp].mean(2)
    legal=(abs(paths[None,:,:]-(p[:,None,:]+offset))>2).sum(2)>=5
    ss=np.where(legal,allsc,-np.inf)
    j=ss.argmax(1);best=ss[np.arange(N),j];bp=ppaths[j]
    assert np.isfinite(best).all() and np.allclose(contrast(Z,bp).mean(1),best,atol=1e-12,rtol=1e-12)
    return allsc,legal,best,paths[j]-HALO
def intersect(scores,legal,parents):
    keep=legal.copy();counts=[keep.sum(1)]
    for sc in scores:
        eligible=np.where(legal,sc,-np.inf)
        ranks=np.maximum(1,np.ceil(legal.sum(1)*.1).astype(int))
        cut=np.array([np.partition(eligible[i],-ranks[i])[-ranks[i]] for i in range(len(legal))])
        keep &= sc>=cut[:,None];counts.append(keep.sum(1))
    assert all(np.all(a>=b) for a,b in zip(counts,counts[1:]))
    worst=np.minimum.reduce([sc-parent[:,None] for sc,parent in zip(scores,parents)])
    margin=np.max(np.where(keep,worst,-np.inf),axis=1)
    margin=np.where(keep.any(1),margin,0.)
    return counts[-1]/np.maximum(1,counts[0]),margin,np.stack(counts,1)
def halo_features(A,B,p,norm,paths):
    assert A.shape==B.shape and A.shape[1:]==(16,208)
    N=len(A);ap=p+HALO;rr=np.arange(16)[None,:];ii=np.arange(N)[:,None]
    oldvalid=(p>=5)&(p<123)
    supports=np.stack([((p-sh>=5)&(p-sh<123)).sum(1) for sh in SHIFTS],1)
    ft=[];art=[];traces=[];winners=[];margins=[];contexts=[];counts=[]
    for raw in (A,B):
        z=norm(raw);fixed=contrast(z,ap);missing=~oldvalid
        restored=np.where(missing.any(1),(fixed*missing).sum(1)/np.maximum(1,missing.sum(1)),fixed.mean(1))
        vals=[];rv=[];pp=[];scorebank=[];mask=None
        for sh in SHIFTS:
            start=16+sh;stop=start+176;path=ap-start
            assert np.array_equal(contrast(z[:,:,start:stop],path),fixed),"Fixed-scale coordinate parity"
            local=norm(raw[:,:,start:stop]);sc=contrast(local,path).mean(1)
            scores,legal,rs,rp=context_scores(local,path,start,paths)
            if mask is None:mask=legal
            else:assert np.array_equal(mask,legal)
            vals.append(sc);rv.append(sc-rs);pp.append(rp);scorebank.append(scores)
        fraction,margin,ct=intersect(scorebank,mask,vals)
        vals=np.stack(vals,1);rv=np.stack(rv,1);pp=np.stack(pp,1)
        near=ap[:,:,None]+np.arange(-2,3)
        residual=np.abs(z[ii[:,:,None],rr[:,:,None],near].argmax(2)-2).mean(1)
        disagree=(np.mean(abs(pp[:,0]-pp[:,1]),1)+np.mean(abs(pp[:,0]-pp[:,2]),1)+np.mean(abs(pp[:,1]-pp[:,2]),1))/3
        ft.extend([fixed.mean(1),np.sort(fixed.reshape(N,4,4).mean(2),axis=1)[:,1],restored,residual,np.median(rv,1),np.ptp(vals,axis=1),np.ptp(rv,axis=1),disagree,missing.mean(1)])
        art.extend([fraction,margin]);counts.append(ct);traces.append(fixed);winners.append(pp);margins.append(rv);contexts.append(vals)
    F=np.column_stack(ft);G=np.column_stack(art)
    assert F.shape==(N,18) and G.shape==(N,4) and np.isfinite(F).all() and np.isfinite(G).all()
    return F,G,np.stack(traces,1),np.stack(winners,1),np.stack(margins,1),np.stack(contexts,1),supports,np.stack(counts,1)
def extract(rows,refs,meta,pre,idx,p,ref,sha,paths):
    N=len(idx);ok=np.zeros(N,bool);F=np.zeros((N,18));G=np.zeros((N,4));tr=np.zeros((N,2,16));rp=np.zeros((N,2,3,16),int)
    mg=np.zeros((N,2,3));cs=np.zeros((N,2,3));su=np.zeros((N,3),int);ct=np.zeros((N,2,4),int)
    hashes=np.full((N,2),"",dtype="<U64");blocks=[];refblocks=[];which=[]
    for j,i in enumerate(idx):
        slot,wi,f0,df=meta[i];a,X=pre[slot][wi];bb,R=pre[slot+1][wi]
        X,R,lo=ref.common(a,bb,X,R)
        assert np.isclose(df,a["df_hz"],rtol=1e-12,atol=1e-12)
        q=(f0-a["channel0_mhz"])*1e6/df-lo;k=int(round(q));assert abs(q-k)<.001 and 0<=k<=X.shape[1]-128
        rawx=X[:,k:k+128].astype(float);rawr=R[:,k:k+128].astype(float)
        xx=np.asarray(rows[i][0],float);rr=np.asarray(refs[i],float)
        inj=np.zeros((16,128),bool);inj[np.arange(16)[rows[i][8]],rows[i][7][rows[i][8]]]=bool(rows[i][1])
        assert np.array_equal(rawx[~inj],xx[~inj]) and np.array_equal(rawr[~inj],rr[~inj])
        if k<HALO or k+128+HALO>X.shape[1]:continue
        A=X[:,k-HALO:k+128+HALO].astype(float).copy();B=R[:,k-HALO:k+128+HALO].astype(float).copy()
        assert A.shape==B.shape==(16,208)
        hashes[j]=[sha(A.tobytes()),sha(B.tobytes())]
        A[:,HALO:HALO+128]=xx;B[:,HALO:HALO+128]=rr
        assert np.array_equal(A[:,HALO:HALO+128],xx) and np.array_equal(B[:,HALO:HALO+128],rr)
        assert np.array_equal(A[:,:HALO],X[:,k-HALO:k]) and np.array_equal(A[:,-HALO:],X[:,k+128:k+128+HALO])
        assert np.array_equal(B[:,:HALO],R[:,k-HALO:k]) and np.array_equal(B[:,-HALO:],R[:,k+128:k+128+HALO])
        ok[j]=True;which.append(j);blocks.append(A);refblocks.append(B)
    for t0 in range(0,len(which),8):
        Q=np.array(which[t0:t0+8]);vals=halo_features(np.stack(blocks[t0:t0+8]),np.stack(refblocks[t0:t0+8]),p[Q],ref.norm,paths)
        F[Q],G[Q],tr[Q],rp[Q],mg[Q],cs[Q],su[Q],ct[Q]=vals
    return {"halo":F,"artmoney":G,"halo_supported":ok,"halo_traces":tr,"halo_rival_paths":rp,"halo_margins":mg,"halo_context_scores":cs,"shift_core_support":su,"intersection_counts":ct,"halo_source_sha":hashes}
def selftest(norm,paths):
    rng=np.random.default_rng(909);A=rng.normal(size=(2,16,208));B=rng.normal(size=(2,16,208));p=np.tile(np.arange(112,128),(2,1));p[1]=np.arange(16)
    A[np.arange(2)[:,None],np.arange(16)[None,:],p+HALO]+=8
    F,G,tr,rp,mg,cs,su,ct=halo_features(A,B,p,norm,paths)
    assert np.all(tr[:,0]>1) and (su[:,1]<16).all()
    sc=[np.array([[0.,2.,3.]]),np.array([[0.,3.,2.]]),np.array([[0.,2.,3.]])]
    frac,margin,ct0=intersect(sc,np.ones((1,3),bool),[np.zeros(1)]*3)
    assert ct0.tolist()==[[3,1,0,0]] and frac[0]==0 and margin[0]==0
    # A measured zero remains eligible if it satisfies the rank criterion.
    _,_,zero=intersect([np.zeros((1,3))]*3,np.array([[True,False,True]]),[np.zeros(1)]*3)
    assert zero.tolist()==[[2,2,2,2]]
    return {"native_edge_support":"PASS","translation":"PASS","successive_intersection":"PASS","measured_zero_retained":"PASS","direct_rival_recalculation":"PASS"}
