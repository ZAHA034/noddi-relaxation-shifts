"""Verbatim record of the shell code that checked b2000 0 values, σ(DN) and S2 (post hoc, descriptive); this file itself was not run."""
import json, math, sys, numpy as np
sys.path.insert(0,".")
import b2k_summary as S
O=S.OUT
H=json.load(open(O/"headers.json"))
rows=[json.loads(l) for l in open(O/"read_sessions.jsonl")]
B=S.load_bank()
out={}
for u in S.VENDOR:
    rs=[r for r in rows if r["unit"]==u]
    cal=S.calib(B,u); Lp,g=S.g_curve(cal,"없음")
    sig=[];fe=[];sw=[];s2=[];exp0=[];N2=[]
    for r in rs:
        L,_=S.leff(r["R"],Lp,g); s=math.sqrt(r["mean_M2"]/(2*L)); sig.append(s)
        fe.append(r["diag"]["frac_excluded"]); N2.append(r["n_b2k"])
        S2=r["sens"].get("S2",{}); sw.append(S2.get("S_W")); 
        if S2.get("sigma_snr"): s2.append(s/S2["sigma_snr"])
        p0=1-math.exp(-0.25/(2*s*s)); exp0.append(1-(1-p0)**r["n_b2k"])   # Rayleigh (L 1) · probability that integer rounding gives a 0 in at least one of the 50 b2000 volumes
    hd=[H["sessions"][r["session"]] for r in rs]
    q=lambda a:(round(float(np.median(a)),3),round(float(np.min(a)),3),round(float(np.max(a)),3))
    out[u]=dict(n=len(rs),dtype=sorted({h["datatype"] for h in hd}),scl=sorted({str((h["scl_slope"],h["scl_inter"])) for h in hd}),
                sigma_DN=q(sig),zero_frac_obs=q(fe),zero_frac_rounding_expected=q(exp0),b0_WM_DN=q([x for x in sw if x]),S2_ratio=q(s2) if s2 else None)
    print(u, json.dumps(out[u],ensure_ascii=False))
json.dump(out,open("<SCRATCH>/zero_check.json","w"),ensure_ascii=False,indent=1)
