#!/usr/bin/env python3
"""Full-likelihood minimization: DESI DR2 BAO + Pantheon+ SN + Planck 2018
(lowl.TT + highl TTTEEE-lite) for LCDM and SEDE on the SAME mochi engine, with the
z_gr_smg=5 fix (CMB now sane across the whole box). Direct scipy multi-start (cobaya's
own Nelder-Mead stalls at the ref simplex). Minimizes {H0, omega_b, omega_cdm, logA, n_s};
tau fixed at 0.0544 (prior-pinned) and A_planck at 1.0025, identically for both models.
Reports chi2 and Delta chi2(SEDE-LCDM), with per-likelihood breakdown."""
import numpy as np, sys, time
from scipy.optimize import minimize
from cobaya.model import get_model
from cobaya.yaml import yaml_load_file

TAU=0.0544; APL=1.0025
NAMES=["H0","omega_b","omega_cdm","logA","n_s"]
BOUNDS=[(60,80),(0.019,0.025),(0.10,0.14),(2.9,3.2),(0.94,0.99)]

def build(tag):
    info=yaml_load_file(f"mcmc/{tag}.yaml"); info["packages_path"]="mcmc/packages"
    info.pop("sampler",None); info.pop("output",None)
    return get_model(info)

def chi2_fn(model):
    def f(x):
        for v,(lo,hi) in zip(x,BOUNDS):
            if not (lo<=v<=hi): return 1e7
        p=dict(zip(NAMES,x)); p["tau_reio"]=TAU; p["A_planck"]=APL
        try:
            ll=model.loglikes(p,as_dict=True)[0]
            c2=sum(-2*v for v in ll.values())
            if not np.isfinite(c2) or c2>1e6: return 1e7   # broken-region wall
            return c2
        except Exception:
            return 1e7
    return f

def best(model,label,starts):
    f=chi2_fn(model); res=[]
    for i,s in enumerate(starts):
        t0=time.time()
        r=minimize(f,s,method="Nelder-Mead",options=dict(xatol=2e-4,fatol=2e-3,maxfev=600))
        res.append((r.fun,r.x)); print(f"  [{label}] start {i}: chi2={r.fun:.3f} ({time.time()-t0:.0f}s) x={np.round(r.x,4)}",flush=True)
    res.sort(key=lambda t:t[0]); c2,x=res[0]
    p=dict(zip(NAMES,x)); p["tau_reio"]=TAU; p["A_planck"]=APL
    ll=model.loglikes(p,as_dict=True)[0]
    print(f"[{label}] BEST chi2={c2:.3f} at {dict(zip(NAMES,np.round(x,4)))}",flush=True)
    for k,v in ll.items(): print(f"    chi2__{k} = {-2*v:.3f}",flush=True)
    return c2,x,{k:-2*v for k,v in ll.items()}

# LCDM: 4 params (no z_gr); reuse same NAMES (LCDM ignores smg). start near known bestfit.
lcdm_starts=[[68.0,0.02237,0.1180,3.045,0.9649],[67.0,0.0225,0.1200,3.05,0.968]]
sede_starts=[[68.0,0.02237,0.1170,3.045,0.9649],[69.0,0.0224,0.1210,3.05,0.968]]

print("="*70,flush=True)
mL=build("lcdm"); cL,xL,dL=best(mL,"LCDM",lcdm_starts)
print("-"*70,flush=True)
mS=build("sede"); cS,xS,dS=best(mS,"SEDE",sede_starts)
print("="*70,flush=True)
print(f"Delta chi2 (SEDE - LCDM) [BAO+SN+Planck] = {cS-cL:+.3f}",flush=True)
print(f"  per-probe deltas: "+", ".join(f"{k.split('.')[-1]}={dS.get(k,0)-dL.get(k,0):+.2f}" for k in dL),flush=True)
