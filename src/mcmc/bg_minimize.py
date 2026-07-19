#!/usr/bin/env python3
"""Direct, trustworthy background-level (DESI BAO + Pantheon+ SN) minimization for
LCDM and SEDE on the SAME mochi engine. 3 params: H0, omega_b, omega_cdm.
Reports chi2 and Delta chi2(SEDE-LCDM). Uses a real Nelder-Mead simplex + a
Powell refine, from multiple starts, so we trust the minima."""
import numpy as np
from scipy.optimize import minimize
from cobaya.model import get_model
from cobaya.yaml import yaml_load_file

def build(tag):
    info = yaml_load_file(f"mcmc/{tag}_bg.yaml")
    info["packages_path"] = "mcmc/packages"
    info.pop("sampler", None); info.pop("output", None)
    return get_model(info)

def chi2_fn(model):
    names = ["H0","omega_b","omega_cdm"]
    def f(x):
        p = dict(zip(names, x))
        if not (60<p["H0"]<80 and 0.019<p["omega_b"]<0.025 and 0.10<p["omega_cdm"]<0.14):
            return 1e6
        try:
            ll = model.loglikes(p, as_dict=True)[0]
            tot = sum(v for v in ll.values())
            chi2 = -2.0*tot
            chi2 += ((p["omega_b"]-0.02218)/0.00055)**2   # BBN prior (Schoneberg+ 2021), identical to both models
            return chi2
        except Exception:
            return 1e6
    return f

def best(model, label):
    f = chi2_fn(model)
    starts = [[68.0,0.02237,0.117],[70.0,0.0224,0.12],[66.5,0.0223,0.113],[69.0,0.0225,0.118]]
    results=[]
    for s in starts:
        r = minimize(f, s, method="Nelder-Mead",
                     options=dict(xatol=1e-4, fatol=1e-3, maxfev=1200))
        r2 = minimize(f, r.x, method="Powell", options=dict(xtol=1e-5, ftol=1e-5, maxfev=2000))
        results.append((r2.fun, r2.x))
    results.sort(key=lambda t: t[0])
    chi2, x = results[0]
    print(f"[{label}] chi2_min = {chi2:.3f}  at H0={x[0]:.3f} omega_b={x[1]:.5f} omega_cdm={x[2]:.4f}")
    print(f"[{label}] spread over starts: {[f'{r[0]:.3f}' for r in results]}")
    # per-likelihood breakdown at bestfit
    ll = model.loglikes(dict(zip(["H0","omega_b","omega_cdm"], x)), as_dict=True)[0]
    for k,v in ll.items(): print(f"    chi2__{k} = {-2*v:.3f}")
    print(f"    chi2__BBN(omega_b) = {((x[1]-0.02218)/0.00055)**2:.3f}")
    return chi2, x

print("="*66)
mL = build("lcdm"); cL,xL = best(mL, "LCDM")
print("-"*66)
mS = build("sede"); cS,xS = best(mS, "SEDE")
print("="*66)
print(f"Delta chi2 (SEDE - LCDM) [BAO+SN background] = {cS-cL:+.3f}")
