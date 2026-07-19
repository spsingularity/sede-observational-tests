#!/usr/bin/env python3
"""
TIER-2 Stage-2 (finite-k): drive the SEDE FPAB fixed point through the mochi-class
modified-gravity Boltzmann code and extract the observables that DISTINGUISH SEDE
from LCDM -- the ones the compressed Stage-1.5 fit could not capture.

Uses the corpus's own SEDE->mochi stable-basis export (fpab_sede_b020_{stable,rho}.dat,
alpha_B0 = b*Omega_X = 0.139982) driven with the exact hi_class keywords from
fpab_variance_fixed_point.py. Requires the .venv classy = mochi MG wrapper.

Outputs (vs the frozen FPAB v0.3 benchmark table):
  sigma8(0), fsigma8(z);  P(k)/P_LCDM(k,0) -> mu_inf, k50, k90;  C_l^phiphi & C_l^TT
  ratios;  and confirmation that mochi's built-in ghost/gradient/tachyon tests PASS
  (the run only completes if the SEDE EFT is stable -- an independent check of the
  corpus's stability claim, through the actual Boltzmann code).
"""
import numpy as np
from classy import Class

import os
MOCHI = os.environ.get("MOCHI_DIR") or os.path.expanduser("~/Projects/mochi-class")
# corpus benchmark cosmology: h=0.68, Om=0.30, A_s=2.10e-9, n_s=0.965
COSMO = dict(h=0.68, omega_b=0.02237, omega_cdm=0.30*0.68**2 - 0.02237,
             A_s=2.10e-9, n_s=0.965, N_ur=3.046, output='mPk,mTk,tCl,pCl,lCl',
             lensing='yes', l_max_scalars=2500, P_k_max_1_Mpc=10.0)
Z_FS8 = [0.0, 0.3, 0.5, 1.0]

def run_lcdm():
    c=Class(); p=dict(COSMO); p['P_k_max_1/Mpc']=p.pop('P_k_max_1_Mpc')
    c.set(p); c.compute(); return c

def run_sede():
    c=Class(); p=dict(COSMO); p['P_k_max_1/Mpc']=p.pop('P_k_max_1_Mpc')
    p.update({
        'Omega_Lambda':0,'Omega_fld':0,'Omega_smg':-1,
        'gravity_model':'stable_params',
        'smg_file_name':'stable_params_input/fpab_sede_b020_stable.dat',
        'parameters_smg':'1.39982e-1',           # alpha_B(a=1)
        'method_gr_smg':'on','z_gr_smg':99,
        'expansion_model':'rho_de','expansion_smg':0.69991,
        'expansion_file_name':'stable_params_input/fpab_sede_b020_rho.dat',
        'output_background_smg':3,
        'skip_stability_tests_smg':'no','cs2_safe_smg':0,'D_safe_smg':0,
        'ct2_safe_smg':0,'M2_safe_smg':0,'a_min_stability_test_smg':0.01,
        'skip_math_stability_smg':'no','exp_rate_smg':1,
        'pert_initial_conditions_smg':'zero','method_qs_smg':'fully_dynamic',
    })
    c.set(p); c.compute(); return c

def pk_ratio(cs, cl, z, kh):
    h=COSMO['h']; k=kh*h
    return np.array([cs.pk(ki,z) for ki in k])/np.array([cl.pk(ki,z) for ki in k])

def fsigma8(c, z):
    # scale-independent proxy from classy; for MG the corpus quotes the ~large-scale value
    return c.scale_independent_growth_factor_f(z)*c.sigma8()*c.scale_independent_growth_factor(z)/c.scale_independent_growth_factor(0)

print("="*78); print("TIER-2 Stage-2 finite-k: SEDE FPAB through mochi-class MG Boltzmann"); print("="*78)
import os
if not os.path.isdir(os.path.join(MOCHI, "stable_params_input")):
    raise SystemExit("Set MOCHI_DIR to the built mochi-class fork (see src/README.md)")
os.makedirs(os.path.join(os.path.dirname(__file__),"output"), exist_ok=True); os.chdir(MOCHI)   # so relative stable_params_input/ paths resolve

print("  computing LCDM ..."); L=run_lcdm()
print("  computing SEDE (stable_params, alpha_B0=0.140) -- runs iff EFT is stable ...")
S=run_sede()
print("  >> SEDE run COMPLETED: mochi's ghost/gradient/tachyon/math stability tests PASSED")
print(f"     (independent confirmation of SEDE stability through the actual Boltzmann code)\n")

print(f"  {'quantity':28} {'LCDM':>10} {'SEDE':>10} {'FPAB v0.3 target':>18}")
print("  "+"-"*68)
s8L,s8S=L.sigma8(),S.sigma8()
print(f"  {'sigma8(0)':28} {s8L:10.5f} {s8S:10.5f} {'0.81088':>18}")
for z in Z_FS8:
    print(f"  {'fsigma8(z=%.1f)'%z:28} {fsigma8(L,z):10.5f} {fsigma8(S,z):10.5f} {('0.4222 @0' if z==0 else ''):>18}")

# --- P(k) enhancement (observable) ---
kh=np.logspace(-3.3, 0.7, 300)
rP=pk_ratio(S,L,0.0,kh)-1.0
aH0=COSMO['h']*100/299792.458                     # 1/Mpc (= H0/c)
print("  "+"-"*68)
print(f"  P/P_LCDM(k=0.1 h/Mpc, z=0)   = {np.interp(0.1,kh,rP)+1:.4f}   (target 1.0261)")
print(f"  P/P_LCDM plateau (k=5 h/Mpc) = {rP[-1]+1:.4f}   (integrated growth, < mu_inf)")

# --- mu(k): the FORCE response, from the modified Poisson for psi ---
# mu(k) proportional to (k^2 psi / delta_m); ratio to LCDM cancels normalizations.
def mu_of_k(cs, cl, z):
    tS,tL=cs.get_transfer(z),cl.get_transfer(z)
    kS=tS['k (h/Mpc)']; kL=tL['k (h/Mpc)']
    muS=np.abs(tS['psi'])/np.abs(tS['d_m'])
    muL=np.interp(kS, kL, np.abs(tL['psi'])/np.abs(tL['d_m']))
    return kS, muS/muL
kk,muk=mu_of_k(S,L,0.0)
muk=muk/muk[0]                                    # normalize to super-horizon (mu->1)
muinf_bg=float(S.get_background()['mu_inf'][-1])   # mochi's own high-k effective coupling at z=0
r=muk-1.0; plateau=r[-1]
def kcross(frac):
    tgt=frac*plateau; idx=np.where(r>=tgt)[0]
    return kk[idx[0]]*COSMO['h'] if len(idx) else np.nan   # 1/Mpc
k50,k90=kcross(0.5),kcross(0.9)
print(f"  mu_inf (mochi background, z=0)= {muinf_bg:.4f}   (target 1.05065)")
print(f"  mu(k) high-k plateau (transfer)= {muk[-1]:.4f}")
print(f"  k50 (mu half-rise)           = {k50:.2e}/Mpc = {k50/aH0:.2f} H0   (target 0.703 H0)")
print(f"  k90                          = {k90:.2e}/Mpc = {k90/aH0:.2f} H0   (target 2.14 H0)")

# lensing C_l^phiphi and TT ratios
clL,clS=L.lensed_cl(2000),S.lensed_cl(2000)
pL,pS=L.raw_cl(2000)['pp'],S.raw_cl(2000)['pp']
for ell in [10,100,1000]:
    print(f"  C_l^phiphi ratio (l={ell:4d})   = {pS[ell]/pL[ell]:.4f}"
          + ({100:'   (target 1.031)',1000:'  (target 1.021)'}.get(ell,'')))
ttL,ttS=clL['tt'],clS['tt']
print(f"  C_l^TT ratio (l=2, ISW)      = {ttS[2]/ttL[2]:.4f}   (target 0.923)")
# gamma_slip = phi/psi at high k (alpha_M=alpha_T=0 -> no slip)
tS=S.get_transfer(0.0); hik=tS['k (h/Mpc)']>1.0
slip=np.mean(np.abs(tS['phi'][hik])/np.abs(tS['psi'][hik]))
print(f"  gamma_slip = phi/psi (high-k) = {slip:.4f}   (target 1; alpha_M=alpha_T=0)")
print(f"\n  sigma8 ratio SEDE/LCDM       = {s8S/s8L:.4f}   (corpus 0.81088/0.8002 = 1.0133)")
print("\n  NOTE: benchmark cosmology (not marginalized). This is the finite-k prediction")
print("  the compressed Stage-1.5 fit could not produce -- now computed through the real")
print("  MG Boltzmann code. Full posterior still needs likelihood data (cobaya-install).")

# ---------------- figure ----------------
try:
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    OUT = os.path.join(os.path.dirname(__file__), "output", r"sede_stage2_finitek.png")
    fig,ax=plt.subplots(1,2,figsize=(11,4.2))
    ax[0].semilogx(kk,(muk-1)*100,'b-',lw=2,label=r'$\mu(k)-1$ (force response)')
    ax[0].semilogx(kh,rP*100,'g--',lw=1.8,label=r'$P/P_{\Lambda{\rm CDM}}-1$')
    ax[0].axhline((muinf_bg-1)*100,color='b',ls=':',lw=1,label=r'$\mu_\infty-1=%.1f\%%$'%((muinf_bg-1)*100))
    ax[0].set(xlabel=r'$k\ [h/{\rm Mpc}]$',ylabel=r'% deviation from $\Lambda$CDM (z=0)',
              title='SEDE finite-k force + power response'); ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)
    ells=np.arange(2,2000)
    ax[1].semilogx(ells,(pS[2:2000]/pL[2:2000]-1)*100,'m-',lw=1.8,label=r'$C_\ell^{\phi\phi}$ (lensing)')
    ax[1].semilogx(ells,(ttS[2:2000]/ttL[2:2000]-1)*100,'c-',lw=1.2,alpha=.8,label=r'$C_\ell^{TT}$')
    ax[1].axhline(0,color='k',lw=.5); ax[1].set_ylim(-10,10)
    ax[1].set(xlabel=r'multipole $\ell$',ylabel='% difference from $\\Lambda$CDM',
              title='SEDE CMB lensing + ISW'); ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
    fig.tight_layout(); fig.savefig(OUT,dpi=140); print(f"\n  figure -> {OUT}")
except Exception as e: print(f"  (figure skipped: {e})")
