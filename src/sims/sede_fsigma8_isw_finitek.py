#!/usr/bin/env python3
"""
Tier-1 #1: finite-k fsigma8(k,z) and the ISW signature of SEDE through mochi-class.
Upgrades D-6's quasi-static estimate (delta_DE/delta_m ~ 7e-5 at RSD) to the actual
scale-dependent growth from the MG Boltzmann code. Deliverables:
  - fsigma8(k,z), SEDE/LCDM ratio, at RSD (k~0.1) vs large scales (k~0.01) vs near-horizon
  - the low-l ISW C_l^TT ratio (where the gate-matter coupling shows up)
The claim under test (D-6): the modification is horizon-confined -> RSD fsigma8 barely moves,
ISW/low-l moves at the few-% level.
"""
import numpy as np, os
from classy import Class
_MOCHI = os.environ.get("MOCHI_DIR") or os.path.expanduser("~/Projects/mochi-class")
if not os.path.isdir(os.path.join(_MOCHI, "stable_params_input")):
    raise SystemExit("Set MOCHI_DIR to the built mochi-class fork (see src/README.md)")
os.chdir(_MOCHI)   # so relative stable_params_input/ paths resolve

COS=dict(h=0.68,omega_b=0.02237,omega_cdm=0.30*0.68**2-0.02237,A_s=2.10e-9,n_s=0.965,N_ur=3.046,
         output='mPk,tCl,lCl',lensing='yes',l_max_scalars=1500)
ZP="0,0.1,0.2,0.3,0.4,0.5,0.7,0.9,1.0,1.1"        # fine z grid for df/dz
def run(sede):
    c=Class(); p=dict(COS); p['P_k_max_1/Mpc']=2.0; p['z_pk']=ZP
    if sede:
        p.update({'Omega_Lambda':0,'Omega_fld':0,'Omega_smg':-1,'gravity_model':'stable_params',
          'smg_file_name':'stable_params_input/fpab_sede_b020_stable.dat','parameters_smg':'1.39982e-1',
          'method_gr_smg':'on','z_gr_smg':5,'expansion_model':'rho_de','expansion_smg':0.69991,  # z_gr_smg=5 not 99: 99 excites a spurious k~0.02 scalar oscillation -> fake ISW (see TIER2_head_to_head_result.md)
          'expansion_file_name':'stable_params_input/fpab_sede_b020_rho.dat','output_background_smg':3,
          'skip_stability_tests_smg':'no','method_qs_smg':'fully_dynamic','pert_initial_conditions_smg':'zero'})
    c.set(p); c.compute(); return c

def fsigma8_k(c, kh, zc):
    """fsigma8(k,z): f(k,z)=-(1+z)/2 d lnP/dz ; sigma8(z) from code; return f*sigma8*sqrt(P/P8-normalized)."""
    h=COS['h']; k=kh*h
    zs=np.array([0.0,0.1,0.2,0.3,0.4,0.5,0.7,0.9,1.0,1.1])
    lnP=np.array([np.log(c.pk(k,z)) for z in zs])
    from scipy.interpolate import CubicSpline
    dP=CubicSpline(zs,lnP)
    f=-(1+zc)/2*dP(zc,1)                              # f(k,z)=dlnD/dlna, D=sqrt(P)
    s8z=c.sigma8()*np.sqrt(c.pk(k,zc)/c.pk(k,0.0))    # sigma8(z) scaled by the k-mode growth (approx tracer)
    return f*s8z

L=run(False); S=run(True)
print("="*74); print("Tier-1 #1: SEDE finite-k fsigma8 + ISW through mochi-class"); print("="*74)
print(f"  sigma8(0): LCDM={L.sigma8():.4f}  SEDE={S.sigma8():.4f}")
print()
print("  fsigma8(k,z) SEDE/LCDM - 1  [%]   (the RSD-scale deviation the modification produces):")
print(f"  {'k[h/Mpc]':>9} | " + " | ".join(f"z={z}" for z in [0.0,0.3,0.5,1.0]))
print("  "+"-"*54)
for kh in [0.005,0.01,0.05,0.1,0.2]:
    row=[]
    for z in [0.0,0.3,0.5,1.0]:
        r=fsigma8_k(S,kh,z)/fsigma8_k(L,kh,z)-1
        row.append(f"{r*100:+6.2f}")
    tag="  <-- RSD scale" if kh==0.1 else ("  <-- ~horizon" if kh<=0.01 else "")
    print(f"  {kh:9.3f} | " + " | ".join(row) + tag)
print()
# ISW: low-l TT ratio
clL=L.raw_cl(1200)['tt']; clS=S.raw_cl(1200)['tt']
print("  ISW signature -- C_l^TT ratio SEDE/LCDM (low-l = ISW-dominated):")
for l in [2,5,10,20,50,200]:
    print(f"    l={l:4d}: {clS[l]/clL[l]:.4f}")
print()
print("="*74); print("  VERDICT (test of the D-6 horizon-confinement claim)"); print("="*74)
r_rsd=fsigma8_k(S,0.1,0.5)/fsigma8_k(L,0.1,0.5)-1
r_hor=fsigma8_k(S,0.005,0.5)/fsigma8_k(L,0.005,0.5)-1
print(f"  fsigma8 deviation at RSD (k=0.1, z=0.5)  = {r_rsd*100:+.2f}%")
print(f"  fsigma8 deviation near horizon (k=0.005) = {r_hor*100:+.2f}%")
print(f"  ISW (l=2) deviation                      = {(clS[2]/clL[2]-1)*100:+.1f}%")
print("  -> if RSD deviation << horizon/ISW deviation, the modification is horizon-confined,")
print("     confirming D-6: the gate-matter coupling shows up in ISW/low-l, not in RSD fsigma8.")
