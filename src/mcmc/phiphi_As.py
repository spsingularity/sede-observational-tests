#!/usr/bin/env python3
"""phiphi test, A_s-profiled. Fix (H0, omega_b, omega_cdm, n_s, tau) at each model's
published sec-5 best fit; minimize A_s alone over the CMB primary (lowlTT + TTTEEE-lite)
so both models sit at their damping-tail optimum (reproducing TTTEEE ~583 / ~580, the
paper's -3.3). Then evaluate Planck 2018 lensing (phiphi, native) at that optimum.

Run from src/mcmc/ with MOCHI_DIR set."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np
from scipy.optimize import minimize_scalar
from cobaya.model import get_model

PACKAGES = os.environ.get("COBAYA_PACKAGES", os.path.join(os.path.dirname(os.path.abspath(__file__)), "packages"))
PRIMARY = {"planck_2018_lowl.TT": None, "planck_2018_highl_plik.TTTEEE_lite": None}
ALL = {**PRIMARY, "planck_2018_lensing.native": None}

SEDE_THEORY = {"classy": {"path": "global", "extra_args": {
    "N_ur": 3.046, "non_linear": "halofit", "Omega_Lambda": 0, "Omega_fld": 0,
    "Omega_smg": -1, "gravity_model": "stable_params",
    "smg_file_name": "inputs/fpab_sede_b020_stable.dat", "parameters_smg": "1.39982e-1",
    "method_gr_smg": "on", "z_gr_smg": 5, "expansion_model": "rho_de",
    "expansion_smg": 0.69991, "expansion_file_name": "inputs/fpab_sede_b020_rho.dat",
    "output_background_smg": 3, "skip_stability_tests_smg": "no",
    "method_qs_smg": "fully_dynamic", "pert_initial_conditions_smg": "zero"}}}
LCDM_THEORY = {"classy": {"path": "global", "extra_args": {"N_ur": 3.046, "non_linear": "halofit"}}}

FIXED = {  # H0, omega_b, omega_cdm, n_s (tau=0.0544, A_planck=1.0025)
    "LCDM": dict(H0=68.53, omega_b=0.0225, omega_cdm=0.1185, n_s=0.9682),
    "SEDE": dict(H0=68.86, omega_b=0.0225, omega_cdm=0.1193, n_s=0.9664)}


def build(theory, likes):
    return get_model({"theory": theory, "likelihood": likes, "packages_path": PACKAGES,
        "params": {"H0": {"prior": {"min": 60, "max": 80}},
                   "omega_b": {"prior": {"min": 0.019, "max": 0.025}},
                   "omega_cdm": {"prior": {"min": 0.10, "max": 0.14}},
                   "A_s": {"prior": {"min": 1.0e-9, "max": 3.0e-9}},
                   "n_s": {"prior": {"min": 0.94, "max": 0.99}},
                   "tau_reio": 0.0544, "A_planck": 1.0025}})


def chi2(model, fixed, As9):
    pt = {**fixed, "A_s": As9 * 1e-9}
    ll, _ = model.loglikes(pt, as_dict=True)
    return {k: -2.0 * v for k, v in ll.items()}


res = {}
for tag, theory in (("LCDM", LCDM_THEORY), ("SEDE", SEDE_THEORY)):
    mp = build(theory, PRIMARY)
    f = FIXED[tag]
    obj = lambda As9: sum(chi2(mp, f, As9).values())
    r = minimize_scalar(obj, bounds=(1.9, 2.3), method="bounded",
                        options={"xatol": 2e-3})
    As_bf = r.x
    ma = build(theory, ALL)
    c = chi2(ma, f, As_bf)
    res[tag] = dict(As=As_bf, **c)
    print(f">>> {tag}: A_s_bf={As_bf:.4f}e-9  "
          f"lowlTT={c['planck_2018_lowl.TT']:.3f}  "
          f"TTTEEE={c['planck_2018_highl_plik.TTTEEE_lite']:.3f}  "
          f"phiphi={c['planck_2018_lensing.native']:.3f}", flush=True)

print("=" * 68)
g = lambda t, key: next(v for k, v in res[t].items() if key in k)
dprim = (g("SEDE", "lowl") + g("SEDE", "TTTEEE")) - (g("LCDM", "lowl") + g("LCDM", "TTTEEE"))
dphi = g("SEDE", "lensing") - g("LCDM", "lensing")
print(f"Delta chi2_primary (lowlTT+TTTEEE)      = {dprim:+.3f}   [paper -3.3]")
print(f"Delta chi2_phiphi                        = {dphi:+.3f}")
print(f"Delta chi2_primary+phiphi (augmented)    = {dprim + dphi:+.3f}")
print(f"phiphi survival: {100*(dprim+dphi)/dprim:.0f}% of the primary preference remains"
      if dprim < 0 else "n/a")
