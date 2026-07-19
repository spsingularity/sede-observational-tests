#!/usr/bin/env python3
"""The decisive next test (obstest §7/§8): does the Planck lensing (phiphi)
likelihood support SEDE's extra TTTEEE smoothing (physical lensing) or penalize
it (an absorbed A_lens systematic)?

We evaluate chi2 from lowl TT + high-l TTTEEE-lite + Planck 2018 lensing (native,
CMB-marginalized, whose lensing amplitude is ~1.0) at each model's production
posterior-mean point, and compare chi2_phiphi(SEDE) vs chi2_phiphi(LCDM).

  Delta chi2_phiphi > 0  => phiphi disfavours SEDE's boosted C_l^phiphi
                            => the -3.3 TTTEEE gain is (partly) an absorbed A_lens systematic
  Delta chi2_phiphi <~ 0 => phiphi is consistent with SEDE's extra lensing
                            => the gain is physical

Run from src/mcmc/ (inputs/*.dat are referenced relatively).
"""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "4")
from cobaya.model import get_model

PACKAGES = os.environ.get("COBAYA_PACKAGES",
                          os.path.join(os.path.dirname(os.path.abspath(__file__)), "packages"))

LIKES = {
    "planck_2018_lowl.TT": None,
    "planck_2018_highl_plik.TTTEEE_lite": None,
    "planck_2018_lensing.native": None,
}

SEDE_THEORY = {
    "classy": {"path": "global", "extra_args": {
        "N_ur": 3.046, "non_linear": "halofit",
        "Omega_Lambda": 0, "Omega_fld": 0, "Omega_smg": -1,
        "gravity_model": "stable_params",
        "smg_file_name": "inputs/fpab_sede_b020_stable.dat",
        "parameters_smg": "1.39982e-1", "method_gr_smg": "on", "z_gr_smg": 5,
        "expansion_model": "rho_de", "expansion_smg": 0.69991,
        "expansion_file_name": "inputs/fpab_sede_b020_rho.dat",
        "output_background_smg": 3, "skip_stability_tests_smg": "no",
        "method_qs_smg": "fully_dynamic", "pert_initial_conditions_smg": "zero",
    }}}
LCDM_THEORY = {"classy": {"path": "global",
                          "extra_args": {"N_ur": 3.046, "non_linear": "halofit"}}}

# CMB-primary BEST-FIT points (paper sec 5 table: H0, omega_cdm, n_s), with the
# minimization's fixed logA=3.045 (A_s=1e-10 exp(3.045)=2.101e-9), omega_b box-edge
# 0.0225, tau=0.0544 -- the settings under which we reproduced TTTEEE=582.7.
import numpy as np
POINTS = {
    "LCDM": dict(H0=68.53, omega_b=0.0225, omega_cdm=0.1185, A_s=2.101e-9,
                 n_s=0.9682, tau_reio=0.0544),
    "SEDE": dict(H0=68.86, omega_b=0.0225, omega_cdm=0.1193, A_s=2.101e-9,
                 n_s=0.9664, tau_reio=0.0544),
}

def run(tag, theory):
    info = {
        "theory": theory, "likelihood": LIKES, "packages_path": PACKAGES,
        "params": {
            "H0": {"prior": {"min": 60, "max": 80}},
            "omega_b": {"prior": {"min": 0.019, "max": 0.025}},
            "omega_cdm": {"prior": {"min": 0.1, "max": 0.14}},
            "A_s": {"prior": {"min": 1.0e-9, "max": 3.0e-9}},
            "n_s": {"prior": {"min": 0.94, "max": 0.99}},
            "tau_reio": {"prior": {"min": 0.02, "max": 0.1}},
            "A_planck": 1.0025,   # plik_lite calibration nuisance, held fixed (paper sec 2)
        },
    }
    model = get_model(info)
    pt = POINTS[tag]
    loglikes, _ = model.loglikes(pt, as_dict=True)
    chi2 = {k: -2.0 * v for k, v in loglikes.items()}
    return chi2

if __name__ == "__main__":
    print("=" * 70)
    res = {}
    for tag, theory in (("LCDM", LCDM_THEORY), ("SEDE", SEDE_THEORY)):
        print(f"\n>>> {tag} at {POINTS[tag]}")
        res[tag] = run(tag, theory)
        for k, v in res[tag].items():
            print(f"    chi2[{k}] = {v:.3f}")
        print(f"    chi2_total = {sum(res[tag].values()):.3f}")
    print("\n" + "=" * 70)
    def g(tag, key):
        return next(v for k, v in res[tag].items() if key in k)
    dphi = g("SEDE", "lensing") - g("LCDM", "lensing")
    dtt  = (g("SEDE","lowl")+g("SEDE","TTTEEE")) - (g("LCDM","lowl")+g("LCDM","TTTEEE"))
    print(f"Delta chi2_phiphi  (SEDE-LCDM) = {dphi:+.3f}")
    print(f"Delta chi2_TT+TTTEEE           = {dtt:+.3f}")
    print(f"Delta chi2_TT+TTTEEE+phiphi    = {dtt+dphi:+.3f}")
    print("interpretation: dphi>0 => phiphi penalizes SEDE (absorbed A_lens);"
          " dphi<~0 => physical.")
