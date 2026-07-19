#!/usr/bin/env python3
"""Rigorous phiphi test (obstest 'decisive next test').

Step 1: for each model, find the CMB-PRIMARY (lowl TT + high-l TTTEEE-lite) best
        fit by Nelder-Mead over {H0, omega_b, omega_cdm, A_s, n_s} (tau, A_planck
        fixed as in the paper). -> Delta chi2_primary (should reproduce ~ -3.3).
Step 2: at each primary best fit, ADD the Planck 2018 lensing (phiphi, native,
        CMB-marginalized) likelihood. -> Delta chi2_phiphi and the augmented
        Delta chi2_primary+phiphi.
Step 3: re-minimize INCLUDING phiphi -> Delta chi2_joint(with phiphi): how much of
        the primary preference survives once the direct lensing datum is in the fit.

Interpretation: if adding phiphi drives Delta chi2 toward 0 / positive, the primary
preference was (partly) an absorbed A_lens systematic; if it stays ~ -3, it is physical.

Run from src/mcmc/ with MOCHI_DIR set.
"""
import os, sys, time
os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np
from scipy.optimize import minimize
from cobaya.model import get_model

PACKAGES = os.environ.get("COBAYA_PACKAGES",
                          os.path.join(os.path.dirname(os.path.abspath(__file__)), "packages"))
TAU, APLANCK = 0.0544, 1.0025

PRIMARY = {"planck_2018_lowl.TT": None,
           "planck_2018_highl_plik.TTTEEE_lite": None}
PHIPHI = {"planck_2018_lensing.native": None}

SEDE_THEORY = {"classy": {"path": "global", "extra_args": {
    "N_ur": 3.046, "non_linear": "halofit",
    "Omega_Lambda": 0, "Omega_fld": 0, "Omega_smg": -1,
    "gravity_model": "stable_params",
    "smg_file_name": "inputs/fpab_sede_b020_stable.dat",
    "parameters_smg": "1.39982e-1", "method_gr_smg": "on", "z_gr_smg": 5,
    "expansion_model": "rho_de", "expansion_smg": 0.69991,
    "expansion_file_name": "inputs/fpab_sede_b020_rho.dat",
    "output_background_smg": 3, "skip_stability_tests_smg": "no",
    "method_qs_smg": "fully_dynamic", "pert_initial_conditions_smg": "zero"}}}
LCDM_THEORY = {"classy": {"path": "global",
               "extra_args": {"N_ur": 3.046, "non_linear": "halofit"}}}

# order: H0, omega_b, omega_cdm, A_s(1e-9), n_s
X0 = {"LCDM": np.array([68.53, 0.02250, 0.1185, 2.10, 0.9682]),
      "SEDE": np.array([68.86, 0.02250, 0.1193, 2.10, 0.9664])}
BOUNDS = [(60, 80), (0.019, 0.025), (0.10, 0.14), (1.6, 2.6), (0.94, 0.99)]
SIMPLEX = np.array([0.3, 0.0003, 0.001, 0.03, 0.003])   # initial step per param


def build(theory, likes):
    info = {"theory": theory, "likelihood": likes, "packages_path": PACKAGES,
            "params": {
                "H0": {"prior": {"min": 60, "max": 80}},
                "omega_b": {"prior": {"min": 0.019, "max": 0.025}},
                "omega_cdm": {"prior": {"min": 0.10, "max": 0.14}},
                "A_s": {"prior": {"min": 1.0e-9, "max": 3.0e-9}},
                "n_s": {"prior": {"min": 0.94, "max": 0.99}},
                "tau_reio": TAU, "A_planck": APLANCK}}
    return get_model(info)


def chi2_vec(model, x):
    pt = {"H0": x[0], "omega_b": x[1], "omega_cdm": x[2],
          "A_s": x[3] * 1e-9, "n_s": x[4]}
    for b, (lo, hi) in zip(x, BOUNDS):
        if not (lo <= b <= hi):
            return 1e10
    try:
        loglikes, _ = model.loglikes(pt, as_dict=True)
    except Exception:
        return 1e10
    return -2.0 * sum(loglikes.values())


def minimize_model(model, x0, label):
    best = None
    for s, jitter in enumerate([0.0, 1.0]):
        start = x0 + jitter * SIMPLEX * np.array([1, -1, 1, -1, 1]) * 0.5
        isim = np.vstack([start] + [start + np.eye(5)[i] * SIMPLEX[i] for i in range(5)])
        r = minimize(lambda x: chi2_vec(model, x), start, method="Nelder-Mead",
                     options={"initial_simplex": isim, "xatol": 1e-4,
                              "fatol": 1e-3, "maxfev": 700})
        print(f"    [{label} start {s}] chi2_min = {r.fun:.3f} at "
              f"H0={r.x[0]:.3f} ob={r.x[1]:.5f} oc={r.x[2]:.5f} "
              f"As={r.x[3]:.4f}e-9 ns={r.x[4]:.5f}  (nfev {r.nfev})", flush=True)
        if best is None or r.fun < best.fun:
            best = r
    return best


def eval_at(theory, likes, x):
    m = build(theory, likes)
    pt = {"H0": x[0], "omega_b": x[1], "omega_cdm": x[2],
          "A_s": x[3] * 1e-9, "n_s": x[4]}
    loglikes, _ = m.loglikes(pt, as_dict=True)
    return {k: -2.0 * v for k, v in loglikes.items()}


def main():
    t0 = time.time()
    THEORY = {"LCDM": LCDM_THEORY, "SEDE": SEDE_THEORY}
    bf, primary = {}, {}
    print("=" * 72 + "\nSTEP 1: CMB-primary (lowlTT + TTTEEE-lite) best fit\n" + "=" * 72)
    for tag in ("LCDM", "SEDE"):
        print(f"\n>>> minimizing {tag} (primary CMB)")
        m = build(THEORY[tag], PRIMARY)
        r = minimize_model(m, X0[tag], tag)
        bf[tag] = r.x
        primary[tag] = r.fun
    dprim = primary["SEDE"] - primary["LCDM"]
    print(f"\nDelta chi2_primary (SEDE-LCDM) = {dprim:+.3f}  "
          f"[paper best-fit TTTEEE term -3.3]")

    print("\n" + "=" * 72 + "\nSTEP 2: add phiphi AT the primary best fit\n" + "=" * 72)
    phi = {}
    for tag in ("LCDM", "SEDE"):
        c = eval_at(THEORY[tag], {**PRIMARY, **PHIPHI}, bf[tag])
        phi[tag] = next(v for k, v in c.items() if "lensing" in k)
        print(f"    {tag}: chi2_phiphi(at primary bf) = {phi[tag]:.3f}   "
              f"[full: {', '.join(f'{k.split(chr(46))[-1]}={v:.2f}' for k,v in c.items())}]")
    dphi = phi["SEDE"] - phi["LCDM"]
    print(f"\nDelta chi2_phiphi (SEDE-LCDM)          = {dphi:+.3f}")
    print(f"Delta chi2_primary                      = {dprim:+.3f}")
    print(f"Delta chi2_primary+phiphi (augmented)   = {dprim + dphi:+.3f}")

    print("\n" + "=" * 72 + "\nSTEP 3: re-minimize INCLUDING phiphi\n" + "=" * 72)
    joint = {}
    for tag in ("LCDM", "SEDE"):
        print(f"\n>>> minimizing {tag} (primary + phiphi)")
        m = build(THEORY[tag], {**PRIMARY, **PHIPHI})
        r = minimize_model(m, bf[tag], tag)
        joint[tag] = r.fun
    djoint = joint["SEDE"] - joint["LCDM"]
    print(f"\nDelta chi2_joint(with phiphi) (SEDE-LCDM) = {djoint:+.3f}")

    print("\n" + "=" * 72 + "\nSUMMARY\n" + "=" * 72)
    print(f"  primary CMB only                : Delta chi2 = {dprim:+.3f}")
    print(f"  + phiphi at primary best fit    : Delta chi2 = {dprim + dphi:+.3f}  "
          f"(phiphi penalty {dphi:+.3f})")
    print(f"  refit with phiphi in the fit    : Delta chi2 = {djoint:+.3f}")
    verdict = ("PHYSICAL (phiphi supports the extra lensing)" if dphi <= 0.5 else
               "SYSTEMATIC-LEANING (phiphi erodes the primary preference)")
    print(f"  verdict: {verdict}")
    print(f"  [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
