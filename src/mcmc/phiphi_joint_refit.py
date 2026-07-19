#!/usr/bin/env python3
"""DEFINITIVE phiphi test (the fix): a BAO+SN+BBN-anchored joint refit, without and
with the Planck lensing (phiphi) likelihood.

This removes the two defects of the earlier estimates:
  - A_s-profiled (-1.9): a chi2 SUM evaluated at the primary optimum (an upper bound, not a fit).
  - full CMB-only refit (-2.2): a valid joint minimum but with omega_cdm drifting to 0.120,
    unpinned by BAO/SN.
Here BAO DR2 + Pantheon+ + BBN pin the background AND both models are properly refit, so the
difference of joint minima is the honest post-phiphi preference.

For each model, Powell-minimize over {H0, omega_b, omega_cdm, A_s, n_s} (tau=0.0544, A_planck=1.0025
fixed; BBN Gaussian prior on omega_b as in the paper) in two cases:
  (J)   BAO + SN + lowlTT + TTTEEE-lite            -> Delta chi2_joint      (reproduce paper -3.4)
  (JL)  + planck_2018_lensing.native (phiphi)      -> Delta chi2_joint+phiphi (THE number)

Run from src/mcmc/ with MOCHI_DIR set and
COBAYA_PACKAGES=<your cobaya packages dir with BAO+SN+Planck primary+lensing>
Writes phiphi_joint_refit_result.json.
"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np
from scipy.optimize import minimize
from cobaya.model import get_model

PACKAGES = os.environ.get("COBAYA_PACKAGES", os.path.join(os.path.dirname(os.path.abspath(__file__)), "packages"))
BBN_MU, BBN_SIG = 0.02218, 0.00055

JOINT = {"bao.desi_dr2.desi_bao_all": None, "sn.pantheonplus": None,
         "planck_2018_lowl.TT": None, "planck_2018_highl_plik.TTTEEE_lite": None}
JL = {**JOINT, "planck_2018_lensing.native": None}

SEDE_THEORY = {"classy": {"path": "global", "extra_args": {
    "N_ur": 3.046, "non_linear": "halofit", "Omega_Lambda": 0, "Omega_fld": 0,
    "Omega_smg": -1, "gravity_model": "stable_params",
    "smg_file_name": "inputs/fpab_sede_b020_stable.dat", "parameters_smg": "1.39982e-1",
    "method_gr_smg": "on", "z_gr_smg": 5, "expansion_model": "rho_de",
    "expansion_smg": 0.69991, "expansion_file_name": "inputs/fpab_sede_b020_rho.dat",
    "output_background_smg": 3, "skip_stability_tests_smg": "no",
    "method_qs_smg": "fully_dynamic", "pert_initial_conditions_smg": "zero"}}}
LCDM_THEORY = {"classy": {"path": "global", "extra_args": {"N_ur": 3.046, "non_linear": "halofit"}}}
THEORY = {"LCDM": LCDM_THEORY, "SEDE": SEDE_THEORY}

X0 = {"LCDM": np.array([68.53, 0.02237, 0.1185, 2.101, 0.9682]),
      "SEDE": np.array([68.86, 0.02237, 0.1193, 2.105, 0.9664])}
BOUNDS = [(60, 80), (0.019, 0.025), (0.10, 0.14), (1.6, 2.6), (0.94, 0.99)]


def build(theory, likes):
    return get_model({"theory": theory, "likelihood": likes, "packages_path": PACKAGES,
        "params": {"H0": {"prior": {"min": 60, "max": 80}},
                   "omega_b": {"prior": {"min": 0.019, "max": 0.025}},
                   "omega_cdm": {"prior": {"min": 0.10, "max": 0.14}},
                   "A_s": {"prior": {"min": 1.0e-9, "max": 3.0e-9}},
                   "n_s": {"prior": {"min": 0.94, "max": 0.99}},
                   "tau_reio": 0.0544, "A_planck": 1.0025}})


def chi2_parts(model, x):
    pt = {"H0": x[0], "omega_b": x[1], "omega_cdm": x[2], "A_s": x[3] * 1e-9, "n_s": x[4]}
    ll, _ = model.loglikes(pt, as_dict=True)
    d = {k: -2.0 * v for k, v in ll.items()}
    d["bbn_prior"] = ((x[1] - BBN_MU) / BBN_SIG) ** 2   # Gaussian BBN penalty on omega_b
    return d


def obj(model, x):
    for v, (lo, hi) in zip(x, BOUNDS):
        if not (lo <= v <= hi):
            return 1e8
    try:
        return sum(chi2_parts(model, x).values())
    except Exception:
        return 1e8


def refit(tag, likes, x0, label):
    m = build(THEORY[tag], likes)
    t = time.time()
    r = minimize(lambda x: obj(m, x), x0, method="Powell", bounds=BOUNDS,
                 options={"xtol": 1e-4, "ftol": 1e-4, "maxfev": 800})
    parts = chi2_parts(m, r.x)
    print(f"  [{label} {tag}] chi2={r.fun:.3f} x=[H0={r.x[0]:.3f} ob={r.x[1]:.5f} "
          f"oc={r.x[2]:.5f} As={r.x[3]:.4f} ns={r.x[4]:.5f}] nfev={r.nfev} ({time.time()-t:.0f}s)",
          flush=True)
    return r.fun, r.x.tolist(), parts


def selftest():
    print("=== SELF-TEST: one joint+phiphi eval per model (validate likelihoods load) ===", flush=True)
    ok = True
    for tag in ("LCDM", "SEDE"):
        try:
            m = build(THEORY[tag], JL)
            p = chi2_parts(m, X0[tag])
            print(f"  {tag}: " + "  ".join(f"{k.split('.')[-1]}={v:.2f}" for k, v in p.items())
                  + f"  TOT={sum(p.values()):.2f}", flush=True)
        except Exception as e:
            ok = False; print(f"  {tag} FAILED: {e}", flush=True)
    return ok


def main():
    t0 = time.time()
    if not selftest():
        print("self-test failed; aborting"); sys.exit(1)
    out = {}
    print("\n=== JOINT REFIT: (J) BAO+SN+lowlTT+TTTEEE  and  (JL) +phiphi ===", flush=True)
    for tag in ("LCDM", "SEDE"):
        fj, xj, cj = refit(tag, JOINT, X0[tag], "J ")
        fl, xl, cl = refit(tag, JL, xj, "JL")
        out[tag] = {"joint": {"chi2": fj, "x": xj, "parts": cj},
                    "joint_phiphi": {"chi2": fl, "x": xl, "parts": cl}}
    dj = out["SEDE"]["joint"]["chi2"] - out["LCDM"]["joint"]["chi2"]
    djl = out["SEDE"]["joint_phiphi"]["chi2"] - out["LCDM"]["joint_phiphi"]["chi2"]
    phi = {t: next(v for k, v in out[t]["joint_phiphi"]["parts"].items() if "lensing" in k) for t in out}
    print("\n=== SUMMARY (BAO/SN/BBN-anchored joint) ===")
    print(f"  Delta chi2_joint  (no phiphi)   = {dj:+.3f}   [paper best-fit -3.4]")
    print(f"  Delta chi2_joint  (with phiphi) = {djl:+.3f}   <-- THE post-phiphi preference")
    print(f"  phiphi penalty (SEDE-LCDM) at joint bf = {phi['SEDE']-phi['LCDM']:+.3f}")
    print(f"  phiphi erosion  = {dj:+.3f} -> {djl:+.3f}")
    out["summary"] = {"delta_joint_no_phiphi": dj, "delta_joint_with_phiphi": djl,
                      "phiphi_penalty_at_joint": phi["SEDE"] - phi["LCDM"], "seconds": time.time() - t0}
    with open("phiphi_joint_refit_result.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"  wrote phiphi_joint_refit_result.json  [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
