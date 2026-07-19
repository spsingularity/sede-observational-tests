#!/usr/bin/env python3
"""phiphi test, FULL multi-parameter refit (pins the A_s-profiled -1.9 exactly).

For each model, Powell-minimize over {H0, omega_b, omega_cdm, A_s, n_s} (tau=0.0544,
A_planck=1.0025 fixed) in two cases:
  (P)  CMB primary only  : lowl TT + TTTEEE-lite          -> Delta chi2_primary (full refit)
  (PL) primary + phiphi  : + planck_2018_lensing.native   -> Delta chi2_joint (full refit)

Case P should reproduce the paper's -3.3. Case PL lets ALL parameters readjust to partly
re-absorb the phiphi penalty, giving the true post-phiphi preference (A_s-profiled gave -1.9).

Start from the phiphi_As best fits (already near-optimal). Powell needs ~100-200 evals/case.
Run from src/mcmc/ with MOCHI_DIR set. Writes phiphi_fullrefit_result.json.
"""
import os, json, time
os.environ.setdefault("OMP_NUM_THREADS", "4")
import numpy as np
from scipy.optimize import minimize
from cobaya.model import get_model

PACKAGES = os.environ.get("COBAYA_PACKAGES", os.path.join(os.path.dirname(os.path.abspath(__file__)), "packages"))
PRIMARY = {"planck_2018_lowl.TT": None, "planck_2018_highl_plik.TTTEEE_lite": None}
PL = {**PRIMARY, "planck_2018_lensing.native": None}

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

# start = phiphi_As best fits: H0, omega_b, omega_cdm, A_s(1e-9), n_s
X0 = {"LCDM": np.array([68.53, 0.0225, 0.1185, 2.1007, 0.9682]),
      "SEDE": np.array([68.86, 0.0225, 0.1193, 2.1054, 0.9664])}
BOUNDS = [(60, 80), (0.019, 0.025), (0.10, 0.14), (1.6, 2.6), (0.94, 0.99)]


def build(theory, likes):
    return get_model({"theory": theory, "likelihood": likes, "packages_path": PACKAGES,
        "params": {"H0": {"prior": {"min": 60, "max": 80}},
                   "omega_b": {"prior": {"min": 0.019, "max": 0.025}},
                   "omega_cdm": {"prior": {"min": 0.10, "max": 0.14}},
                   "A_s": {"prior": {"min": 1.0e-9, "max": 3.0e-9}},
                   "n_s": {"prior": {"min": 0.94, "max": 0.99}},
                   "tau_reio": 0.0544, "A_planck": 1.0025}})


def chi2_dict(model, x):
    pt = {"H0": x[0], "omega_b": x[1], "omega_cdm": x[2], "A_s": x[3] * 1e-9, "n_s": x[4]}
    ll, _ = model.loglikes(pt, as_dict=True)
    return {k: -2.0 * v for k, v in ll.items()}


def obj(model, x):
    for v, (lo, hi) in zip(x, BOUNDS):
        if not (lo <= v <= hi):
            return 1e8
    try:
        return sum(chi2_dict(model, x).values())
    except Exception:
        return 1e8


def refit(tag, likes, x0, label):
    m = build(THEORY[tag], likes)
    t = time.time()
    r = minimize(lambda x: obj(m, x), x0, method="Powell",
                 bounds=BOUNDS, options={"xtol": 1e-4, "ftol": 1e-4, "maxfev": 600})
    c = chi2_dict(m, r.x)
    print(f"  [{label} {tag}] chi2={r.fun:.3f}  x=[H0={r.x[0]:.3f} ob={r.x[1]:.5f} "
          f"oc={r.x[2]:.5f} As={r.x[3]:.4f} ns={r.x[4]:.5f}] nfev={r.nfev} "
          f"({time.time()-t:.0f}s)", flush=True)
    return r.fun, r.x, c


def main():
    t0 = time.time()
    out = {}
    print("=== FULL REFIT: primary (P) and primary+phiphi (PL) ===", flush=True)
    for tag in ("LCDM", "SEDE"):
        fp, xp, cp = refit(tag, PRIMARY, X0[tag], "P ")
        fl, xl, cl = refit(tag, PL, xp, "PL")   # start PL from the primary best fit
        out[tag] = {"primary": {"chi2": fp, "x": xp.tolist(), "parts": cp},
                    "with_phiphi": {"chi2": fl, "x": xl.tolist(), "parts": cl}}
    g = lambda t, case: out[t][case]["parts"]
    dprim = (sum(g("SEDE", "primary").values()) - sum(g("LCDM", "primary").values()))
    # joint chi2 for case PL already includes phiphi:
    djoint = out["SEDE"]["with_phiphi"]["chi2"] - out["LCDM"]["with_phiphi"]["chi2"]
    phi_pl = {t: next(v for k, v in g(t, "with_phiphi").items() if "lensing" in k) for t in out}
    dphi = phi_pl["SEDE"] - phi_pl["LCDM"]
    print("\n=== SUMMARY (full refit) ===")
    print(f"  Delta chi2_primary   (full refit)          = {dprim:+.3f}   [paper -3.3]")
    print(f"  Delta chi2_joint     (primary+phiphi refit) = {djoint:+.3f}   [A_s-profiled gave -1.9]")
    print(f"  phiphi term at the joint best fit (SEDE-LCDM) = {dphi:+.3f}")
    out["summary"] = {"delta_primary_fullrefit": dprim,
                      "delta_joint_with_phiphi_fullrefit": djoint,
                      "delta_phiphi_at_joint": dphi, "seconds": time.time() - t0}
    with open("phiphi_fullrefit_result.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"  wrote phiphi_fullrefit_result.json  [{time.time()-t0:.0f}s]")


if __name__ == "__main__":
    main()
