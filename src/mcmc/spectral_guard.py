"""Per-sample spectral sanity guard for SEDE/smg runs (cobaya external likelihood).

The FPAB stable_params reconstruction has an early-time/low-k tachyonic window in its
scalar sector; with z_gr_smg <= 5 the spectra are converged and sane across the explored
parameter box (see TIER2_head_to_head_result.md), but a production MCMC must not silently
accept a sample whose CMB transfer has gone pathological (the failure mode is a fake ISW:
low-l C_l^TT orders of magnitude too large, and/or negative C_l at high l).

This likelihood returns 0 for a sane spectrum and -inf for a pathological one, so bad
samples are rejected explicitly rather than by a chi2 wall. Add to any SEDE yaml exactly as
the shipped lcdm_mc.yaml / sede_mc.yaml do (run cobaya from src/mcmc/ so python_path: . resolves):

  likelihood:
    spectral_guard.SpectralGuard:
      python_path: .

Criteria (deliberately loose physical windows -- they catch the 10^3-10^34 pathology,
never a legitimate cosmology):
  - D_l^TT(220) in [1000, 20000] muK^2   (first acoustic peak)
  - no negative C_l^TT for 30 <= l <= 2000
  - D_l^TT(2) < 5000 muK^2               (ISW plateau; pathology gives >~ 10^5)
"""
import numpy as np
from cobaya.likelihood import Likelihood

T0_MUK = 2.7255e6


class SpectralGuard(Likelihood):
    def get_requirements(self):
        return {"Cl": {"tt": 2000}}

    def logp(self, **params_values):
        cl = self.provider.get_Cl(ell_factor=False, units="1")["tt"]
        ll = np.arange(len(cl))
        Dl = cl * ll * (ll + 1) / (2 * np.pi) * T0_MUK**2
        lmax = min(2000, len(Dl) - 1)
        if not (1000.0 < Dl[220] < 20000.0):
            return -np.inf
        if (Dl[30:lmax] < 0).any():
            return -np.inf
        if Dl[2] > 5000.0:
            return -np.inf
        return 0.0
