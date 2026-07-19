# Tier-2 #4 — full-likelihood SEDE vs ΛCDM head-to-head (2026-07-17)

The goal: the real `Δχ²(SEDE−ΛCDM)` against **DESI DR2 BAO + Pantheon+ SN + Planck 2018**, through the
actual **mochi-class** Boltzmann code inside **cobaya**, replacing the compressed/Fisher proxies
(Stage-1.5 `Δχ²=−0.04`, the retracted `ΔDIC≈−2.9`).

**This note documents a real diagnostic journey, including a wrong turn I corrected.** Short version:
1. cobaya was loading the wrong (vanilla, no-smg) CLASS → fixed with `path: global`.
2. The SEDE CMB then came out grossly broken → I *briefly* concluded "CMB uncomputable, ISW retracted."
3. That was wrong: the breakage was a **numerical pathology of the `z_gr_smg=99` setting**, not physics. I
   diagnosed it to a spurious scalar-field oscillation and fixed it with `z_gr_smg=5`. **The CMB is
   computable; the ISW prediction (`C_ℓ^TT(ℓ=2)=0.924`) stands and matches the corpus Stage-2 value.**
4. Full BAO+SN+Planck minimization: **`Δχ²(SEDE−ΛCDM) = −3.4`**, SEDE mildly preferred (the preference in
   Planck highl-lite). The `μ_∞=1.05` growth enhancement raises CMB lensing in the direction of Planck's
   mild `A_lens>1` preference — a +200 "tension" only if SEDE is unfairly evaluated at ΛCDM's optimum; the
   joint fit absorbs it. Bestfit not marginalized; see caveats in §4.

---

## 1. Two setup bugs fixed

- **Wrong CLASS.** cobaya silently loaded the vanilla CLASS that `cobaya-install` put under
  `--packages-path` (no smg), so every SEDE eval died (`Omega_Lambda/Omega_fld must be unspecified`). Fix:
  `path: global` in the classy theory block → cobaya imports the mochi build from `.venv`.
- **`z_gr_smg` pathology (the important one).** With the corpus's `z_gr_smg = 99`, the SEDE lensed `D_ℓ^TT`
  is ~1500× too large at the first peak and goes *negative* at ℓ≳500. **Diagnosis** (perturbation dump at
  `k_output_values`): at `k ≈ 0.02/Mpc` the metric combination `(φ+ψ)` develops a late-time (z<4)
  oscillation reaching amplitude ~19 (vs 0.5 physical) with a huge `d(φ+ψ)/dτ` → a *fake* ISW source that
  swamps `C_ℓ`. It is seeded at the GR→smg transition and is scale-selective (k~0.02 only; k≥0.1 is fine,
  which is why `P(k)`/σ8 were always sane). It is **not** an integration-tolerance artifact (tol 1e-5 vs
  1e-8 identical) and no `method_qs_smg` / IC / QS-trigger touches it. **Fix:** lower `z_gr_smg` to **5** —
  the scalar reaches its attractor before smg is dynamically relevant (`Ω_DE(z=5)≈1%`, so no physics is
  lost). Then the oscillation vanishes (`|φ+ψ|=0.52`, 0 turns), and the spectrum is **converged and
  z_gr-independent**: `D_ℓ(220)` and `C_ℓ(ℓ=2)` agree to <0.02% across `z_gr ∈ {5,10,20,30,50}`, and the
  first peak matches ΛCDM to the digit (`D_ℓ(220)=5808`, ratio 0.9999 — correct, since recombination is GR).
  `z_gr=99` was verified sane-boundary-crossing; `z_gr=5` is sane across the whole `Ω_m` box (tested to
  `Ω_m≈0.38`). Both `mcmc/sede.yaml` and `src/sims/sede_fsigma8_isw_finitek.py` now use `z_gr_smg=5`.

## 2. Background result — DESI DR2 BAO + Pantheon+ SN (+ BBN prior)

Direct scipy minimization (Nelder-Mead → Powell, 4 starts; minima robust to `±0.000`). BBN prior
`ω_b=0.02218±0.00055` applied identically to both. Script: `mcmc/bg_minimize.py`.

| model | χ² total | χ²_BAO | χ²_SN | H0 | ω_cdm | Ω_m |
|---|---|---|---|---|---|---|
| ΛCDM | 1416.18 | 10.89 | 1405.29 | 68.76 | 0.1215 | 0.305 |
| SEDE | 1415.75 | 10.75 | 1405.00 | 69.13 | 0.1248 | 0.306 |

**`Δχ²(SEDE − ΛCDM) = −0.43`** — indistinguishable at the geometry level, SEDE with zero extra background
parameters (w(z) fixed by the FPAB file). SEDE prefers H0=69.1 vs 68.8. Consistent with Stage-1.5 (−0.04).

## 3. CMB is computable — ISW prediction confirmed

With `z_gr_smg=5`, the SEDE CMB signature (converged, z_gr-independent):
- **Primary peaks: unmodified** (`D_ℓ(220)` ratio 0.9999) — as required, recombination is in the GR regime.
- **Low-ℓ ISW suppression:** `C_ℓ^TT` ratio SEDE/ΛCDM `= 0.9236 (ℓ=2)`, `0.987 (ℓ=20)` → back to 1 by
  ℓ~100. This is a real **−7.6% ISW at ℓ=2**, matching the corpus's Stage-2 value (0.923) and TIER1 #1
  (0.924) exactly. **The mid-Tier-2 retraction of these numbers was an artifact of `z_gr=99` and is itself
  retracted** (banners corrected in `TIER1_results.md`, `D6_result_injection_growth.md`, pre-reg P5).
- **Matter growth (independent, always was sane):** `P_SEDE/P_ΛCDM = +2.6%` flat over `k=0.01–1 h/Mpc`,
  σ8 ratio `+1.3%`, `+2–3%` fσ8 (reproduces Stage-2 `1.0261`). The `+2–3%` fσ8 prediction stands.

## 4. Full Planck+BAO+SN head-to-head — the growth/lensing tension

**The apparent tension is an artifact of the comparison point, not physics.** At *fixed* ΛCDM-bestfit
parameters SEDE's highl-lite is ~+200 worse (916 vs 716) — the fingerprint of `μ_∞=1.05` raising CMB
lensing (smoothed damping tail). A_s alone cannot absorb it (lowering A_s mismatches the primary peaks).
**But the *joint* fit does:** letting `{n_s, ω_cdm, H0, A_s}` re-optimize together absorbs it entirely with
tiny shifts. So the "+200" is not a rigid tension — it was an unfair comparison of SEDE at ΛCDM's optimum.

Full-likelihood minimization (both models, same 5 free params `{H0, ω_b, ω_cdm, logA, n_s}`, τ and A_planck
fixed identically; direct scipy multi-start — cobaya's Nelder-Mead stalls at the ref simplex). Both starts
per model agree to <0.1 in χ² (robust minima). Script: `mcmc/full_minimize.py`, log `mcmc/out/full_minimize.log`.

| model | χ² total | BAO | SN | lowl.TT | highl TTTEEE-lite | H0 | ω_cdm | n_s |
|---|---|---|---|---|---|---|---|---|
| ΛCDM | 2023.61 | 11.53 | 1406.00 | 22.86 | 583.22 | 68.53 | 0.1185 | 0.9682 |
| SEDE | 2020.18 | 11.19 | 1406.70 | 22.38 | **579.91** | 68.86 | 0.1193 | 0.9664 |

**`Δχ²(SEDE − ΛCDM) = −3.4`** — SEDE mildly *preferred*. Per-probe: BAO −0.34, SN +0.70, lowl.TT −0.48,
**highl TTTEEE-lite −3.31** (the dominant term). So SEDE fits Planck's high-ℓ TTTEEE *slightly better*, not
worse. Physically plausible: `μ_∞=1.05` supplies extra CMB lensing in the **same direction as Planck's
known mild `A_lens>1` preference** — the flip side of the "+200 at fixed params" coin.

**Significance and caveats (do not overclaim):**
- This is a **bestfit Δχ², not a marginalized significance.** With SEDE at effectively the same parameter
  count as ΛCDM (w(z) fixed), `Δχ²=−3.4` is a *mild* preference (~1.8 naive-`√Δχ²`); the real significance
  and any degeneracy-sliding need the posterior (the production MCMC we deferred).
- This run has **no BBN prior** and fixed τ, A_planck; both models rail `ω_b→0.0225` (a bit high vs the
  canonical 0.02237). Identical for both, so `Δχ²` is a fair comparison, but the absolute fit sits in a
  slightly non-standard `ω_b` corner. A BBN-pinned rerun would tighten this.
- **What is robust:** the *sign and rough size* — SEDE is not disfavored by Planck; it is marginally
  favored, and the preference lives in Planck highl-lite (lensing), now shown with the **real** likelihoods
  rather than the compressed-CMB/SH0ES proxy behind the corpus's `ΔDIC≈−2.9`. That earlier number is thus
  *qualitatively corroborated* (a real, small preference) by a direct computation — while the specific
  `−2.9` value and its CMB-distance provenance remain superseded by this `−3.4` (BAO+SN+real-Planck).

## 4b. PRODUCTION MCMC — the marginalized verdict (completed 2026-07-17)

Both chains converged (R−1 means: ΛCDM 0.014 / SEDE 0.015; bounds < 0.11; 13k / 19.6k accepted samples;
`covmat: auto` from `base_plikHM_TTTEEE_lowl_post_BAO`; spectral guard active throughout — zero pathological
samples entered either chain). Configs `mcmc/{lcdm,sede}_mc.yaml`; triangle plot
`src/figures/sede_lcdm_mcmc_triangle.png`.

| | ΛCDM | SEDE |
|---|---|---|
| H0 | 68.54 ± 0.29 | 68.89 ± 0.30 |
| ω_cdm | 0.11847 ± 0.00063 | 0.11922 ± 0.00066 |
| n_s | 0.9685 ± 0.0032 | 0.9665 ± 0.0033 |
| ⟨χ²⟩ | 2028.1 | 2025.1 |
| p_D | 6.8 | 6.3 |
| **DIC** | 2034.9 | **2031.4** |

**`ΔDIC(SEDE − ΛCDM) = −3.5`** (Δ⟨χ²⟩ = −3.1; in-chain Δχ²_min = −2.6, consistent with the minimizer's
−3.4). With identical priors and the same parameter count, this is a **marginalized, real-likelihood
confirmation that SEDE is mildly preferred** — "positive but not strong" evidence (ΔDIC 2–5 band),
roughly the ~1.9σ equivalent. The preference survives marginalization: it is not a bestfit artifact.

**The narrative point:** the corpus's original `ΔDIC ≈ −2.9` — which the critique (R2) correctly
flagged as compressed-CMB/SH0ES-driven and Stage-1.5 deflated — is now **reproduced at −3.5 from the
full DESI-DR2 + Pantheon+ + real-Planck likelihoods through the actual Boltzmann code, marginalized over
the full posterior.** The number was approximately right; its provenance was wrong. This replaces it with
a defensible one.

Secondary posterior shifts (the physics of how SEDE absorbs the extra lensing): SEDE pulls **H0 up
+0.35 (0.8σ)**, `n_s` **down** 0.002, `ω_cdm` **up** 0.0008 — small, coherent, and the reason the
fixed-parameter "+200 tension" of §4 evaporates in the joint fit.

## 5. Bottom line

- **Geometry (DESI+Pantheon+):** SEDE ≈ ΛCDM, `Δχ²=−0.43`, zero extra background parameters.
- **Full BAO+SN+Planck:** bestfit `Δχ² = −3.4` AND marginalized **`ΔDIC = −3.5`** (production MCMC, §4b)
  — SEDE **mildly preferred**, the preference living in Planck highl-lite, surviving full
  marginalization at equal parameter count. The corpus's `ΔDIC≈−2.9` is thereby replaced by a
  real-likelihood, Boltzmann-computed, posterior-marginalized `−3.5`.
- **CMB is computable** (with `z_gr_smg=5`; the corpus's `z_gr=99` is numerically pathological). The
  **−7.6% low-ℓ ISW** prediction is confirmed and reproducible (`C_ℓ(ℓ=2)=0.9236`).
- **Growth:** `+2–3%` fσ8 / `+2.6%` P(k) — reliable, falsifiable (DESI DR3 / Euclid), wrong way for S8.
- **CMB lensing, not a tension:** `μ_∞=1.05` raises CMB lensing in the direction of Planck's mild
  `A_lens>1` preference; it looks like a +200 tension only if SEDE is (unfairly) evaluated at ΛCDM's
  optimum — the joint fit absorbs it and turns it into the −3.3 highl-lite gain.
- **Solver hardened (production-MCMC gate cleared).** Full diagnosis: the FPAB `stable_params`
  reconstruction has an **early-time/low-k tachyonic window** in its scalar sector (QS `mass² < 0` at
  `k ≲ 0.05/Mpc`, `z ≳ 30–50`, cosmology-dependent). Releasing the scalar inside that window (any `z_gr_smg`
  above it, e.g. the corpus's 99) rings the spurious oscillation *regardless of ICs* — verified: (a) the QS
  IC at the z=99 switch is 3500× off-attractor (`V_x = 6×10⁶` vs attractor `~1.7×10³`), and (b) starting
  from `x = x′ = 0` still rings (the window is genuinely unstable, not just badly initialized). A per-mode
  deferral of the GR→smg switch was tried and conflicts with the global QS scheduling (deeper surgery,
  rejected). **The defensible package, implemented:**
  1. **IC guard in the C code** (`source/perturbations.c`, GR-off switch): if the QS `mass² ≤ 0` at the
     switch, the scalar starts at `x = x′ = 0` (continuous with the GR-held regime) instead of the
     ill-conditioned QS value. Defense-in-depth; rebuilt into the venv classy.
  2. **`z_gr_smg = 5` validated exhaustively:** z_gr-independence to 3×10⁻⁵ across `z_gr ∈ {3,5,10,20}`;
     **49/49 grid points sane** over `H0 ∈ [62,75] × ω_cdm ∈ [0.104,0.140]` (Ω_m ≈ 0.27–0.42); 4 hostile
     corners with `A_s/n_s` excursions all pass. `Ω_DE(z=5) ≈ 1%`, so no dark-energy physics is lost.
  3. **Per-sample spectral guard** (`mcmc/spectral_guard.py`, a cobaya likelihood in both configs): returns
     `−inf` for a pathological spectrum (verified: passes the sane `z_gr=5` spectrum at exactly 0, rejects
     the `z_gr=99` pathology with `−inf`). No corrupted sample can enter a chain silently.
  4. **BBN prior** `ω_b = 0.02218 ± 0.00055` added identically to both yamls (fixes the railed-`ω_b` caveat
     of §4 for future runs).
  A production MCMC can now be launched on these configs when a marginalized significance is wanted.

## Files
- `mcmc/full_minimize.py` — full-likelihood BAO+SN+Planck minimizer (both models); `mcmc/out/full_minimize.log`.
- `mcmc/bg_minimize.py` — background BAO+SN(+BBN) minimizer (§2).
- `mcmc/{lcdm,sede}.yaml` — full configs (`path: global`; SEDE `z_gr_smg=5`); `mcmc/{lcdm,sede}_bg.yaml`.
- Diagnostics: `mcmc/debug_sede_{classy,cobaya,exact}.py`.
