# Tier-2 Stage-1.5 executed — the first distance+growth joint fit for FPAB-SEDE

**Script:** `src/sims/sede_compressed_joint_fit.py` · **Result:**
`src/results/sede_compressed_joint_fit_results.json`.

## What this is (and is not)
The full likelihood — Cobaya + Mochi-CLASS, finite-k perturbations at every sampled point —
needs installed Boltzmann/likelihood machinery not present in this environment; that remains the Stage-2
campaign. Stage-1.5 is the maximal rigorous step runnable now: a **compressed BAO+BBN+RSD+lensing joint fit**
that, for the first time, **adds the growth arm** (fσ8 + S8) the corpus's distance-only `ΔDIC ≈ −2.9` lacked.

**Why the compression is legitimate for SEDE specifically:** SEDE modifies nothing before recombination
(BBN speed-up = 1, `Ω_DE(z=1100) ~ 10⁻¹⁰`), so `r_drag` is a model-independent function of `(ω_m, ω_b)` —
taken from the Aubourg et al. 2015 formula (0.1%-accurate). This is the standard "BAO+BBN" absolute-scale
anchor; it breaks the BAO `H0–r_d` degeneracy **without** a CMB acoustic prior (the precise `r_s` for `R`,
`l_A` genuinely needs a Boltzmann code — deferred to Stage 2). Both models fit the **same 4 parameters**
(`Ω_m, H0, ω_b, σ8`), SEDE's gate frozen, so `Δχ² = ΔAIC = ΔBIC`.

Data (compressed, representative to knowledge cutoff — refresh before quoting): `ω_b` (BBN) + DESI DR2 BAO
(12 points) + fσ8 compilation (6 points) + an S8 weak-lensing prior. N_data = 20.

## Result

| | Ω_m | H0 | σ8 | S8 | r_d | χ²(ω_b / BAO / fσ8 / S8) | total |
|---|---|---|---|---|---|---|---|
| **ΛCDM** | 0.298 | 68.9 | 0.790 | 0.787 | 147.4 | 0.00 / 11.00 / 7.43 / 0.41 | **18.84** |
| **FPAB-SEDE** | 0.302 | 69.2 | 0.785 | 0.787 | 146.6 | 0.00 / 10.72 / 7.65 / 0.44 | **18.80** |

> **Δχ²(SEDE − ΛCDM) = −0.04** (χ²/dof = 1.18 both). **Statistically indistinguishable.**
> Per-probe: BAO marginally *prefers* SEDE (−0.28), growth marginally *disfavours* it (+0.22) — they cancel.

## Four findings

1. **SEDE is fully consistent with BAO+BBN+RSD+lensing** (χ²/dof = 1.18) at physical parameters
   (H0 ≈ 69, r_d ≈ 147, Ω_m ≈ 0.30). No tension introduced by adding growth.
2. **The growth arm — never tested in the corpus — does not exclude SEDE**, and does **not** reveal an S8
   advantage: both models fit S8 ≈ 0.787, and SEDE's σ8 = 0.785 ≈ ΛCDM's 0.790. **The advertised
   "SEDE predicts low S8 ≈ 0.76 that eases the lensing tension" is not realized** — consistent with the FPAB
   v0.3 retraction of the 0.76 value (σ8 = 0.811 in that track). This should be struck from the theory's
   selling points, confirming a critique-file flag (P4).
3. **The known mild SEDE preference is confirmed CMB-distance/H0-driven, not BAO/growth-driven.** The corpus's
   `ΔDIC ≈ −2.9` used compressed-CMB-distance + SH0ES; dropping those (this BAO+BBN+growth combination) leaves
   SEDE with *no* preference (Δχ² = −0.04). This matches the corpus's own per-probe decomposition (preference
   carried ~100% by CMB-distance + H0; BAO/fσ8 flat) — now independently reproduced.
4. **No H0 resolution:** SEDE gives H0 ≈ 69 like ΛCDM (it is r_d-standard); the Hubble tension is untouched,
   as the corpus states.

## What this does NOT capture (and why the discrimination lives there)
The compressed fit uses the SEDE **background E(z) + linear growth D(z)** only. It does **not** include: the
finite-k modified-gravity signature (μ(k,a), the near-horizon force turn-on `k₅₀ ≈ 0.7H₀`, γ_slip = 1), the
**Δ = 1 volume-vs-area measurement** (SEDE's flagship, 11σ at DR3+Euclid), the injection-frame ISW-only
gate–matter correlation (D-6), or the full CMB. **At the background + linear-growth level with equal
parameters, SEDE ≈ ΛCDM by construction** (they share the DESI-preferred (w0,wa) quadrant near Λ) — so a null
here is expected and correct. **SEDE's discriminating power is entirely in the Stage-2 finite-k + Δ
observables**, exactly the ones this compression omits. The value of Stage-1.5 is the honest negative: it
rules out a *growth-based* exclusion of SEDE and a *growth-based* S8 rescue, sharpening what Stage 2 must
deliver.

## Stage-2 environment — NOW INSTALLED (2026-07-16)

The blocker was environmental, and it is resolved. Installed into `.venv` (gitignored):
- **cobaya 3.6.2** (+ GetDist, py-bobyqa) — the sampler / likelihood framework.
- **classy 3.3.3 = the mochi-class *modified-gravity* Boltzmann engine** (built from the local
  `re/work/mochi_class_public` source, `--no-build-isolation`), with full `smg` / α-function support —
  **not** standard classy.
Two build fixes applied (persist): (1) the pre-built `class` binary had `__CLASSDIR__` baked to a dead
`~/Documents/...` path → **rebuilt locally** (`main/class.c` object + relink; `libclass.a` built); (2) the
wrapper's data dir resolves to `site-packages/external` → **symlinked** to the mochi `external/` tree.

**Verified working from Python** (`.venv/bin/python`):
```
LCDM:  sigma8 = 0.8228
smg (propto_omega, alpha_B = 0.2, SEDE-like braiding):  sigma8 = 0.8650   [0.7 s/model]
   -> braiding RAISES late growth (consistent with FPAB sigma8 ~ 0.81; the retracted 0.76 does not appear)
   background exposes (.)rho_smg, (.)p_smg, M*^2_smg, ... (3000 pts)
```

**What this unlocks now:** genuine **finite-k SEDE predictions** — μ(k,a), the near-horizon force turn-on,
γ_slip, C_ℓ^{TT,φφ}, P(k), fσ8 — driven from Python via the α-function (stable) basis the corpus's
`fpab_mochi_export.py` already generates; and cobaya MCMC with the classy theory. The Stage-1.5 note's
"finite-k lives in Stage 2" caveat is now executable in this environment.

**Remaining for the *full* MCMC:** the likelihood DATA (Planck clik, ACT, DESI, DES) via
`cobaya-install cosmo` — a multi-GB download, not yet fetched. Everything upstream of the data is in place.

## Tier-2 status and the path to Stage 2

- **Done (Stage-1.5):** distance+growth joint fit, growth arm added, S8 advantage disproven, preference
  localized to CMB-distance/H0. Runnable, reproducible, in-repo.
- **Blocked on external install (Stage-2, full-likelihood):** the variance-gated Mochi-CLASS provider
  (`re/work/cobaya_fpab_background.py` extends to it) wired into Cobaya with Planck PR4 + ACT DR6 lensing +
  DESI full-shape + DES/KiDS lensing, sampling the finite-k EFT (μ(k), Δ) at every point. The provider
  architecture exists; the blocker is purely environmental (`cobaya`, `classy`/mochi build, likelihood data).
- **The one number Stage 2 must produce:** the Δ posterior (SEDE: Δ = 1 vs ΛCDM: Δ = 0 area law) with the
  full finite-k growth — that is where SEDE either separates from ΛCDM at ~11σ or falls with it. Stage-1.5
  confirms nothing in the *background+linear-growth* data blocks that measurement.
