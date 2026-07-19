# The Planck-lensing (φφ) test — result

**Question (obstest §5/§7/§8):** is SEDE's extra high-ℓ TTTEEE smoothing physical lensing,
or is it absorbing the Planck A_lens>1 lite-likelihood systematic? The direct Planck 2018
lensing (φφ) likelihood — whose reconstructed amplitude sits at ≈1.0, unlike the TT-derived
A_lens≈1.2 — decides.

## Definitive test: BAO/SN/BBN-anchored joint refit (`phiphi_joint_refit.py`)

The clean estimator: refit **both** models over {H0, ω_b, ω_cdm, A_s, n_s} to the full
**DESI DR2 BAO + Pantheon+ SN + BBN prior + lowl TT + TTTEEE-lite** data, once **without** and
once **with** the Planck 2018 lensing likelihood (`planck_2018_lensing.native`, CMB-marginalized).
BAO+SN+BBN pin the background density (ω_cdm stays ≈0.119, H0≈68.5/68.8 — no drift), and both
models are properly refit, so the difference of joint minima is the honest post-φφ preference.
τ=0.0544, A_planck=1.0025 fixed. Powell; from `src/mcmc/`, `COBAYA_PACKAGES` = a full cobaya
packages dir (Planck lowl+plik_lite+lensing, DESI DR2 BAO, Pantheon+).

| | χ²_joint (no φφ) | χ²_joint (+φφ) | χ²_φφ at joint bf |
|---|---|---|---|
| ΛCDM | 2023.91 | 2032.75 | 8.82 |
| SEDE | 2020.40 | 2030.44 | 9.99 |

- **Δχ²_joint (no φφ) = −3.5** — reproduces the best-fit −3.4 and marginalized DIC −3.5.
- **Δχ²_joint (with φφ) = −2.3** — **the honest post-φφ preference.**
- φφ penalty (SEDE−ΛCDM) at the joint best fit = **+1.2**; the direct lensing datum erodes the
  preference from −3.5 to −2.3 (≈⅓; 66% survives). SEDE fits φφ worse because μ_∞=1.05 boosts
  C_ℓ^φφ above the measured amplitude.

χ²_φφ(ΛCDM) ≈ 8.8 at the joint best fit is consistent with (though not a proof of) a well-calibrated
∼9-bandpower likelihood. The φφ likelihood is CMB-marginalized, so it shares some Planck TT
information with the primary fit and is not a fully independent datum.

## Earlier estimates (superseded by the joint refit above)

Two cruder estimates were computed first; the joint refit supersedes both.

| estimate | Δχ²_primary | Δχ²_post-φφ | φφ penalty | issue |
|---|---|---|---|---|
| A_s-profiled at primary optimum (`phiphi_As.py`) | −3.17 | −1.90 | +1.27 | χ² SUM at the primary optimum, not a joint fit (upper bound) |
| CMB-primary-only full refit (`phiphi_fullrefit.py`) | −2.55 | −2.20 | +0.30 | ω_cdm drifts to 0.120, unpinned by BAO/SN |

The A_s-profiled penalty (+1.27) is close to the anchored-joint +1.17, confirming it was the more
faithful of the two; the CMB-only +0.30 is an artifact of the unphysical ω_cdm drift (which BAO/SN
forbid). The properly anchored **joint refit (−2.3)** is the number the paper quotes.

## Interpretation

The decisive test lands on the *partly-systematic* side: roughly a third of the high-ℓ TTTEEE gain is
**consistent with** the absorbed A_lens systematic rather than physical lensing (the test cannot
distinguish an absorbed systematic from SEDE's μ_∞ simply over-predicting real lensing). The residual
(≈−2.3, 66%) is a genuine but weaker pull. This consolidates the "statistically consistent with ΛCDM,
mild pull" reading; the direct lensing datum erodes the preference and never strengthens it.

**Caveats.** φφ is the CMB-marginalized native likelihood (shares Planck TT information); τ, A_planck
fixed; single Boltzmann fork; a full marginalized DIC with φφ in the chains (vs this best-fit refit)
is the remaining refinement.
