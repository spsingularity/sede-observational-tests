# Tier-1 calculations executed (2026-07-17)

Three desk-doable calculations run with the installed engine. Honest outcomes: **#1 produced a real new
prediction and corrected an earlier overclaim of mine; #2 came out deflationary (the 3-loop does not tighten,
for a physical reason); #3 turned out not to be a fixed calculation at all.** No overclaiming — the results
are what the numbers say.

---

## #1 — finite-k `fσ8` + ISW through mochi-class  → NEW PREDICTION, and it CORRECTS D-6

> ✅ **CONFIRMED — a brief retraction was itself retracted (Tier-2, `TIER2_head_to_head_result.md`,
> 2026-07-17).** Both halves of this result stand. The **fσ8** part is from the sub-horizon `P(k)` (`+2.6%`
> flat, `σ8 +1.3%`, matches Stage-2 `1.0261`). The **ISW** part (`C_ℓ^TT ℓ=2 = 0.924`, `−7.6%`) is **real
> and reproducible** — but *only* with `z_gr_smg ≤ 5–50`. Mid-Tier-2 I briefly retracted it after runs at
> the corpus's `z_gr_smg = 99` gave a broken CMB (`D_ℓ^TT` ~1500× too large, negative at high ℓ). That was
> diagnosed to a **spurious scalar-field oscillation** excited at the GR→smg transition at `k ≈ 0.02/Mpc`
> (a *fake* ISW), NOT physics: lowering `z_gr_smg` to 5 (scalar on its attractor before smg is dynamical)
> gives a **converged, z_gr-independent** spectrum with `C_ℓ^TT(ℓ=2) = 0.9236` — exactly the Stage-2 value.
> The ISW prediction is intact; the mochi configs (`mcmc/sede.yaml`, `sede_fsigma8_isw_finitek.py`) now use
> `z_gr_smg = 5`.

**Script:** `src/sims/sede_fsigma8_isw_finitek.py`.

| k [h/Mpc] | z=0 | z=0.3 | z=0.5 | z=1.0 |
|---|---|---|---|---|
| 0.005 (~horizon) | +2.76% | +2.15% | +1.96% | +1.72% |
| 0.05 | +2.93% | +2.27% | +2.05% | +1.73% |
| **0.10 (RSD)** | **+2.92%** | **+2.26%** | **+2.04%** | **+1.73%** |
| 0.20 | +2.92% | +2.27% | +2.04% | +1.73% |

`ISW C_ℓ^TT` ratio: `ℓ=2: 0.924`, `ℓ=5: 0.940`, `ℓ=10: 0.966` → back to 1 by `ℓ~200`.

**The correction to D-6 (a self-catch).** D-6 claimed the modification is "horizon-confined, RSD-safe,
`δfσ8 ~ 7×10⁻⁵`." That is **wrong for `fσ8`.** D-6 conflated two different quantities:
- the **DE clustering** `δ_DE/δ_m` — which *is* horizon-suppressed (~`7×10⁻⁵` at RSD; dark energy doesn't
  cluster), and
- the **matter growth rate** — which feels the effective gravitational coupling `μ_∞ = 1.05` at *all*
  sub-horizon scales, and is **not** horizon-confined.

The observable `fσ8` is the second one. So SEDE predicts a **`+2–3%` `fσ8` enhancement, essentially
scale-independent across `k = 0.005–0.2 h/Mpc`** (i.e. right through the RSD range), plus a **`−7.6%` ISW
suppression at `ℓ=2`.**

**Consequences:**
1. **SEDE is NOT RSD-invisible.** The `+2–3%` `fσ8` is within reach of DESI DR3 + Euclid (current `fσ8`
   errors ~5–8%). This is a genuine, testable prediction — better than the "invisible" D-6 story.
2. **It goes the wrong way for the S8/growth tension** (higher `fσ8`, not lower) — quantifying the P4/P5
   finding that SEDE *raises* growth. The retracted `σ8=0.76`-helps-S8 story is now contradicted at the
   `fσ8` level too.
3. **Pre-registration P5 is corrected:** the claim "horizon-only, strictly nothing at RSD scales" is **false**
   — there is a scale-independent `+2–3%` RSD signal from `μ_∞`. What *is* horizon-confined is the DE
   clustering (`δ_DE`), not the growth signature. P5 should read: "`+2–3%` scale-independent `fσ8` enhancement
   (from `μ_∞`) + a `−7.6%` low-ℓ ISW; the *DE-perturbation* `δ_DE` is horizon-confined."

This is the most valuable Tier-1 outcome: a real prediction, and the removal of an overclaim I made in D-6.

## #2 — 3-loop mirror-QCD running  → deflationary (does not tighten)

**Script:** `src/sims/mirror_qcd_3loop.py`. With `α_h(M_Pl) = α_s(M_Pl) = 1/53.0` and `N_f=6`:

| loop | `Λ_pole` |
|---|---|
| 1 | 28 MeV |
| 2 | 105 MeV |
| 3 | (spurious large-α fixed point — β develops a zero at `a≈0.28`, outside perturbative control) |
| 4 | 133 MeV |

**The 3-loop does not tighten the band; it reveals why.** `N_f=6 SU(3)` sits close enough to strong coupling
(and to the conformal-window edge, `N_f* ~ 8–12`) that the confinement scale is **loop-order-sensitive** —
1-loop → 4-loop span a factor ~5 (28–133 MeV), and the 3-loop β-function has a spurious IR fixed point at
`α≈3.5` where perturbation theory has already failed. So the mirror-`N_f=6` scale is bracketed to the **right
decade** (~30–130 MeV, consistent with the transmutation `μ`), but perturbation theory **cannot pin it** to
the bare `μ=28.5` vs the UZ-corrected `52 MeV`. The physical value needs the lattice `T_c/Λ_MS̄` ratio (an
`O(1)` nonperturbative input). Honest verdict: the zero-parameter mirror-unification lands `μ` in the right
decade — a real order-of-magnitude success — but the precise value is a lattice question, not a perturbative
one. (This is a *negative* refinement of the earlier "band 28–105" — 3-loop confirms the band is irreducibly
`O(1)`-wide from perturbation theory near the conformal edge.)

## #3 — S₃ `Pin⁺ ν mod 16` anomaly  → not a fixed calculation

> ✅ **UPDATE (2026-07-17, `PIN16_assignment_constructed.md`).** The assignment has now been **exhibited**:
> all ~16 CP-relevant Weyl fermions organize into CP-massable pairs under the report's own spurion
> structure, each pair contributes `+1−1 = 0`, so **ν = 0 mod 16 automatically** — the per-pair sign
> freedom never changes the sum. Status: `[OPEN — bookkeeping]` → `[CONSTRUCTED]` (conditional only on
> CP-real dressings at the symmetric point, which exist by per-pair phase absorption). The verdict below
> ("not a *calculation*") stands — it is a construction — but the existence question is now answered.

The gauged-CP `Pin⁺` `ν mod 16` index is, per the S3_GAUGING report's own words, **"bookkeeping"** — it depends
on the **intrinsic CP phases assigned to each Majorana/vectorlike fermion** (Φ_V/Φ_D-inos, the ℤ₃₁ spectators
C/D/Y₁₋₄, and the R3 8-Weyl spectator set). Those phases are a **model-building choice** in the gauged-CP
embedding, not fixed by the theory. So `ν` is **assignment-dependent**: the honest statement is that it is a
discrete constraint (`ν ≡ 0 mod 16`) to be *satisfied by choosing* the CP phases over a finite set — an
existence question, not a computed number. I can frame it, but **calling it "calculated" would be exactly the
kind of overclaim this whole audit targets.** Verdict: it stays `[OPEN — bookkeeping, satisfiable-by-choice]`,
not closed. (Resolving it firmly needs the full CP-phase assignment table, which is a construction, not a
calculation.)

---

## Tier-1 scoreboard

| # | result | status |
|---|---|---|
| **#1 fσ8/ISW** | `+2–3%` scale-independent `fσ8` (RSD-visible) + `−7.6%` ISW; **corrects D-6** and P5 | ✅ real new prediction |
| **#2 3-loop μ** | band `28–133 MeV`, not tightened; near conformal edge → lattice needed | ✅ done, deflationary |
| **#3 Pin⁺ anomaly** | assignment-dependent bookkeeping, not a fixed number | ⚠ not calculable as posed |

**Net:** one genuine new prediction (that also removes one of my own overclaims), one honest negative
refinement, and one item correctly identified as un-calculable-as-posed. The `fσ8` result (#1) is the keeper —
it makes SEDE *more* testable (a `+2–3%` RSD signal DESI/Euclid can see) and it is propagated as a correction
to `D6_result_injection_growth.md` and pre-registration P5.
