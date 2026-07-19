# Reproducing the analysis

The analysis needs three components, none machine-specific:

1. **The mochi-class fork** (the modified-gravity Boltzmann code). Build the pinned fork
   (`74d4b05` solver guard + `b365166` frozen FPAB inputs):
   ```
   pip install cython numpy            # build deps
   cd <mochi-class fork>
   rm -f python/classy.cpp; rm -rf python/build   # clear stale artifacts (known gotcha)
   pip install . --no-build-isolation --force-reinstall --no-deps
   python -c "import classy"           # must expose gravity_model / stable_params / Omega_smg
   ```
   Then point the scripts at it: `export MOCHI_DIR=<mochi-class fork>`.

2. **cobaya + likelihood data**: `cobaya-install mcmc/likes_install.yaml -p mcmc/packages` (run from
   `src/`). This installs everything the scripts use: plik_lite TTTEEE + lowl TT, **Planck 2018 lensing
   (φφ)**, DESI DR2 BAO and Pantheon+ SN (Pantheon+ pulls a network download). The scripts default to
   `packages_path = src/mcmc/packages` (both `full_minimize.py` and the `phiphi_*` scripts), so installing
   into that path makes them run as shipped; otherwise `export COBAYA_PACKAGES=<your packages dir>` and the
   `phiphi_*` scripts pick it up. The `sims/` growth and ISW predictions need none of this, only `MOCHI_DIR`.

3. **The frozen SEDE inputs** are vendored here in `mcmc/inputs/`
   (`fpab_sede_b020_{stable,rho}.dat`); the cobaya yamls reference them relatively, so run
   cobaya from `src/mcmc/`.

## What runs

- `sims/sede_fsigma8_isw_finitek.py` — the growth (fσ8, σ8) and ISW predictions through
  mochi-class (reproduces C_ℓ^TT(2) ratio ≈ 0.924, −7.6% ISW, +1.3% σ8, +2% scale-independent
  fσ8). Needs only `MOCHI_DIR`.
- `mcmc/{lcdm,sede}_mc.yaml` — the production MCMC configs (cobaya). Run from `src/mcmc/`.
- `mcmc/full_minimize.py` — the best-fit minimizer.
- `mcmc/phiphi_joint_refit.py` — **the definitive Planck-lensing (φφ) test** (§5): a BAO/SN/BBN-anchored
  joint refit of both models without and with φφ, giving Δχ² = −3.5 → −2.3. Run from `src/mcmc/`; needs
  the full data set from step 2. `phiphi_As.py` (A_s-profiled) and `phiphi_fullrefit.py` (CMB-only refit)
  are the earlier, superseded estimates; see `results/phiphi_test_result.md`.

## Not yet committed

The production chain files, the minima, and the chain→(⟨χ²⟩, p_D, DIC) reduction script are to
be committed with the tagged release so the headline ΔDIC is auditable (see the paper's
Reproducibility section). `CHECKSUMS.txt` lists their SHA-256.
