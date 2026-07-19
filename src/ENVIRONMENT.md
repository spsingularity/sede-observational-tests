# Pinned environment

- python 3.12; cobaya 3.6.2; getdist (cobaya dep); numpy/scipy per requirements.txt
- Boltzmann engine: mochi-class fork, commits `74d4b05` (QS-IC solver guard) + `b365166`
  (frozen FPAB-SEDE inputs: stable_params .dat + inifiles). Build:
  `rm python/classy.cpp && rm -rf python/build && pip install . --no-build-isolation`
  (the stale-cpp removal is required; symlink site-packages/external -> <fork>/external).
- Likelihood data: `cobaya-install` per likes_install.yaml (DESI DR2 BAO, Pantheon+,
  Planck 2018 lowl.TT + plik_lite TTTEEE via clipy; validated to 4e-9 on its test point).
- cobaya must load the fork's classy: `path: global` in the theory block (already in the yamls).
- SEDE runs REQUIRE z_gr_smg = 5 (see manuscript section 3) and the spectral guard likelihood.
- Chain integrity: src/CHECKSUMS.txt (sha256 of the production chains and minima).
