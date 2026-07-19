#!/usr/bin/env python3
"""Regenerate fsigma8_prediction.png from the committed fσ8 prediction data.

Reconstructed 2026-07-21: the original plotting script for this figure was never committed
(commit 9429568 in the hub added only the PNG + the results JSON). The computation pipeline is
research/hub/unified/sims/sede_fsigma8_isw_finitek.py, whose output is the JSON read here
(figures/data/fsigma8_prediction.json), so this reproduces the exact published curves and fixes
the panel-(b) title that overran the axis.

Run:  python make_fsigma8_figure.py   ->  paper/figures/fsigma8_prediction.png (+ figures/ copy)
"""
import os, json
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, "data", "fsigma8_prediction.json")))
z = d["z"]

fig, ax = plt.subplots(1, 2, figsize=(11, 4.0), constrained_layout=True)

# --- (a) fσ8(z) growth-rate history: ΛCDM vs SEDE at k = 0.1 h/Mpc ---
ax[0].plot(z, d["fs8_lcdm_k0p1"], color="0.35", lw=1.8, label=r"$\Lambda$CDM")
ax[0].plot(z, d["fs8_sede_k0p1"], color="#c0392b", lw=1.8, ls="--", label="SEDE")
ax[0].set_title("(a) growth-rate history")
ax[0].set_xlabel(r"$z$")
ax[0].set_ylabel(r"$f\sigma_8(z)$   ($k=0.1\,h\,\mathrm{Mpc}^{-1}$)")
ax[0].legend(loc="upper right", fontsize=9)

# --- (b) SEDE enhancement Δfσ8/fσ8 [%] — scale-independent (three k overlap) ---
styles = [("0.05", "#e69138", "-", 2.4), ("0.1", "#c0392b", "--", 2.0), ("0.2", "#7d3c98", ":", 2.0)]
for k, c, ls, lw in styles:
    ax[1].plot(z, [100.0 * r for r in d["ratio"][k]], color=c, ls=ls, lw=lw,
               label=rf"$k = {k}\,h\,\mathrm{{Mpc}}^{{-1}}$")
ax[1].axhspan(0.0, 1.0, color="0.85", alpha=0.6)
ax[1].text(0.55, 0.5, "DESI DR3 / Euclid\nper-bin sensitivity", fontsize=8.5, color="0.4",
           ha="left", va="center")
ax[1].set_title("(b) SEDE enhancement — scale-independent")
ax[1].set_xlabel(r"$z$")
ax[1].set_ylabel(r"$\Delta f\sigma_8/f\sigma_8$   [%]")
ax[1].set_ylim(0.0, 3.5)
ax[1].legend(loc="upper right", fontsize=9)

for out in [os.path.join(HERE, "..", "paper", "figures", "fsigma8_prediction.png"),
            os.path.join(HERE, "fsigma8_prediction.png")]:
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=170)
    print("wrote", os.path.normpath(out))
plt.close(fig)
