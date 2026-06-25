"""
plot_ratio.py -- generate figures for exact axis-aligned box families.

Outputs (in ../writeup/figures/):
  ratio_vs_NA.pdf   -- N_B / N_A vs N_A for several box aspect ratios.
  excess_vs_NA.pdf  -- excess vs N_A^(2/3), showing boundary-order growth.

Run:  python3 plot_ratio.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
import fluorite as F

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

FIG_DIR = os.path.join(os.path.dirname(__file__), "..", "writeup", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

FAMILIES = [
    (r"cube $m{\times}m{\times}m$", "C0", "o", lambda m: (m, m, m), range(1, 41)),
    (r"rod $m{\times}m{\times}4m$", "C3", "^", lambda m: (m, m, 4 * m), range(1, 21)),
    (r"slab $m{\times}4m{\times}4m$", "C2", "s", lambda m: (m, 4 * m, 4 * m), range(1, 16)),
]

series = []
for label, color, marker, dims, values in FAMILIES:
    rows = []
    for m in values:
        a, b, c = dims(m)
        na, nb, excess = F.box_cluster(a, b, c)
        rows.append((m, a, b, c, na, nb, excess))
    series.append((label, color, marker, dims, rows))

# --------------------------------------------------------------------------
# Figure 1: ratio vs N_A
# --------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6.2, 4.3))
for label, color, marker, dims, rows in series:
    NA = [row[4] for row in rows]
    ratio = [row[5] / row[4] for row in rows]
    ax.plot(NA, ratio, marker + "-", ms=3, lw=1.0, color=color, label=label)

max_NA = max(row[4] for _, _, _, _, rows in series for row in rows)
NA_line = np.linspace(5, max_NA, 400)
ax.plot(NA_line, F.asymptotic_ratio(F.C_CUBE, NA_line), ":", color="C0", lw=1.4,
        label=rf"cube asymptote, $2+{F.C_CUBE:.2f}\,N_A^{{-1/3}}$")
ax.axhline(2, color="k", lw=0.7, alpha=0.5)
ax.set_xlabel(r"$N_A$  (number of A atoms)")
ax.set_ylabel(r"$N_B / N_A$")
ax.set_xscale("log")
ax.set_ylim(1.9, 8.5)
ax.set_title("Saturated fluorite boxes")
ax.legend(loc="upper right", fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "ratio_vs_NA.pdf"))
print("wrote", os.path.join(FIG_DIR, "ratio_vs_NA.pdf"))

# --------------------------------------------------------------------------
# Figure 2: excess vs N_A^(2/3)
# --------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6.2, 4.3))
for label, color, marker, dims, rows in series:
    x = [row[4] ** (2 / 3) for row in rows]
    y = [row[6] for row in rows]
    ax.plot(x, y, marker + "-", ms=3, lw=1.0, color=color, label=label)

xx = np.linspace(1, max_NA ** (2 / 3), 200)
for alpha, beta, gamma, color in [(1, 1, 1, "C0"), (1, 1, 4, "C3"), (1, 4, 4, "C2")]:
    coeff = F.box_asymptotic_c(alpha, beta, gamma)
    ax.plot(xx, coeff * xx, ":", color=color, lw=1.2,
            label=rf"$c={coeff:.2f}$ for {alpha}:{beta}:{gamma}")

ax.set_xlabel(r"$N_A^{2/3}$")
ax.set_ylabel(r"excess  $= N_B - 2 N_A$")
ax.set_title(r"Box excess is boundary-order:  $\propto N_A^{2/3}$")
ax.legend(loc="upper left", fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIG_DIR, "excess_vs_NA.pdf"))
print("wrote", os.path.join(FIG_DIR, "excess_vs_NA.pdf"))