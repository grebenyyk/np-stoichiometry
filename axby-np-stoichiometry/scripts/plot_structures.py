"""
plot_structures.py -- ratio N_B/N_A vs N_A for all five structures (cubes),
using closed forms so we can go to large n.  Shows the common N^(-1/3)
approach to the bulk ratio y with structure-dependent coefficient c.

Output: figures/ratio_vs_NA.pdf
Run:    python3 plot_structures.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
import structures as S
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(FIG, exist_ok=True)

# (label, N_A_fn, N_B_fn, y, c, color)
specs = [
    ("fluorite",   S.fluorite_NA_cube,   S.fluorite_NB_cube,    2, S.C_FLUORITE,   "C0"),
    ("NaCl",       lambda n: S.nacl_NA(n,n,n),       lambda n: S.nacl_NB(n,n,n),       1, S.C_NACL,       "C1"),
    ("CsCl",       lambda n: S.cscl_NA(n,n,n),       lambda n: S.cscl_NB(n,n,n),       1, S.C_CSCL,       "C2"),
    ("zincblende", lambda n: S.zincblende_NA(n,n,n), lambda n: S.zincblende_NB(n,n,n), 1, S.C_ZINCBLENDE, "C3"),
    ("rutile",     S.rutile_NTl_cube,    S.rutile_NO_cube,      2, S.C_RUTILE,     "C4"),
]

ns = list(range(1, 121))
fig, ax = plt.subplots(figsize=(6.6, 4.5))
for name, NAf, NBf, y, c, col in specs:
    NA = [NAf(n) for n in ns]
    r = [NBf(n) / NAf(n) for n in ns]
    ax.plot(NA, r, "-", lw=1.1, color=col, label=f"{name}  ($y$={y}, $c$={c:.2f})")
NAline = np.linspace(5, max(NA), 400)
for name, NAf, NBf, y, c, col in specs:
    ax.plot(NAline, y + c * NAline ** (-1/3), ":", color=col, lw=0.9, alpha=0.6)
ax.axhline(2, color="k", lw=0.6, alpha=0.4)
ax.axhline(1, color="k", lw=0.6, alpha=0.4)
ax.set_xscale("log")
ax.set_ylim(0.95, 8.5)
ax.set_xlabel(r"$N_A$  (A atoms)")
ax.set_ylabel(r"$N_B / N_A$")
ax.set_title("Saturated nanoparticle stoichiometry: 5 structure types")
ax.legend(loc="upper right", fontsize=8)
fig.tight_layout()
out = os.path.join(FIG, "ratio_vs_NA.pdf")
fig.savefig(out)
print("wrote", out)
