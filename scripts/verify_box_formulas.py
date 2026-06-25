"""
verify_box_formulas.py -- verify the exact axis-aligned box formulas.

This replaces the earlier projected-area/Wulff check.  The central supported
claim is finite and axis-aligned:

    excess = N_B - 2 N_A
           = (ab+ac+bc) + (a+b+c) + (1 - U - eps)

for an a x b x c box of B-cubes.  The final section also prints a small
oblique-boundary example showing why this file does not claim a universal
projected-area theorem.

Run:  python3 verify_box_formulas.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from fluorite import (box_cluster, cluster_counts, NA_box_cf, NB_box_cf,
                      excess_box_cf, U_box, _epsilon, C_CUBE,
                      asymptotic_ratio)

PASS, FAIL = 0, 0

def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
    else:
        FAIL += 1
        print("  FAIL:", msg)

print("=" * 72)
print("1. Sample boxes: exact decomposition")
print("=" * 72)
print(f"{'box':>10} {'excess':>7} {'ab+ac+bc':>9} {'a+b+c':>7} {'corner':>7}")
for (a, b, c) in [(2, 2, 2), (4, 2, 1), (3, 3, 3), (5, 4, 3),
                  (6, 6, 6), (7, 5, 3), (1, 1, 1), (2, 3, 4)]:
    na, nb, excess = box_cluster(a, b, c)
    surface = a * b + a * c + b * c
    edge = a + b + c
    corner = 1 - U_box(a, b, c) - _epsilon(a, b, c)
    check(excess == surface + edge + corner, f"decomposition {a}x{b}x{c}")
    print(f"  {a}x{b}x{c:<3} {excess:>7} {surface:>9} {edge:>7} {corner:>7}")

print()
print("=" * 72)
print("2. Exhaustive box formulas for a,b,c in 1..13")
print("=" * 72)
p0 = PASS
for a in range(1, 14):
    for b in range(1, 14):
        for c in range(1, 14):
            na, nb, excess = box_cluster(a, b, c)
            check(na == NA_box_cf(a, b, c), f"N_A {a}x{b}x{c}")
            check(nb == NB_box_cf(a, b, c), f"N_B {a}x{b}x{c}")
            check(excess == excess_box_cf(a, b, c), f"excess {a}x{b}x{c}")
print(f"  checked {PASS - p0} identities")

print()
print("=" * 72)
print("3. Cube asymptotic: N_B/N_A = 2 + c N_A^(-1/3) + o")
print("=" * 72)
print(f"  c_cube = 3*2^(2/3) = {C_CUBE:.4f}")
for n in [20, 50, 100]:
    na, nb, _ = box_cluster(n, n, n)
    actual = nb / na
    pred = asymptotic_ratio(C_CUBE, na)
    print(f"  cube n={n:>3}: N_A={na:>6}  actual={actual:.5f}  pred={pred:.5f}")
check(abs(box_cluster(100, 100, 100)[1] / box_cluster(100, 100, 100)[0]
          - asymptotic_ratio(C_CUBE, box_cluster(100, 100, 100)[0])) < 0.01,
      "cube ratio matches asymptotic coefficient at n=100")

print()
print("=" * 72)
print("4. Scope check: oblique boundary is not covered by the box formula")
print("=" * 72)

def tetra_cluster(m):
    sites = []
    for i in range(m + 1):
        for j in range(m + 1 - i):
            for k in range(m + 1 - i - j):
                if (i + j + k) % 2 == 0:
                    sites.append((i, j, k))
    return cluster_counts(sites)

print(f"{'m':>4} {'N_A':>8} {'excess':>8} {'excess/m^2':>12} {'naive proj coeff':>16}")
for m in [40, 80, 100]:
    na, nb, excess = tetra_cluster(m)
    # The naive full-cube projected-area argument would give 1.5*m^2 here.
    print(f"{m:>4} {na:>8} {excess:>8} {excess/(m*m):>12.4f} {1.5:>16.4f}")
check(abs(tetra_cluster(100)[2] / (100 * 100) - 1.5) > 0.1,
      "oblique tetrahedron differs from naive projected-area coefficient")

print()
print("=" * 72)
print(f"RESULT: {PASS} checks passed, {FAIL} failed")
print("=" * 72)
sys.exit(1 if FAIL else 0)