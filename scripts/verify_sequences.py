"""
verify_sequences.py -- check that the closed forms match brute-force
enumeration, and that the cubic-shell sequences match the claimed OEIS data.

Run:  python3 verify_sequences.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from fluorite import (cube_cluster, box_cluster,
                      NA_cube_cf, NB_cube_cf, excess_cube_cf,
                      NA_box_cf, NB_box_cf, excess_box_cf, U_box, _epsilon,
                      cubic_shell_table, OEIS_HITS, OEIS_NOVEL)

PASS, FAIL = 0, 0

def check(cond, msg):
    global PASS, FAIL
    if cond:
        PASS += 1
    else:
        FAIL += 1
        print("  FAIL:", msg)

print("=" * 72)
print("1. Cubic shells: brute force vs closed form (n = 1..40)")
print("=" * 72)
for n in range(1, 41):
    na, nb, e = cube_cluster(n)
    check(na == NA_cube_cf(n), f"N_A n={n}: bf={na} cf={NA_cube_cf(n)}")
    check(nb == NB_cube_cf(n), f"N_B n={n}: bf={nb} cf={NB_cube_cf(n)}")
    check(e == excess_cube_cf(n), f"excess n={n}: bf={e} cf={excess_cube_cf(n)}")
print(f"  checked n=1..40  ({PASS} ok so far)")

print()
print("=" * 72)
print("2. General boxes: brute force vs closed form")
print("=" * 72)
p0 = PASS
for a in range(1, 14):
    for b in range(1, 14):
        for c in range(1, 14):
            na, nb, e = box_cluster(a, b, c)
            check(na == NA_box_cf(a, b, c), f"N_A box {a}x{b}x{c}")
            check(nb == NB_box_cf(a, b, c), f"N_B box {a}x{b}x{c}")
            check(e == excess_box_cf(a, b, c), f"excess box {a}x{b}x{c}")
print(f"  checked all a,b,c in 1..13  ({PASS - p0} identities)")

print()
print("=" * 72)
print("3. Bond-counting identity:  excess = (1/4) sum_B (4 - d_B)  [exact]")
print("=" * 72)
# recompute d_B for a cube by counting, for each B, how many in-box A-cubes
# share it.  Interior B has d_B=4; surface B has d_B<4.
def excess_via_bonds(a, b, c):
    A = set(black := __import__('fluorite').black_cubes_in_box(a, b, c))
    from collections import defaultdict
    d = defaultdict(int)
    for cube in black:
        for corner in __import__('fluorite').cube_corners(cube):
            d[corner] += 1
    return sum(4 - v for v in d.values()) // 4, len(d)
for (a, b, c) in [(2, 2, 2), (4, 2, 1), (3, 3, 3), (5, 4, 3)]:
    na, nb, e = box_cluster(a, b, c)
    e_bond, nb_bond = excess_via_bonds(a, b, c)
    check(e == e_bond, f"bond identity {a}x{b}x{c}: excess={e} bond={e_bond}")
    check(nb == nb_bond, f"B-count {a}x{b}x{c}: {nb} vs {nb_bond}")

print()
print("=" * 72)
print("4. Non-universality: same N_A, different N_B")
print("=" * 72)
na1, nb1, _ = box_cluster(2, 2, 2)
na2, nb2, _ = box_cluster(4, 2, 1)
print(f"  2x2x2 cube : N_A={na1}  N_B={nb1}")
print(f"  4x2x1 slab : N_A={na2}  N_B={nb2}")
check(na1 == na2 and nb1 != nb2, "same N_A must give different N_B")

print()
print("=" * 72)
print("5. OEIS data: cubic-shell sequences vs published first terms")
print("=" * 72)
rows = cubic_shell_table(20)
NA_seq = ",".join(str(r[1]) for r in rows)
NB_seq = ",".join(str(r[2]) for r in rows)
E_seq = ",".join(str(r[3]) for r in rows)
# A036486 published data (offset 0): 0,1,4,14,32,63,108,172,256,365,500,666,...
A036486 = "0," + NA_seq
check(A036486.startswith("0,1,4,14,32,63,108,172,256,365,500,666,864,1099"),
      "N_A matches A036486")
print(f"  N_A    = {NA_seq[:80]}...")
print(f"  N_B    = {NB_seq[:80]}...")
print(f"  excess = {E_seq[:80]}...")

print()
print("=" * 72)
print(f"RESULT: {PASS} checks passed, {FAIL} failed")
print("=" * 72)
sys.exit(1 if FAIL else 0)
