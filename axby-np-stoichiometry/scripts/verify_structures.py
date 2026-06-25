"""
verify_structures.py -- verify closed forms vs brute force for all five
structures, confirm the rutile cube formula 8 n^2 - n to large n, read off
asymptotic coefficients, and print OEIS-ready cube sequences.

Run:  python3 verify_structures.py
"""
import os, sys, math
sys.path.insert(0, os.path.dirname(__file__))
import structures as S

PASS, FAIL = 0, 0
def check(cond, msg):
    global PASS, FAIL
    if cond: PASS += 1
    else: FAIL += 1; print("  FAIL:", msg)

print("=" * 72)
print("1. CsCl:  N_A=abc, N_B=(a+1)(b+1)(c+1), excess=(a+1)(b+1)(c+1)-abc")
print("=" * 72)
p0 = PASS
for a in range(1, 12):
    for b in range(1, 12):
        for c in range(1, 12):
            na, nb, e = S.cscl_box(a, b, c)
            check(na == S.cscl_NA(a, b, c), f"CsCl N_A {a}x{b}x{c}")
            check(nb == S.cscl_NB(a, b, c), f"CsCl N_B {a}x{b}x{c}")
            check(e == S.cscl_excess(a, b, c), f"CsCl excess {a}x{b}x{c}")
print(f"  CsCl boxes 1..11 cubed: ok ({PASS - p0} identities)")

print()
print("=" * 72)
print("2. NaCl:  face-pair formula (corrected for dimension parity)")
print("=" * 72)
p0 = PASS
for a in range(1, 13):
    for b in range(1, 13):
        for c in range(1, 13):
            na, nb, e = S.nacl_box(a, b, c)
            check(na == S.nacl_NA(a, b, c), f"NaCl N_A {a}x{b}x{c}")
            check(nb == S.nacl_NB(a, b, c), f"NaCl N_B {a}x{b}x{c}")
            check(e == S.nacl_excess(a, b, c), f"NaCl excess {a}x{b}x{c}")
print(f"  NaCl boxes 1..12 cubed: ok ({PASS - p0} identities)")

print()
print("=" * 72)
print("3. Zincblende:  excess = (ab+ac+bc + a+b+c + corner)/2,")
print("   corner = 0 if a,b,c share parity else -3")
print("=" * 72)
p0 = PASS
for a in range(1, 17):
    for b in range(1, 17):
        for c in range(1, 17):
            na, nb, e = S.zincblende_box(a, b, c)
            check(na == S.zincblende_NA(a, b, c), f"zb N_A {a}x{b}x{c}")
            check(nb == S.zincblende_NB(a, b, c), f"zb N_B {a}x{b}x{c}")
            check(e == S.zincblende_excess(a, b, c), f"zb excess {a}x{b}x{c}")
# cube formula
for n in range(1, 30):
    na, nb, e = S.zincblende_box(n, n, n)
    check(e == S.zincblende_excess_cube(n), f"zb cube excess n={n}")
    check(nb == S.zincblende_NB_cube(n), f"zb cube N_B n={n}")
print(f"  Zincblende boxes 1..16 cubed + cubes 1..29: ok ({PASS - p0} identities)")

print()
print("=" * 72)
print("4. Rutile cube:  excess = 8 n^2 - n  (verified to large n)")
print("=" * 72)
p0 = PASS
for n in range(1, 21):
    na, nb, e = S.rutile_box(n, n, n)
    check(na == S.rutile_NTl_cube(n), f"rutile N_Ti n={n}: {na} vs {S.rutile_NTl_cube(n)}")
    check(e == S.rutile_excess_cube(n), f"rutile excess n={n}: {e} vs {S.rutile_excess_cube(n)}")
    check(nb == S.rutile_NO_cube(n), f"rutile N_O n={n}: {nb} vs {S.rutile_NO_cube(n)}")
print(f"  Rutile cubes n=1..20: excess = 8 n^2 - n exactly  ({PASS - p0} identities)")
print("  (general tetragonal box a x b x c is NOT pursued as a closed form here;")
print("   the cube already shows the O 'general position' does not block a clean count.)")

print()
print("=" * 72)
print("5. Asymptotic ratio  N_B/N_A = y + c * N_A^(-1/3)  (cubes)")
print("=" * 72)
print(f"  {'structure':>11} {'y':>3} {'surf':>5} {'rho':>5} {'c':>6}   {'n=20':>9} {'n=50':>9} {'n=100':>9}")
def ratio_at(boxfn, n):
    na, nb, _ = boxfn(n, n, n); return nb / na
specs = [
    ("fluorite",   S.fluorite_box,   2, S.C_FLUORITE,   3.0, 0.5),
    ("CsCl",       S.cscl_box,       1, S.C_CSCL,       3.0, 1.0),
    ("NaCl",       S.nacl_box,       1, S.C_NACL,       3.0, 0.5),
    ("zincblende", S.zincblende_box, 1, S.C_ZINCBLENDE, 1.5, 0.5),
    ("rutile",     S.rutile_box,     2, S.C_RUTILE,     8.0, 2.0),
]
for name, fn, y, c, surf, rho in specs:
    print(f"  {name:>11} {y:>3} {surf:>5} {rho:>5} {c:>6.3f}   "
          f"{ratio_at(fn,20):>9.4f} {ratio_at(fn,50):>9.4f} {ratio_at(fn,100):>9.4f}")
print("  (c = surf / rho^(2/3). All approach the bulk ratio y from above; same")
print("   N^(-1/3) exponent, different coefficient c.)")

print()
print("=" * 72)
print("6. OEIS-ready cube sequences (n = 1..12)")
print("=" * 72)
def seq(fn, idx=1):
    return ",".join(str(fn(n, n, n)[idx]) for n in range(1, 13))
def seqc(fn):
    return ",".join(str(fn(n)) for n in range(1, 13))
print(f"  fluorite   N_A   = {seq(S.fluorite_box, 0)}")
print(f"  fluorite   N_B   = {seq(S.fluorite_box, 1)}")
print(f"  fluorite   exc   = {seq(S.fluorite_box, 2)}")
print(f"  CsCl       N_A   = {seq(S.cscl_box, 0)}   (= n^3)")
print(f"  CsCl       N_B   = {seq(S.cscl_box, 1)}   (= (n+1)^3)")
print(f"  CsCl       exc   = {seq(S.cscl_box, 2)}   (= 3n^2+3n+1)")
print(f"  NaCl       N_A   = {seq(S.nacl_box, 0)}")
print(f"  NaCl       N_B   = {seq(S.nacl_box, 1)}")
print(f"  NaCl       exc   = {seq(S.nacl_box, 2)}")
print(f"  zincblende N_A   = {seq(S.zincblende_box, 0)}")
print(f"  zincblende N_B   = {seq(S.zincblende_box, 1)}")
print(f"  zincblende exc   = {seq(S.zincblende_box, 2)}   (= 3n(n+1)/2)")
print(f"  rutile     N_Ti  = {seq(S.rutile_box, 0)}   (= 2n^3)")
print(f"  rutile     N_O   = {seq(S.rutile_box, 1)}")
print(f"  rutile     exc   = {seq(S.rutile_box, 2)}   (= 8n^2 - n)")

print()
print("=" * 72)
print(f"RESULT: {PASS} checks passed, {FAIL} failed")
print("=" * 72)
sys.exit(1 if FAIL else 0)
