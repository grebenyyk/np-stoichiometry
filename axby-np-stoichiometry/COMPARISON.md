# Saturate-A stoichiometry across binary structure types

A companion to the fluorite-only study in the parent directory (`..`).  Same
disciplined scope: **exact axis-aligned box / cube counting only**, no
projected-area theorem, no Wulff / sphere / ellipsoid claims.  Closed forms are
stated only where verified against brute force; the one case expected to
resist a closed form (rutile) does not, and that is reported honestly.

Structures studied: **rock-salt (NaCl), CsCl, zincblende (ZnS)**
(cubic AB, two-sublattice) and **rutile (TiO₂)** (tetragonal AB₂, multi-site),
with **fluorite** carried as the reference.

## 1.  The model and a universal identity

For a structure AB_y, **saturate A** = keep every A's full B-coordination;
B-set = ⋃(A's B-neighbours).  `excess = N_B − y·N_A`.  The fluorite bond
identity generalises verbatim:

> **excess = (1/CN(B)) · Σ_B (CN(B) − d_B)  ≥  0**,   d_B = #in-cluster A-neighbours of B.

So **every** structure is B/anion-rich under saturate-A — the relation to
fluorite is structural, not incidental.  Verified by the bond count for all five.

## 2.  Unifying picture (cubic two-sublattice family)

In B-lattice units each cubic type is *"A on Z³ (full or one checkerboard
colour), B = A + V"* — they differ only by (A-density, neighbour set V):

| structure   | A-density | V (B-neighbours of A)        | CN(A) | CN(B) | y |
|-------------|-----------|------------------------------|-------|-------|---|
| fluorite    | ½ (even)  | 8 cube corners (±1,±1,±1)    | 8     | 4     | 2 |
| CsCl        | 1 (full)  | 8 cube corners (±1,±1,±1)    | 8     | 8     | 1 |
| NaCl        | ½ (even)  | 6 axis (±1,0,0)…             | 6     | 6     | 1 |
| zincblende  | ½ (even)  | 4 tetrahedron (even # of −)  | 4     | 4     | 1 |

(NaCl is the one where A and B live on the *same* lattice, different colours;
the other three have B at body-centre interstitials.)  Rutile does not fit this
table (Ti on 2 sublattices, O on 4 general-position sites) and is treated
separately.

## 3.  Closed forms — cube, verified

All five admit a clean closed form for the n×n×n saturated cube (verified:
fluorite/CsCl/NaCl/zincblende boxes to 11–16 cubed, rutile cube to n=20;
**21 583 brute-vs-formula checks pass, 0 fail**):

| structure   | N_A            | N_B                          | excess                      | period |
|-------------|----------------|------------------------------|-----------------------------|--------|
| fluorite    | ⌈n³/2⌉         | (n+1)³ − U(parity)           | 3n²+3n + parity             | 2      |
| CsCl        | n³             | (n+1)³                       | 3n²+3n+1                    | 1      |
| NaCl        | ⌈n³/2⌉         | ⌊n³/2⌋ + Σ face-pairs        | 3n² + 2·[n odd]             | 2      |
| zincblende  | ⌈n³/2⌉         | ⌊(n+1)³/2⌋                   | 3n(n+1)/2                   | 1      |
| rutile      | 2n³            | 4n³+8n²−n                    | 8n²−n                       | 1      |

General axis-aligned boxes are also tractable: CsCl exact (`(a+1)(b+1)(c+1)`);
NaCl and zincblende pick up a small parity "corner" term (NaCl: face-pair
parity; zincblende: `corner = 0` if a,b,c share parity else `−3`,
`excess = (ab+ac+bc+a+b+c+corner)/2`).  The rutile *general* tetragonal box
(a×b×c) is **not** derived here — only its cube.

**Answer to "are they all trivial with closed forms?": for axis-aligned cubes,
yes — all five, including the multi-site rutile.**  The expected hard case
(rutile, O on general positions with parameter u≈0.305) is *not* hard: for
generic u the O sites never coincide, so the count depends only on the fixed
Ti–O neighbour topology, not on u, and reduces to `8n²−n`.  One might expect
incommensurability to force a large quasi-polynomial period; it does not,
because counting a union of translated finite sets needs no commensurability.

## 4.  Relations to fluorite (the explicit mapping)

- **NaCl ≈ fluorite.** Same `N_A = ⌈n³/2⌉`, same leading excess `3n²`, same
  asymptotic coefficient (c ≈ 4.762).  Differs *only* in the edge term:
  fluorite `3n²+3n` (the +3n comes from the 8 *diagonal* neighbours inflating
  edges), NaCl `3n²` (the 6 *axis* neighbours inflate faces only — **no edge
  term**).  NaCl is "fluorite without the edge term."
- **CsCl = fluorite with the checkerboard filled in.** Same V (8 corners), but
  A full (density 1) vs fluorite's ½.  So `N_A = n³` (vs ⌈n³/2⌉) and
  `N_B = (n+1)³` with **no parity correction U** (fluorite had `(n+1)³−U`).
- **Zincblende = half of fluorite's B-shell.** V is 4 of fluorite's 8 vectors
  (one tetrahedron).  Its `N_B = ⌊(n+1)³/2⌋ = A036487`, the **floor companion**
  of fluorite's `N_A = ⌈n³/2⌉ = A036486` (the two sequences are cross-referenced
  in the OEIS).  Its excess `3n(n+1)/2` is half of fluorite's leading+edge,
  with no parity term.
- **Rutile** stands apart (tetragonal, AB₂, multi-site) but still yields a
  period-1 cube form.

## 5.  Asymptotics — same exponent, different coefficient

Every structure has `excess ∝ N_A^{2/3}` (a 3-D boundary term), so

> **N_B / N_A  =  y  +  c · N_A^{−1/3}  +  o(N_A^{−1/3})**,   c structure-dependent.

| structure   | y | surface coeff (excess/n²) | A-density | c = surf/ρ^{2/3} |
|-------------|---|---------------------------|-----------|------------------|
| fluorite    | 2 | 3                         | ½         | 3·2^{2/3} ≈ 4.762 |
| NaCl        | 1 | 3                         | ½         | 3·2^{2/3} ≈ 4.762 |
| CsCl        | 1 | 3                         | 1         | 3                |
| zincblende  | 1 | 1.5                       | ½         | 1.5·2^{2/3} ≈ 2.381 |
| rutile      | 2 | 8                         | 2         | 4·2^{1/3} ≈ 5.040 |

**Same N^{−1/3} exponent as fluorite throughout; the coefficient differs.**
NaCl coincides with fluorite exactly; CsCl is more stoichiometric (higher
A-density); zincblende the most stoichiometric (CN=4, fewest B stick out);
rutile the most B-rich.  (Rutile's "n×n×n" is in *cell* units — physically a
tetragonal prism since a≠c.)

## 6.  OEIS status (live-checked)

| structure   | N_A | N_B | excess |
|-------------|-----|-----|--------|
| fluorite    | A036486 | novel | novel |
| CsCl        | A000578 (n³) | A000578 ((n+1)³) | **A003215** (centered hexagonal) |
| NaCl        | A036486 | **novel** | **novel** |
| zincblende  | A036486 | **A036487** (⌊n³/2⌋) | A045943 (3× triangular) |
| rutile      | (2n³) | **novel** | **A139274** (n(8n−1)) |

The excess sequences are mostly already catalogued (CsCl→hexagonal,
zincblende→triangular, rutile→n(8n−1)); only **NaCl excess** is absent.
Three sequences are genuinely novel: **NaCl N_B**, **NaCl excess**, **rutile
N_O** (= 4n³+8n²−n).  NaCl is the most common structure type, which makes its
two novel sequences the natural submission candidates — but, as in the fluorite
project, these are elementary parity/quasi-polynomial sequences whose value
rests on the crystallographic interpretation, not on the numbers themselves.

## 7.  Candidate observations worth flagging (honest)

These are structural observations from the counting, **not** theorems, and not
deep discoveries — but they are the most interesting things to come out:

1. **The (A-density, V) unification** of the four cubic types, with fluorite/
   NaCl/zincblende sharing `A = ⌈n³/2⌉`.
2. **The ceil/floor duality**: fluorite `N_A = ⌈n³/2⌉` (A036486) and zincblende
   `N_B = ⌊(n+1)³/2⌋` (A036487) are companion OEIS sequences — a tidy link
   between the cation count of one structure and the anion count of another.
3. **NaCl has no edge term** (axis-6 neighbours) where fluorite has +3n
   (diagonal-8): the neighbour *geometry* controls the subleading term.
4. **Rutile is clean despite O on general positions** — generic-position
   incommensurability does not obstruct a closed count.  This contradicts the
   naive expectation that multi-site / general-position
   structures would resist closed forms.

## 8.  What is NOT claimed

- No universal projected-area / Wulff / sphere theorem (those were removed from
  the fluorite project and are not reinstated here).
- Closed forms are verified only for **axis-aligned boxes/cubes** (and the
  rutile *cube* specifically); oblique or curved boundaries are out of scope.
- The asymptotic coefficients are observed leading terms, not proven bounds.
- Rutile depends on the standard neighbour topology (u≈0.305); a different
  topology (e.g. a phase transition in u) would change the count.

## Reproduce

```bash
cd scripts
python3 verify_structures.py        # 21 583 checks pass
python3 oeis_check_structures.py    # live OEIS (needs internet)
python3 plot_structures.py          # ratio-vs-N_A figure
python3 structures.py               # brute-force tables
```
