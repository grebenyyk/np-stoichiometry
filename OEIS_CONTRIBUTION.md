# OEIS contribution strategy

## Candidate sequences

For the `n x n x n` saturated fluorite nanocube, `n = 1, 2, 3, ...`, with the
checkerboard phase chosen so that the origin cube is included:

| # | sequence | closed form | first terms |
| --- | --- | --- | --- |
| 1 | **B-count, full** | `(n+1)^3 - 2(1 + (-1)^n)` | 8, 23, 64, 121, 216, 339, 512, 725, 1000, 1327, ... |
| 2 | **B-excess, full** | `3n(n+1) - (3/2)(1 + (-1)^n)` | 6, 15, 36, 57, 90, 123, 168, 213, 270, 327, ... |
| 3 | B-count, even shells | `(2m+1)^3 - 4` | 23, 121, 339, 725, 1327, 2193, ... |
| 4 | B-excess, even shells | `3(n^2+n-1)` at `n=2m` | 15, 57, 123, 213, 327, ... |

Direct OEIS searches on 2026-06-25 returned no match for these four sequences.
That should be treated as a practical cataloguing fact, not as mathematical
novelty.  Re-run `python3 scripts/oeis_check.py` before submitting.

Related known entries:

- `N_A` full: [A036486](https://oeis.org/A036486), `ceil(n^3/2)`.
- `N_A` odd shells: [A050492](https://oeis.org/A050492), thickened cube numbers.
- `N_B` odd shells: [A016743](https://oeis.org/A016743), even cubes `(2m)^3`.
- B-excess odd shells: [A152746](https://oeis.org/A152746), six times hexagonal numbers.

## Recommendation

**Submit the full B-excess sequence first.**  It is the cleanest entry: it says
exactly how many extra B atoms are required beyond the bulk `AB2` ratio in the
saturated nanocube model.  The formula is short, the interpretation is concrete,
and it cross-references naturally to A036486.

**Treat the full B-count as an optional companion.**  It is useful, but it is
also just `2*A036486 + excess`, so it may be viewed as derivative.  If the
excess entry is accepted, the B-count can either be a separate linked entry or a
formula/comment on the same submission, depending on editor feedback.

**Do not submit the even-shell subsequences separately.**  They are simply the
even halves of the two full sequences.  They add little independent value and
are more likely to be rejected as redundant.

## Frank assessment

These are thin sequences.  They are parity-modulated polynomials, not deep
number-theoretic objects.  Submitted bare, they are likely to look like near
cubes and near quadratics with an arbitrary parity correction.

They become worth submitting only because the definition is natural: count B
atoms, or B-excess over `2N_A`, in a saturated fluorite nanocube.  That is a
legitimate combinatorial/crystallographic interpretation, especially because
A036486 already records the matching cation count.

Expected acceptance is **modest, not guaranteed**.  A reasonable mental model is
roughly 40-60%, depending on how compactly the entry is written and whether the
editor sees the fluorite interpretation as useful enough.  The OEIS value is
findability, not mathematical depth.

## What to say

Use the exact finite-counting story:

- Saturated fluorite `AB2` nanocube of side `n` in B-lattice units.
- B sites form a simple-cubic lattice; A sites are body centres of one
  checkerboard colour of B-cubes.
- Every selected A-cube contributes all eight corner B sites; shared corners are
  counted once.
- `N_A(n)=ceil(n^3/2)` is A036486.
- The excess is `a(n)=N_B(n)-2N_A(n)`.
- For boxes, the exact decomposition is
  `N_B-2N_A=(ab+ac+bc)+(a+b+c)+(1-U-eps)`.

Also state the checkerboard phase convention.  For odd `n`, the opposite colour
inside the same enclosing `n x n x n` box gives a different finite cluster; this
project uses the phase containing the origin cube, matching the `ceil(n^3/2)`
convention.

## What not to say

Do not pitch this as a universal projected-area theorem, a sphere/ellipsoid
asymptotic, or a Wulff-shape result.  The exact result here is for finite
axis-aligned boxes.  Oblique or curved boundaries interact with the checkerboard
A-sublattice differently and need separate treatment.

Avoid strong words like "deep" or "new theorem" for the sequences.  A better
tone is: elementary sequence with a concrete crystallographic counting model.

## Draft OEIS fields

### B-excess entry, preferred

- **Name:** Excess number of B atoms over twice A atoms in a saturated fluorite
  `AB2` nanocube of side `n`.
- **Data:** `6, 15, 36, 57, 90, 123, 168, 213, 270, 327, 396, 465, 546, 627,
  720, 813, 918, 1023, 1140, 1257, 1386, 1515, 1656, 1797, 1950, 2103, 2268,
  2433, 2610, 2787`
- **Formula:** `a(n)=3*n*(n+1)-3*(1+(-1)^n)/2`, for `n>=1`.
- **Equivalently:** `a(2m-1)=6m(2m-1)`, `a(2m)=3(4m^2+2m-1)`.
- **Comment:** B sites form a simple-cubic lattice; A sites are body centres of
  one checkerboard colour of B-cubes.  The saturated nanocube contains all B
  corners incident to the selected A-cubes, counted once.
- **Crossrefs:** A036486, A050492, A016743, A152746, and the optional B-count
  companion if submitted.
- **Offset:** `1,1`.

### B-count entry, optional companion

- **Name:** Number of B atoms in a saturated fluorite `AB2` nanocube of side
  `n`.
- **Data:** `8, 23, 64, 121, 216, 339, 512, 725, 1000, 1327, 1728, 2193, 2744,
  3371, 4096, 4909, 5832, 6855, 8000, 9257, 10648, 12163, 13824, 15621, 17576,
  19679, 21952, 24385, 27000, 29787`
- **Formula:** `a(n)=(n+1)^3-2*(1+(-1)^n)`, for `n>=1`.
- **Comment:** With `b(n)=ceil(n^3/2)` from A036486 and excess `e(n)`, this is
  `a(n)=2*b(n)+e(n)`.
- **Offset:** `1,1`.

## Bottom line

Submit the excess sequence if you want the fluorite counting result findable.
Submit the B-count only as a companion or comment.  Skip the even-shell entries.
The exact box formula is the real content; the OEIS entries are a useful index
card for it.
