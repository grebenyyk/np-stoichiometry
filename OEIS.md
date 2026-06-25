# OEIS connections

Lattice counts in this repository connect to the following OEIS sequences.
Identifications are verified by live lookup (`scripts/oeis_check.py`,
`axby-np-stoichiometry/scripts/oeis_check_structures.py`) and by brute-force
enumeration.

## Fluorite (CaF2, AB2) nanocube, side n

| quantity | formula | OEIS |
|---|---|---|
| N_A (cations) | ceil(n^3/2) | A036486 |
| N_A, odd shells (n = 2m-1) | 4m^3 - 6m^2 + 3m | A050492 (thickened cubes) |
| N_B, odd shells | (2m)^3 | A016743 (even cubes) |
| excess, odd shells | 6m(2m-1) | A152746 (6 x hexagonal) |
| N_B (full) | (n+1)^3 - 2(1 + (-1)^n) | not in OEIS |
| excess (full) | 3n(n+1) - (3/2)(1 + (-1)^n) | not in OEIS |

## Multi-structure extension (axby-np-stoichiometry/)

| structure | quantity | formula | OEIS |
|---|---|---|---|
| CsCl | N_B | (n+1)^3 | A000578 (cubes) |
| CsCl | excess | 3n^2 + 3n + 1 | A003215 (centered hexagonal) |
| NaCl | N_B | floor(n^3/2) + 3n^2 + 3*(n mod 2) | not in OEIS |
| NaCl | excess | 3n^2 + 1 - (-1)^n | not in OEIS |
| zincblende | N_B | floor((n+1)^3/2) | A036487 (floor(n^3/2)) |
| zincblende | excess | 3n(n+1)/2 | A045943 (3 x triangular) |
| rutile | N_O | 4n^3 + 8n^2 - n | not in OEIS |
| rutile | excess | 8n^2 - n | A139274 (n(8n-1)) |

Sequences marked "not in OEIS" returned no match on live lookup. The NaCl
excess sequence, `a(n) = 3n^2 + 1 - (-1)^n` (0, 5, 12, 29, 48, 77, 108, ...),
is the most direct candidate for a new entry; its first 55 terms and a
verified generating function are produced by `axby-np-stoichiometry/scripts/`.
