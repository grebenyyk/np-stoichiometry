# Stoichiometry of a saturated fluorite nanoparticle

A self-contained lattice-counting study of the *A*:*B* ratio in an idealised
fluorite (CaF2, AB2) nanoparticle.  The model assumes that every cation A,
including surface cations, keeps its full eightfold coordination by anions B:
no surface reconstruction, charge relaxation, or chemistry beyond the fixed
fluorite adjacency graph.

This repository also contains a multi-structure extension — rock-salt, CsCl,
zincblende, rutile — in [axby-np-stoichiometry/](axby-np-stoichiometry/); see
its [COMPARISON.md](axby-np-stoichiometry/COMPARISON.md).

The useful output is an exact set of formulas for saturated axis-aligned boxes,
especially the `n x n x n` nanocube.  The resulting integer sequences are elementary, but they have a clear crystallographic interpretation.

## Key results

1. **Fluorite recast.**  The B sublattice is simple cubic (spacing `a/2`).
   A sites are the body centres of one checkerboard colour of the B-cubes, so
   `CN(A)=8` and `CN(B)=4`.

2. **Exact B-rich identity.**  For any finite saturated A-cluster `S`,

   ```text
   N_B - 2 N_A = (1/4) sum_B (4 - d_B) >= 0.
   ```

   Here `d_B` is the number of selected A-cubes incident to B.  Every finite
   saturated particle is therefore B-rich; the excess is a boundary term.

3. **Form factor matters.**  `N_B` is not determined by `N_A` alone.  For
   example, both a `2 x 2 x 2` cube and a `4 x 2 x 1` slab have `N_A=4`, but
   their B-counts are `23` and `26`.

4. **Exact axis-aligned box formula.**  For an `a x b x c` box of B-cubes,

   ```text
   N_A = ceil(abc/2)
   N_B = (a+1)(b+1)(c+1) - U
   N_B - 2 N_A = (ab+ac+bc) + (a+b+c) + (1 - U - eps)
   ```

   where `eps=1` when `a,b,c` are all odd and `0` otherwise, and `U` is a
   parity-dependent count of enclosing-box vertices not touched by any selected
   A-cube.  The leading `ab+ac+bc` term is the sum of the three coordinate
   projected areas for this axis-aligned box family only.

5. **Scope caveat.**  The project does not claim a universal projected-area
   theorem, sphere/ellipsoid coefficient, or Wulff shape.  Those tempting
   continuum formulas ignore the fact that A occupies only one checkerboard
   colour of B-cubes; oblique and curved boundaries need separate analysis.

6. **Cubic-shell sequences / OEIS.**  With the checkerboard phase chosen so the
   origin cube is included,

   ```text
   N_A(n) = ceil(n^3/2)
   N_B(n) = (n+1)^3 - 2(1 + (-1)^n)
   excess(n) = 3n(n+1) - (3/2)(1 + (-1)^n)
   ```

   `N_A` and the odd-shell subsequences are already in the OEIS.  Direct OEIS
   searches currently return no match for the full `N_B` and full `excess`
   sequences; see [OEIS.md](OEIS.md) for the OEIS sequence identifications.

## Layout

```text
np-stoichiometry/
├── README.md                   # this file (fluorite study)
├── OEIS.md                     # OEIS sequence identifications
├── writeup/
│   ├── main.tex                # LaTeX writeup of the exact counting result
│   └── figures/                # generated plots for exact box families
├── scripts/
│   ├── fluorite.py             # core library: counts, closed forms, OEIS data
│   ├── verify_sequences.py     # sequence formulas vs brute force
│   ├── verify_box_formulas.py  # exact box decomposition and scope checks
│   ├── oeis_check.py           # live OEIS API lookup
│   ├── plot_ratio.py           # generate figures for axis-aligned families
│   └── run_all.sh              # run checks + figures + data table
├── data/
│   └── sequences.csv           # n, N_A, N_B, excess, ... for n = 1..30
└── axby-np-stoichiometry/      # multi-structure extension (NaCl, CsCl, ZnS, TiO2)
    ├── COMPARISON.md           # full writeup across structure types
    └── scripts/                # enumerators + closed forms for 5 structures
```

## Reproduce

```bash
cd scripts
bash run_all.sh
python3 verify_sequences.py
python3 verify_box_formulas.py
python3 oeis_check.py          # needs internet

# build the PDF:
cd ../writeup && pdflatex main.tex && pdflatex main.tex
```

Requires Python 3 (stdlib only, plus `matplotlib` for the figures) and a TeX
distribution for the writeup.
