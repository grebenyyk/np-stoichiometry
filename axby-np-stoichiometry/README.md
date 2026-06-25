# Saturate-A stoichiometry across binary structure types

A companion to the fluorite study in the [parent directory](..).
Same question (excess B on the surface of a saturated nanoparticle), applied to
**rock-salt (NaCl), CsCl, zincblende (ZnS)**, and **rutile (TiO₂)**, with
fluorite carried as the reference.  Same disciplined scope: exact axis-aligned
box/cube counting only — no projected-area theorem, no Wulff / sphere claims.
Full writeup: [COMPARISON.md](COMPARISON.md).

## Headline results (all verified — 21 583 brute-vs-formula checks pass)

1. **Universal bond identity** (generalises fluorite):
   `excess = (1/CN(B))·Σ_B(CN(B)−d_B) ≥ 0` — every structure is B-rich under
   saturate-A.

2. **All five have clean cube closed forms** (the expected hard case, rutile,
   is clean too — O on general positions does not obstruct the count):

   | structure | excess(cube) | period |
   |---|---|---|
   | fluorite | 3n²+3n+parity | 2 |
   | CsCl | 3n²+3n+1 | 1 |
   | NaCl | 3n²+2·[n odd] | 2 |
   | zincblende | 3n(n+1)/2 | 1 |
   | rutile | 8n²−n | 1 |

3. **Same asymptotics as fluorite** in exponent (`N_B/N_A = y + c·N_A^{−1/3}`),
   different coefficient: fluorite=NaCl (c≈4.76), CsCl (3), zincblende (≈2.38),
   rutile (≈5.04).

4. **OEIS**: most sequences are known (CsCl excess→A003215, zincblende
   N_B→A036487, zincblende excess→A045943, rutile excess→A139274).  **Novel**:
   NaCl N_B, NaCl excess, rutile N_O.

5. **Relations to fluorite**: NaCl = fluorite without the edge term; CsCl =
   fluorite with the checkerboard filled; zincblende N_B = ⌊(n+1)³/2⌋ (the
   floor-companion A036487 of fluorite's N_A = A036486).

## Layout

```
axby-np-stoichiometry/
├── COMPARISON.md                 # full writeup (relations, forms, OEIS, asymptotics)
├── README.md
├── figures/ratio_vs_NA.pdf       # 5-structure ratio plot
└── scripts/
    ├── structures.py             # enumerators + closed forms (5 structures)
    ├── verify_structures.py      # brute vs closed form + asymptotics + OEIS seqs
    ├── oeis_check_structures.py  # live OEIS lookup
    └── plot_structures.py        # ratio-vs-N_A figure
```

## Reproduce

```bash
cd scripts
python3 verify_structures.py        # 21 583 checks pass, 0 fail
python3 oeis_check_structures.py    # live OEIS (needs internet)
python3 plot_structures.py
```
Requires Python 3 (stdlib + matplotlib).
