"""
structures.py -- saturate-A stoichiometry for several binary A_x B_y structure
types, as exact axis-aligned box / cube counts.

This is a companion to the fluorite-only project.  The same disciplined scope
applies: we ONLY do exact lattice counting for axis-aligned boxes (and, for
rutile, tetragonal boxes).  No projected-area theorem, no Wulff / sphere /
ellipsoid claims.  Closed forms are stated only when verified against brute
force; anything that resists a simple form is reported as such.

================================================================
The unifying picture (cubic two-sublattice family)
================================================================
In B-lattice units, every cubic structure here is

        A-sites  =  a sublattice of Z^3   (full, or one checkerboard colour)
        B-sites  =  A + V                 (V = the B-neighbour vectors of A)

so the four cubic types differ only by (A-density, V):

  structure   A-density   V (B-neighbours of A)            CN(A) CN(B)  y
  ---------   ----------  -------------------------------  ----- ----- ---
  fluorite    1/2 (even)  8 cube corners  ( +/-1,+/-1,+/-1 )   8    4    2
  CsCl        1   (full)  8 cube corners  ( +/-1,+/-1,+/-1 )   8    8    1
  NaCl        1/2 (even)  6 axis          ( +/-1,0,0 etc )     6    6    1
  zincblende  1/2 (even)  4 tetrahedron   ( even # of minuses) 4    4    1

(We scale the lattice by 2 so all coordinates are integers: A lives on
(2i,2j,2k)-type points, B on odd-offset points.  This changes nothing in the
counts.)

Saturation model: given a finite A-set S (an a x b x c box of cells),
        B(S) = union over A in S of (A + V),
i.e. every A keeps its full B-coordination.  excess = N_B - y * N_A.

General bond identity (holds for ALL structures, same form as fluorite):
        excess = (1/CN(B)) * sum_B (CN(B) - d_B)  >= 0,
so every structure is B/anion-rich under saturate-A.

Rutile (TiO2, AB2) is treated separately: Ti on two sublattices, O on four
general-position sites (parameter u ~ 0.305).  Its counts do NOT reduce to a
low-period quasi-polynomial -- the demonstration that tractability is not
universal.
"""

from itertools import product
from fractions import Fraction

# --------------------------------------------------------------------------
# Generic saturate-A counter
# --------------------------------------------------------------------------

def saturate_count(A_sites, V, y):
    """Given A-sites (list of integer tuples) and neighbour vectors V,
    return (N_A, N_B, excess = N_B - y*N_A).  B = union of A+V."""
    B = set()
    for a in A_sites:
        for v in V:
            B.add((a[0] + v[0], a[1] + v[1], a[2] + v[2]))
    N_A = len(A_sites)
    N_B = len(B)
    return N_A, N_B, N_B - y * N_A


# neighbour-vector sets V (in the scaled-by-2 integer lattice)
V_CORNERS_8 = [(s1, s2, s3) for s1 in (1, -1) for s2 in (1, -1) for s3 in (1, -1)]
V_AXIS_6 = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
V_TETRA_4 = [(1, 1, 1), (-1, -1, 1), (-1, 1, -1), (1, -1, -1)]  # even # of minuses


# --------------------------------------------------------------------------
# A-site generators for an a x b x c box (cell indices i in [0,a), etc.)
# --------------------------------------------------------------------------

def A_full(a, b, c):
    """Full SC: A = (2i,2j,2k), all i,j,k.  (CsCl)"""
    return [(2 * i, 2 * j, 2 * k) for i in range(a) for j in range(b) for k in range(c)]


def A_even(a, b, c):
    """Checkerboard on the x2 frame: A = (2i,2j,2k) with i+j+k even.
    (fluorite, zincblende -- B sits at body-centre interstitials.)"""
    return [(2 * i, 2 * j, 2 * k)
            for i in range(a) for j in range(b) for k in range(c)
            if (i + j + k) % 2 == 0]


def A_even_unit(a, b, c):
    """Checkerboard on the unit lattice: A = (i,j,k) with i+j+k even.  (NaCl --
    here A and B live on the SAME simple-cubic lattice, different colours, so B
    is reached by axis steps of one lattice spacing, not body-centre offsets.)"""
    return [(i, j, k)
            for i in range(a) for j in range(b) for k in range(c)
            if (i + j + k) % 2 == 0]


# --------------------------------------------------------------------------
# Per-structure brute-force counts
# --------------------------------------------------------------------------

def fluorite_box(a, b, c):
    """Fluorite AB2: A even, V = 8 corners, y = 2.  (reference, re-derived here)"""
    return saturate_count(A_even(a, b, c), V_CORNERS_8, y=2)


def cscl_box(a, b, c):
    """CsCl AB: A full, V = 8 corners, y = 1."""
    return saturate_count(A_full(a, b, c), V_CORNERS_8, y=1)


def nacl_box(a, b, c):
    """NaCl AB: A even (unit lattice), V = 6 axis, y = 1."""
    return saturate_count(A_even_unit(a, b, c), V_AXIS_6, y=1)


def zincblende_box(a, b, c):
    """Zincblende AB: A even, V = 4 tetrahedron, y = 1."""
    return saturate_count(A_even(a, b, c), V_TETRA_4, y=1)


# --------------------------------------------------------------------------
# Closed forms (verified in verify_structures.py)
# --------------------------------------------------------------------------

def _ceil(x):
    return -(-x // 2) * 2 // 2  # not used; use the (n+1)//2 idiom below


def ceil_half(n):
    """ceil(n/2) for integer n >= 0."""
    return (n + 1) // 2


def floor_half(n):
    return n // 2


# --- Fluorite (reference; matches the fluorite-only project) ---
def fluorite_NA(a, b, c):
    return ceil_half(a * b * c)   # = ceil(abc/2)


def fluorite_NB_cube(n):
    """(n+1)^3 - 2(1 + (-1)^n)  (i.e. -4 on even n, 0 on odd n)."""
    return (n + 1) ** 3 - 2 * (1 + ((-1) ** n))


def fluorite_NA_cube(n):
    return ceil_half(n ** 3)


def fluorite_excess_cube(n):
    return 3 * n * (n + 1) - 3 * ((1 + ((-1) ** n)) // 2)


# --- CsCl: N_A = abc, N_B = (a+1)(b+1)(c+1), excess = (a+1)(b+1)(c+1) - abc ---
def cscl_NA(a, b, c):
    return a * b * c


def cscl_NB(a, b, c):
    return (a + 1) * (b + 1) * (c + 1)


def cscl_excess(a, b, c):
    return (a + 1) * (b + 1) * (c + 1) - a * b * c   # = ab+ac+bc + a+b+c + 1


# --- NaCl: N_A = ceil(abc/2);  N_B = floor(abc/2) + sum of 6 face contributions.
#     Each pair of faces normal to dimension d (d=a,b,c) with face-area F:
#       left face (neighbour at coord 0, even)  -> ceil(F/2)
#       right face (neighbour at coord d-1)      -> ceil(F/2) if d odd, else floor(F/2)
#     i.e. face_pair(d,F) = 2*ceil(F/2) if d odd, else F.
#     excess = N_B - ceil(abc/2).  (Cube: 3n^2 [+2 if n odd]; NO edge term.) ---
def _face_pair(dim, area):
    left = ceil_half(area)
    right = ceil_half(area) if (dim % 2 == 1) else floor_half(area)
    return left + right


def nacl_NA(a, b, c):
    return ceil_half(a * b * c)


def nacl_NB(a, b, c):
    return (floor_half(a * b * c)
            + _face_pair(a, b * c) + _face_pair(b, a * c) + _face_pair(c, a * b))


def nacl_excess(a, b, c):
    return nacl_NB(a, b, c) - nacl_NA(a, b, c)


# --- Zincblende: discovered empirically from brute force, then verified.
#     Cube:   excess = 3 n (n+1) / 2          (period 1, no parity term)
#     Box:    excess = (ab + ac + bc + a + b + c + corner) / 2,
#             where corner = 0 if a,b,c share parity, else -3.
#             (Equivalent: subtract 3/2 per "mixed-parity" box.) ---
def _zb_corner(a, b, c):
    return 0 if (a % 2 == b % 2 == c % 2) else -3


def zincblende_NA(a, b, c):
    return ceil_half(a * b * c)


def zincblende_excess(a, b, c):
    return (a * b + a * c + b * c + a + b + c + _zb_corner(a, b, c)) // 2


def zincblende_NB(a, b, c):
    return zincblende_NA(a, b, c) + zincblende_excess(a, b, c)


def zincblende_excess_cube(n):
    return 3 * n * (n + 1) // 2


def zincblende_NB_cube(n):
    return ceil_half(n ** 3) + zincblende_excess_cube(n)


# --------------------------------------------------------------------------
# Asymptotic coefficient c  ( ratio = y + c * N_A^(-1/3) )  for cubes
# --------------------------------------------------------------------------
import math


def asymptotic_c(density_A, surface_per_unit):
    """For a cube of side n with A-density `density_A` (fraction of Z^3 cells)
    and leading excess ~ `surface_per_unit` * n^2, the ratio coefficient is
        c = surface_per_unit / (density_A)^(2/3) * (1)   [since N_A ~ density*n^3]
    Actually ratio-2... = excess/N_A ~ (surf*n^2)/(density*n^3) = (surf/density)/n
    and n ~ (N_A/density)^(1/3), so c = surf/density * density^(1/3) = surf/density^(2/3).
    """
    return surface_per_unit / (density_A ** (2.0 / 3.0))


# Observed / derived leading surface coefficients (excess ~ surf * n^2 for cube):
SURF_FLUORITE = 3.0    # 3n^2 (+3n edge)
SURF_CSCL = 3.0        # 3n^2 (+3n+1)
SURF_NACL = 3.0        # 3n^2 (no edge term)
SURF_ZINCBLENDE = 1.5  # (3/2) n^2  (CN=4: only 4 B stick out per A)
SURF_RUTILE = 8.0      # 8n^2 (-n)  -- tetragonal cube, 2 Ti per cell

C_FLUORITE = asymptotic_c(0.5, SURF_FLUORITE)     # 3 * 2^(2/3) ~ 4.762
C_CSCL = asymptotic_c(1.0, SURF_CSCL)             # 3
C_NACL = asymptotic_c(0.5, SURF_NACL)             # 3 * 2^(2/3) ~ 4.762
C_ZINCBLENDE = asymptotic_c(0.5, SURF_ZINCBLENDE) # 1.5 * 2^(2/3) ~ 2.381
C_RUTILE = asymptotic_c(2.0, SURF_RUTILE)         # 8 / 2^(2/3) = 4*2^(1/3) ~ 5.040


# --------------------------------------------------------------------------
# Rutile (TiO2, AB2) -- the hard multi-site case
# --------------------------------------------------------------------------
# Tetragonal: a=b=1, c=0.644, u=0.305 (Ti at 2a, O at 4f general position).
# Ti sublattices: (0,0,0) and (1/2,1/2,1/2).
# O positions (Wyckoff 4f): (u,u,0), (-u,-u,0), (1/2+u,1/2-u,1/2), (1/2-u,1/2+u,1/2)
# Neighbour vectors V_Ti are determined numerically (6 nearest O per Ti), then
# the saturate count is exact given those vectors.

RUTILE_C = Fraction(161, 250)   # 0.644
RUTILE_U = Fraction(61, 200)    # 0.305
RUTILE_A = Fraction(1, 1)


def _rutile_O_positions():
    """All O fractional coords in the unit cell (4 sites)."""
    u, c = RUTILE_U, RUTILE_C
    half = Fraction(1, 2)
    return [(u, u, 0), (-u, -u, 0),
            (half + u, half - u, half),
            (half - u, half + u, half)]


def _rutile_Ti_positions():
    half = Fraction(1, 2)
    return [(0, 0, 0), (half, half, half)]


def _nearest_O_vectors(ti_frac, nvoronoi=6, cell_range=2):
    """Return the `nvoronoi` nearest O-vectors (as Fraction 3-tuples) to a Ti
    at fractional coord `ti_frac`, scanning O in cells [−cell_range,cell_range]."""
    O_base = _rutile_O_positions()
    a, c = RUTILE_A, RUTILE_C
    cands = []
    for (di, dj, dk) in product(range(-cell_range, cell_range + 1), repeat=3):
        for (ox, oy, oz) in O_base:
            fx = ox + di
            fy = oy + dj
            fz = oz + dk
            # physical displacement from ti_frac
            dx = (fx - ti_frac[0]) * a
            dy = (fy - ti_frac[1]) * a
            dz = (fz - ti_frac[2]) * c
            dist2 = dx * dx + dy * dy + dz * dz
            vec = (dx, dy, dz)
            cands.append((dist2, vec))
    cands.sort(key=lambda t: t[0])
    return [v for _, v in cands[:nvoronoi]]


def rutile_V():
    """Neighbour vectors for the two Ti sublattices (6 O each), as Fraction tuples."""
    Ti = _rutile_Ti_positions()
    return [_nearest_O_vectors(ti) for ti in Ti]


def rutile_box(na, nb, nc):
    """Tetragonal na x nb x nc box of cells; saturate Ti.  Returns (N_Ti, N_O, excess).
    excess = N_O - 2*N_Ti  (rutile is AB2).  Brute force, exact given V."""
    V = rutile_V()  # list of two 6-vector lists
    Ti_base = _rutile_Ti_positions()
    a, c = RUTILE_A, RUTILE_C
    A_sites = []
    for i in range(na):
        for j in range(nb):
            for k in range(nc):
                for (tx, ty, tz) in Ti_base:
                    px = (i + tx) * a
                    py = (j + ty) * a
                    pz = (k + tz) * c
                    A_sites.append((px, py, pz))
    # attach the correct V per Ti sublattice (alternates with cell parity of Ti_base)
    B = set()
    idx = 0
    for i in range(na):
        for j in range(nb):
            for k in range(nc):
                for s, (tx, ty, tz) in enumerate(Ti_base):
                    px = (i + tx) * a
                    py = (j + ty) * a
                    pz = (k + tz) * c
                    for v in V[s]:
                        B.add((px + v[0], py + v[1], pz + v[2]))
                    idx += 1
    N_A = len(A_sites)
    N_B = len(B)
    return N_A, N_B, N_B - 2 * N_A


# Rutile cube closed form (empirical, verified to large n in verify_structures.py):
#   N_Ti = 2 n^3,   excess = 8 n^2 - n,   N_O = 4 n^3 + 8 n^2 - n.
# Despite O sitting on general positions (parameter u), the count is clean
# because for generic u the O sites never coincide, so the count depends only
# on the (fixed) Ti-O neighbour topology, not on u itself.
def rutile_NTl_cube(n):
    return 2 * n ** 3


def rutile_excess_cube(n):
    return 8 * n * n - n


def rutile_NO_cube(n):
    return 4 * n ** 3 + 8 * n * n - n


# --------------------------------------------------------------------------
# Cube-shell tables
# --------------------------------------------------------------------------

def cube_table(struct_box, n_max):
    return [((n,) + struct_box(n, n, n)) for n in range(1, n_max + 1)]


if __name__ == "__main__":
    print("Neighbour vectors V for each cubic structure (scaled-by-2 lattice):")
    print("  fluorite  V8 =", V_CORNERS_8)
    print("  CsCl      V8 =", V_CORNERS_8, " (same V, but A full not checkerboard)")
    print("  NaCl      V6 =", V_AXIS_6)
    print("  zincblende V4=", V_TETRA_4)
    print()
    print("n x n x n saturated cubes  (N_A, N_B, excess):")
    print(f"{'n':>3} {'fluorite':>18} {'CsCl':>14} {'NaCl':>14} {'zincblende':>16}")
    for n in range(1, 9):
        f = fluorite_box(n, n, n)
        cs = cscl_box(n, n, n)
        na = nacl_box(n, n, n)
        zb = zincblende_box(n, n, n)
        print(f"{n:>3}  {f[0]:>4},{f[1]:>5},{f[2]:>4}   "
              f"{cs[0]:>3},{cs[1]:>4},{cs[2]:>3}   "
              f"{na[0]:>3},{na[1]:>4},{na[2]:>3}   "
              f"{zb[0]:>3},{zb[1]:>4},{zb[2]:>3}")
    print()
    print("rutile (TiO2) tetragonal n x n x n box (N_Ti, N_O, excess = N_O-2N_Ti):")
    for n in range(1, 7):
        rt = rutile_box(n, n, n)
        print(f"  n={n}: N_Ti={rt[0]:>5}  N_O={rt[1]:>6}  excess={rt[2]:>5}  ratio={rt[1]/rt[0]:.4f}")
