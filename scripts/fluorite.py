"""
fluorite.py -- combinatorial model of the stoichiometry of a saturated
fluorite (AB2) nanoparticle.

Recast of the fluorite structure
--------------------------------
In the conventional cubic cell (side a) of fluorite CaF2:
  * A (cation) occupies the FCC positions;
  * B (anion) occupies all 8 tetrahedral holes.

Equivalently (this is the key simplification used throughout):

  * B forms a SIMPLE-CUBIC lattice with spacing a/2
        B = (a/2) Z^3 + (a/4)(1,1,1).
  * A sits at the body centre of every SECOND B-cube, in a 3D checkerboard:
        A = centre of cube [i,i+1]x[j,j+1]x[k,k+1]   with  i+j+k even.
  * Hence  CN(A) = 8  (A is cube-cornered by 8 B),
          CN(B) = 4  (each B is shared by 4 A-cubes).

We work in B-lattice units (B spacing = 1, B-cube volume = 1), so the
number densities are  rho_A = 1/2 ,  rho_B = 1 .

The saturation model ("no surface artifacts")
---------------------------------------------
A particle is specified by a finite set S of A-sites (the cation cluster,
i.e. the shape region).  The B-set is then

        B(S) = union over A in S of the 8 corner-B's of A's cube,

i.e. EVERY A -- including surface ones -- keeps its full 8-fold B shell.
Nothing is added, removed, or reconstructed.  This is exactly the idealised
"surface A atoms saturated with B, no surface artifacts" limit.

All exact box and cube counts below are verified against brute-force
enumeration in verify_sequences.py and verify_box_formulas.py.
"""

from itertools import product
import math

# --------------------------------------------------------------------------
# Brute-force enumeration (ground truth)
# --------------------------------------------------------------------------

def black_cubes_in_box(a, b, c):
    """A-sites = black cubes (i+j+k even) inside [0,a)x[0,b)x[0,c)."""
    return [(i, j, k)
            for i in range(a) for j in range(b) for k in range(c)
            if (i + j + k) % 2 == 0]


def cube_corners(cube):
    """The 8 B-lattice corners of a unit B-cube indexed by `cube`."""
    i, j, k = cube
    return [(i + dx, j + dy, k + dz)
            for dx in (0, 1) for dy in (0, 1) for dz in (0, 1)]


def cluster_counts(A_sites):
    """Return (N_A, N_B, excess = N_B - 2 N_A) for a given A-set."""
    B = set()
    for cube in A_sites:
        B.update(cube_corners(cube))
    N_A = len(A_sites)
    N_B = len(B)
    return N_A, N_B, N_B - 2 * N_A


def cube_cluster(n):
    """n x n x n cube, B-axis aligned (brute force)."""
    return cluster_counts(black_cubes_in_box(n, n, n))


def box_cluster(a, b, c):
    """a x b x c box, B-axis aligned (brute force)."""
    return cluster_counts(black_cubes_in_box(a, b, c))


def sphere_cluster(R):
    """A-sites whose cube centre lies within radius R of the origin
    (centred at (1/2,1/2,1/2) so the origin is a cube corner).  Brute force."""
    m = int(math.ceil(R)) + 1
    R2 = R * R
    sites = []
    for i in range(-m, m + 1):
        for j in range(-m, m + 1):
            for k in range(-m, m + 1):
                if (i + j + k) % 2 != 0:
                    continue
                if (i + 0.5) ** 2 + (j + 0.5) ** 2 + (k + 0.5) ** 2 <= R2:
                    sites.append((i, j, k))
    return cluster_counts(sites)


# --------------------------------------------------------------------------
# Closed forms
# --------------------------------------------------------------------------

def NA_box_cf(a, b, c):
    """N_A = ceil(a b c / 2).  (Number of even-sum triples in an a x b x c box.)"""
    return (a * b * c + (1 if (a % 2 and b % 2 and c % 2) else 0)) // 2


def _epsilon(a, b, c):
    """epsilon = 1 if a,b,c all odd (so abc odd), else 0."""
    return 1 if (a % 2 and b % 2 and c % 2) else 0


def U_box(a, b, c):
    """Number of the 8 box *vertices* that are NOT a corner of any in-box
    A-cube (a parity-dependent corner correction, in {0,1,2,3,4}).

    A box vertex (x,y,z), x in {0,a} etc., is incident to a unique in-box
    cube (min(x,a-1), min(y,b-1), min(z,c-1)); it is "lost" iff that cube is
    WHITE (index sum odd).  Enumerating the 8 vertices gives the seven
    parity conditions below (the (0,0,0) vertex is never lost).
    """
    def iv(p):
        return 1 if p else 0
    return (iv(a % 2 == 0) + iv(b % 2 == 0) + iv(c % 2 == 0)
            + iv((a + b) % 2 == 1) + iv((a + c) % 2 == 1) + iv((b + c) % 2 == 1)
            + iv((a + b + c) % 2 == 0))


def NB_box_cf(a, b, c):
    """N_B = (a+1)(b+1)(c+1) - U(a,b,c)."""
    return (a + 1) * (b + 1) * (c + 1) - U_box(a, b, c)


def excess_box_cf(a, b, c):
    """Exact excess for an a x b x c box, decomposed as

        excess = (ab+ac+bc)            [surface = sum of 3 projections]
               + (a + b + c)           [edge term  = 1/4 of total edge length]
               + (1 - U - eps)         [corner / lattice-parity term]

    where eps = [a,b,c all odd].
    """
    eps = _epsilon(a, b, c)
    surface = a * b + a * c + b * c
    edge = a + b + c
    corner = 1 - U_box(a, b, c) - eps
    return surface + edge + corner


def NA_cube_cf(n):
    return NA_box_cf(n, n, n)


def NB_cube_cf(n):
    """N_B(n) = (n+1)^3 - 2(1 + (-1)^n)   (i.e. -4 on even n, 0 on odd n)."""
    return (n + 1) ** 3 - 2 * (1 + ((-1) ** n))


def excess_cube_cf(n):
    """excess(n) = 3 n (n+1) - (3/2)(1 + (-1)^n)
       = 3 n (n+1)        if n odd
       = 3 (n^2 + n - 1)  if n even
    """
    return 3 * n * (n + 1) - 3 * ((1 + ((-1) ** n)) // 2)


# --------------------------------------------------------------------------
# Axis-aligned asymptotic helpers
# --------------------------------------------------------------------------


def asymptotic_ratio(shape_c, N_A):
    """N_B / N_A ~ 2 + c * N_A^(-1/3) for fixed axis-aligned box shapes."""
    return 2.0 + shape_c * N_A ** (-1.0 / 3.0)


def box_asymptotic_c(alpha, beta, gamma):
    """Leading ratio coefficient for boxes m*alpha by m*beta by m*gamma."""
    surface = alpha * beta + alpha * gamma + beta * gamma
    volume_density = alpha * beta * gamma / 2.0
    return surface / (volume_density ** (2.0 / 3.0))


C_CUBE = box_asymptotic_c(1, 1, 1)             # ~ 4.7622


# --------------------------------------------------------------------------
# OEIS sequence material (cubic shells, n = 1, 2, 3, ...)
# --------------------------------------------------------------------------

def cubic_shell_table(n_max):
    """Rows: (n, N_A, N_B, excess) for the n x n x n saturated cube."""
    rows = []
    for n in range(1, n_max + 1):
        na, nb, e = cube_cluster(n)
        rows.append((n, na, nb, e))
    return rows


# OEIS identifications (rechecked live by oeis_check.py when the network works):
OEIS_HITS = {
    "N_A_full":       ("A036486", "ceiling(n^3/2)",                 "1,4,14,32,63,108,172,256,365,500,..."),
    "N_A_odd_shells": ("A050492", "thickened cube numbers 4m^3-6m^2+3m", "1,14,63,172,365,666,1099,..."),
    "N_B_odd_shells": ("A016743", "even cubes (2m)^3",              "8,64,216,512,1000,1728,..."),
    "excess_odd":     ("A152746", "6 * hexagonal numbers 6m(2m-1)", "6,36,90,168,270,396,..."),
}

OEIS_NOVEL = {
    "N_B_full":       ("(n+1)^3 - 2(1+(-1)^n)",        "8,23,64,121,216,339,512,725,1000,1327,..."),
    "N_B_even_shells":("(2m+1)^3 - 4",                 "23,121,339,725,1327,2193,..."),
    "excess_full":    ("3n(n+1) - (3/2)(1+(-1)^n)",    "6,15,36,57,90,123,168,213,270,327,..."),
    "excess_even":    ("3(n^2+n-1) at n=2m",           "15,57,123,213,327,..."),
}


if __name__ == "__main__":
    print("n x n x n saturated fluorite cube  (brute force == closed form):")
    print(f"{'n':>3} {'N_A':>6} {'N_B':>7} {'excess':>7}   "
          f"{'NA_cf':>6} {'NB_cf':>7} {'E_cf':>7}")
    for n in range(1, 11):
        na, nb, e = cube_cluster(n)
        print(f"{n:>3} {na:>6} {nb:>7} {e:>7}   "
              f"{NA_cube_cf(n):>6} {NB_cube_cf(n):>7} {excess_cube_cf(n):>7}")
