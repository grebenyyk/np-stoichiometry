"""
oeis_check_structures.py -- live OEIS lookup for the cube sequences of the
new structures (CsCl, NaCl, zincblende, rutile), to see which are already
catalogued and which are novel.

Run:  python3 oeis_check_structures.py   (needs internet)
"""
import sys, os, json, urllib.request

sys.path.insert(0, os.path.dirname(__file__))
import structures as S

OEIS = "https://oeis.org/search?fmt=json&q="

def oeis_search(q):
    try:
        url = OEIS + urllib.parse.quote(q)
        req = urllib.request.Request(url, headers={"User-Agent": "axby-np-stoichiometry/1.0 (research)"})
        with urllib.request.urlopen(req, timeout=25) as r:
            d = json.loads(r.read().decode())
    except Exception as e:
        print(f"    (network error: {e})")
        return None
    return d if d is not None else []

def seq(fn, idx):
    return ",".join(str(fn(n, n, n)[idx]) for n in range(1, 13))

queries = {
    "CsCl   N_B":    seq(S.cscl_box, 1),        # (n+1)^3
    "CsCl   excess": seq(S.cscl_box, 2),        # 3n^2+3n+1
    "NaCl   N_B":    seq(S.nacl_box, 1),
    "NaCl   excess": seq(S.nacl_box, 2),
    "zblende N_B":   seq(S.zincblende_box, 1),
    "zblende excess":seq(S.zincblende_box, 2),  # 3n(n+1)/2
    "rutile N_O":    seq(S.rutile_box, 1),
    "rutile excess": seq(S.rutile_box, 2),      # 8n^2 - n
}

print("Live OEIS lookup for new-structure cube sequences\n" + "=" * 72)
for name, q in queries.items():
    res = oeis_search(q)
    if res is None:
        verdict = "NETWORK ERROR"
    elif len(res) == 0:
        verdict = "NO MATCH  -> novel candidate"
    else:
        r = res[0]
        verdict = f"A{r['number']:06d}  {r['name'][:55]}"
    print(f"  {name:16s} {verdict}")
    print(f"    query: {q}")
