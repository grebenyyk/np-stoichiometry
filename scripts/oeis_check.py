"""
oeis_check.py -- query the OEIS JSON API live to confirm which of the
cubic-shell sequences are already catalogued and which are novel.

Run:  python3 oeis_check.py
(requires internet access to https://oeis.org)
"""
import sys, os, json, urllib.request

sys.path.insert(0, os.path.dirname(__file__))
from fluorite import cubic_shell_table, OEIS_HITS, OEIS_NOVEL

OEIS = "https://oeis.org/search?fmt=json&q="

def oeis_search(q):
    """Return list of result dicts, or [] if no match (OEIS returns null)."""
    try:
        url = OEIS + urllib.parse.quote(q)
        req = urllib.request.Request(url, headers={"User-Agent": "np-stoichiometry/1.0 (research)"})
        with urllib.request.urlopen(req, timeout=20) as r:
            d = json.loads(r.read().decode())
    except Exception as e:
        print(f"    (network error: {e})")
        return None
    return d if d is not None else []

rows = cubic_shell_table(15)
seqs = {
    "N_A full (n=1..)":   ",".join(str(r[1]) for r in rows),
    "N_B full (n=1..)":   ",".join(str(r[2]) for r in rows),
    "excess full (n=1..)":",".join(str(r[3]) for r in rows),
    "N_A odd shells":     ",".join(str(rows[i][1]) for i in range(0, 15, 2)),
    "N_B odd shells":     ",".join(str(rows[i][2]) for i in range(0, 15, 2)),
    "N_B even shells":    ",".join(str(rows[i][2]) for i in range(1, 15, 2)),
    "excess odd shells":  ",".join(str(rows[i][3]) for i in range(0, 15, 2)),
    "excess even shells": ",".join(str(rows[i][3]) for i in range(1, 15, 2)),
}

print("Live OEIS lookup (https://oeis.org)\n" + "=" * 72)
for name, q in seqs.items():
    res = oeis_search(q)
    if res is None:
        verdict = "NETWORK ERROR"
    elif len(res) == 0:
        verdict = "NO MATCH  -> novel candidate"
    else:
        r = res[0]
        verdict = f"A{r['number']:06d}  {r['name'][:60]}"
    print(f"  {name:20s} {verdict}")
    print(f"    query: {q}")

print("\n" + "=" * 72)
print("Hard-coded identifications (from fluorite.py):")
print("  HIT  ->", OEIS_HITS)
print("  NEW  ->", OEIS_NOVEL)
