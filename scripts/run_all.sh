#!/usr/bin/env bash
# run_all.sh -- run every verification, generate figures + data table.
# Usage:  bash run_all.sh
set -u
cd "$(dirname "$0")"
PY=python3

echo "################################################################"
echo "# 1. Sequence closed-forms vs brute force + OEIS data + bonds  #"
echo "################################################################"
$PY verify_sequences.py || echo "(verify_sequences failed)"
echo
echo "################################################################"
echo "# 2. Exact axis-aligned box formulas                            #"
echo "################################################################"
$PY verify_box_formulas.py || echo "(verify_box_formulas failed)"
echo
echo "################################################################"
echo "# 3. Live OEIS lookup (needs internet)                         #"
echo "################################################################"
$PY oeis_check.py || echo "(oeis_check skipped/failed -- offline?)"
echo
echo "################################################################"
echo "# 4. Figures for the writeup                                    #"
echo "################################################################"
$PY plot_ratio.py || echo "(plot_ratio failed -- matplotlib?)"
echo
echo "################################################################"
echo "# 5. Data table  ->  ../data/sequences.csv                      #"
echo "################################################################"
$PY - <<'PY'
import sys, os, csv
sys.path.insert(0, ".")
from fluorite import cubic_shell_table, NA_cube_cf, NB_cube_cf, excess_cube_cf
rows = cubic_shell_table(30)
out = os.path.join("..", "data", "sequences.csv")
with open(out, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["n","N_A","N_B","excess","N_A_closed","N_B_closed","excess_closed",
                "N_B_minus_2N_A","ratio_NBoverNA"])
    for n,na,nb,e in rows:
        w.writerow([n,na,nb,e,NA_cube_cf(n),NB_cube_cf(n),excess_cube_cf(n),
                    e, round(nb/na,6)])
print("wrote", out, "with", len(rows), "rows")
PY
echo
echo "done."
