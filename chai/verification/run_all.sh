#!/bin/sh
# Runs the checks from scratch in this directory.
#   sh run_all.sh          about 15 seconds; uses the shipped results.json / results_consts_*.json (from --full runs) for the supplement
#   sh run_all.sh --full   adds the five exhaustive binary closures (about 5 minutes) and rewrites results.json
# Any failing script stops the run with a non-zero exit status (set -e).
set -e
cd "$(dirname "$0")"
python3 -I chai_checker.py "$@"           # default: results_default.json/.md; with --full: results.json/.md; also unary_catalogue.json
python3 -I classification_table.py         # classification_127.csv (two methods, asserted equal)
python3 -I independent_free_constants.py   # independent count 40 / 95 and per-subset comparison with classification_127.csv
python3 -I check_c9_chain.py               # polarization, difference lift, G3, three-slice C9
python3 -I independent_check.py            # independent re-implementation of terms and unary sizes
python3 -I verify_direction_rotation.py    # v10r3: quarter turn, folding, congruences, sign rule, constants case (Prop. 6), XY3/G3, orbits (stdlib only)
python3 -I chai_verify_v1_increments.py ../paper/Chai2026v1_EN.tex   # Chai2026v1: derivation, printed tables, constructions on every binary function (stdlib only)
python3 -I make_tables.py                  # unary_catalogue.csv, S5_unary_table.tex
python3 -I make_supplement.py              # chai2026v1_supplement.tex (asserts that all result files are present and passed)
# pdflatex chai2026v1_supplement.tex (twice) builds the PDF; it is shipped in ../paper/
sh make_sha.sh                             # SHA256SUMS for every file in this directory
