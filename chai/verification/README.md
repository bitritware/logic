# Chai verification package

Code and data that accompany the paper *Chai: A Signed Three-Valued Algebra from Parameterized Boolean Negation*.
Every number in the paper that comes from a finite computation is produced here.

## What is checked

On D = {-1, 0, +1} the checker verifies, by exhaustive computation or by explicit terms:

- the parity lemma, the polarization identity AND = DIF(OR, XOR), and the involution (SYN, DIF);
- the explicit terms of the paper (indicator terms, the terms k_c, MAG, the connectives of Table III);
- the five binary operations generate every zero-fixing function (all 6561 binary ones are rebuilt, plus random ones of arity 3 and 4);
- DIF alone generates exactly the 81 odd binary functions;
- the classification of all 127 nonempty subsets of the seven operations: 40 are complete without nullary constants, 95 with them. Two methods are used (certificates, and the Slupecki criterion with unary closure);
- a regression test: the counting procedure gives other numbers when the signature is changed;
- the difference lift, the 2862 orbits, the nine-trit words, and the shortest terms of the 27 unary functions;
- (v10r3) the quarter turn of the truth tables (`verify_direction_rotation.py`, standard library only): the cycle SYN -> DIF -> -SYN -> -DIF, OR <-> XOR, AND -> -AND; the lemma R(T_F)=T_{F o rho} for all 19683 binary F; the complex table Z=T_SYN+i T_DIF (exact Gaussian integers); the folding MAG (compatibility of its kernel with each operation, congruences, subalgebras, automorphisms); Boolean restriction and set identities; C4 and dihedral orbit counts; the 17 sign-linear functions; the sign rule sgn(ab)=sgn(|a+b|-|a-b|); the constants case (95 of 127: the terms and invariants of Proposition 6, an independent Slupecki certificate for every subset, agreement with `classification_127.csv`); the XY3/G3 tables of S17; and negative controls (wrong sign, wrong direction, reflection taken for rotation) that must come out false;
- the chain polarization, difference lift, G3 branches and the three-slice C9 construction (`check_c9_chain.py`, a few seconds): the printed truth tables, all 262144 mask pairs for n=2, and the term expansion of all 19683 binary functions.

In the output files the key `free_constants` and the closures named `with free constants` mean the case "0, +1, -1 are available as nullary operations".

## Requirements

Python 3.10 or newer and numpy 1.26 or newer (`pip install -r requirements.txt`). Results in this package were produced with Python 3.13 and numpy 2.5; the versions are recorded in `results.json`.

## Running

```
python3 check_c9_chain.py           # a few seconds; the chain polarization -> difference lift -> G3 -> C9
python3 verify_direction_rotation.py  # a few seconds; writes direction_rotation_results.json/.log (exit status 1 on any failure)
python3 chai_checker.py            # about 15 seconds; writes results_default.json/.md and unary_catalogue.json
python3 chai_checker.py --full     # adds the five exhaustive binary closures (about 5 minutes); writes results.json/.md
python3 classification_table.py     # the 127 subsets by two methods; writes classification_127.csv
python3 chai_checker.py --full-consts --only=SYN+OR   # closure with nullary constants for one generator set (a few minutes each)
sh run_all.sh                       # default checks, classification, independent scripts, tables, supplement source (add --full for exhaustive closures)
```

The `--only` names used here are DIF, SYN+OR, SYN+XOR, SYN+AND, SYN+ANT and SYN+NOT.
A default run never overwrites `results.json`; only `--full` does. The supplement generator reads `results.json` and the six `results_consts_*.json` files and stops if a closure entry is missing.

A run succeeded if the last line of the output is `ALL ASSERTIONS PASSED`. Every item is an assertion and the script stops at the first failure.

## Files

| file | content |
|---|---|
| `chai_checker.py` | the checker (single file) |
| `independent_check.py` | independent re-implementation of the terms and of the sizes of the 27 unary functions; prints `catalogue mismatches ... 0` |
| `independent_free_constants.py` | independent count of the complete subsets: 40 without, 95 with nullary constants; compares every subset with `classification_127.csv` and prints `mismatches 0` |
| `verify_direction_rotation.py` | v10: independent check of the quarter turn, complex table, folding, congruence, Boolean/set, orbit, sign-rule, constants (Proposition 6) and XY3/G3 statements (S14-S17 of the supplement); 115 tests including negative controls |
| `direction_rotation_results.json`, `.log` | its machine-readable results and readable log (Python version, command, passed/failed counts) |
| `check_c9_chain.py` | audit of polarization, difference lift, G3 and the three-slice C9 construction; ends with `ALL ASSERTIONS PASSED` |
| `make_tables.py` | builds `unary_catalogue.csv` and `S5_unary_table.tex` from `unary_catalogue.json` |
| `classification_table.py` | classification of all 127 subsets by two methods (certificates of incompleteness; Slupecki criterion with unary closure), asserted equal; writes `classification_127.csv` |
| `classification_127.csv` | per subset: complete without / with nullary constants, certificate or reason if incomplete |
| `make_supplement.py` | builds `chai2026v1_supplement.tex` (sections S0-S17) from the result files; the compiled PDF is in `../paper/` |
| `results.json`, `results.md` | results of the `--full` run (used by the supplement) |
| `results_default.json`, `.md` | results of the default run (written when the default run is repeated) |
| `results_consts_*.json`, `.md` | results of the closures with nullary constants, one file per generator set |
| `unary_catalogue.json`, `.csv` | the 27 unary functions: code (prefix 0t), value vector, size, lexicographically smallest minimal term over the seven operations and over {DIF, NOT} |
| `full_run.log`, `consts_*.log` | console output of the long runs |
| `SHA256SUMS` | hashes of all files in this directory except `SHA256SUMS` itself |

Check the hashes with `sha256sum -c SHA256SUMS`. The `results*.json` files contain the Python and numpy versions, so a run on another machine will have different environment entries; compare the results, not whole files.

## Conventions

- Tree size of a term: number of operator nodes, repeated subterms counted each time.
- Ties between minimal terms are broken by the lexicographically smallest printed term, so the output is stable.
- The prefix `0t` marks a ternary literal, as `0b` and `0x` do for binary and hexadecimal. Digits 0, 1, 2 stand for 0, +1, -1.

## Rights

Provisional license: Apache-2.0 (see `../LICENSE`), to be confirmed by the author's patent agent before the first public release. The paper's relation to the author's published patent documents is stated in the paper itself.
