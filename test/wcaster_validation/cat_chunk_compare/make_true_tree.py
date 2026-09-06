#!/usr/bin/env python3
"""make_true_tree.py — Generate example/cat/true_species_tree.nwk once.

Reuses the existing CASTER cat output (results_20260722_1/caster.tre), strips
branch lengths / bootstrap support, keeps real species names. Future cat tests
reference this file instead of re-running CASTER.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / 'lib'))
from tree_utils import clean_newick

HERE = Path(__file__).resolve().parent
PROJ = HERE.parents[2]                          # .../project/caster_trial
SRC = PROJ / 'example/cat/results_20260722_1/caster.tre'
OUT = PROJ / 'example/cat/true_species_tree.nwk'

nw = SRC.read_text().strip()
clean = clean_newick(nw) + ';\n'
OUT.write_text(clean)
print(f'wrote {OUT}:')
print(clean)
