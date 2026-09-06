#!/usr/bin/env python3
"""compare_topology.py — cat chunk test: RF between chunk_wtrial, chunk_wcaster, and true tree."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / 'lib'))
from tree_utils import rf_distance

HERE = Path(__file__).resolve().parent
PROJ = HERE.parents[2]
RES = HERE / 'results'

wt = (RES / 'chunk_wtrial.tre').read_text().strip()
wc = (RES / 'chunk_wcaster.tre').read_text().strip()
truth = (PROJ / 'example/cat/true_species_tree.nwk').read_text().strip()

print(f'chunk_wtrial:   {wt}')
print(f'chunk_wcaster:  {wc}')
print(f'true species tree: {truth}')

d_wt_wc = rf_distance(wt, wc)
d_wt_t = rf_distance(wt, truth)
d_wc_t = rf_distance(wc, truth)

print(f'\nRF(chunk_wtrial, chunk_wcaster) = {d_wt_wc}')
print(f'RF(chunk_wtrial, true)           = {d_wt_t}')
print(f'RF(chunk_wcaster, true)          = {d_wc_t}')

ok = d_wt_wc == 0 and d_wt_t == 0 and d_wc_t == 0
print(f'\nVERDICT: {"PASS" if ok else "FAIL"}  (all RF==0)')
sys.exit(0 if ok else 1)
