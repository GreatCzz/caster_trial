#!/usr/bin/env python3
"""compare_topology.py — RF distance between chunk_wtrial and chunk_wcaster output trees."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / 'lib'))
from tree_utils import rf_distance

HERE = Path(__file__).resolve().parent
RES = HERE / 'results'

wt = (RES / 'chunk_wtrial.tre').read_text().strip()
wc = (RES / 'chunk_wcaster.tre').read_text().strip()
rf = rf_distance(wt, wc)

print(f'chunk_wtrial:   {wt}')
print(f'chunk_wcaster:  {wc}')
print(f'\nRF distance: {rf}')
ok = rf == 0
print(f'VERDICT: {"PASS" if ok else "FAIL"}  (RF==0 => identical topology)')
sys.exit(0 if ok else 1)
