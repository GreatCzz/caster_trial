#!/usr/bin/env python3
"""Test2 verify: alignment_wcaster vs alignment_wtrial on 4 taxa.

With only 4 taxa the sole quartet contains ref, so global (wcaster) and ref-restricted
(wtrial) scoring must agree on every evaluated topology. We compare:
  1. per-species weights (identical)
  2. every logged Score line (identical search scores under same data)
  3. final tree RF == 0
"""
import re, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / 'lib'))
from tree_utils import rf_distance

HERE = Path(__file__).resolve().parent
RES = HERE / 'results'

def parse_scores(log):
    return [float(m) for m in re.findall(r'^\s*-?\s*Score:\s+([0-9.eE+\-]+)', open(log).read(), re.M)]

def parse_weights(log):
    w = {}
    for line in open(log):
        m = re.search(r'^\s*-\s*(\S+)\s*=\s*([0-9.eE+\-]+)', line)
        if m and m.group(1) != 'weights':
            w[m.group(1)] = float(m.group(2))
    return w

wc_scores = parse_scores(RES / 'alignment_wcaster.log')
wt_scores = parse_scores(RES / 'alignment_wtrial.log')
wc_w = parse_weights(RES / 'alignment_wcaster.log')
wt_w = parse_weights(RES / 'alignment_wtrial.log')

# 1. weights identical
print('== weights ==')
print(f'  wcaster: {wc_w}')
print(f'  wtrial:  {wt_w}')
w_ok = all(abs(wc_w.get(k, -1) - v) < 1e-12 for k, v in wt_w.items()) and set(wc_w) == set(wt_w)

# 2. all Score lines equal
print(f'\n== scores == wcaster: {len(wc_scores)} lines, wtrial: {len(wt_scores)} lines')
n = min(len(wc_scores), len(wt_scores))
maxdiff = max((abs(a - b) for a, b in zip(wc_scores, wt_scores)), default=0)
print(f'  max |diff| over {n} aligned lines: {maxdiff:.6g}')
score_ok = len(wc_scores) == len(wt_scores) and all(abs(a - b) < 1e-9 for a, b in zip(wc_scores, wt_scores))

# 3. RF final tree
rf = rf_distance((RES / 'alignment_wcaster.tre').read_text().strip(),
                 (RES / 'alignment_wtrial.tre').read_text().strip())
print(f'\nRF distance: {rf}')

ok = w_ok and score_ok and rf == 0
print(f'\nVERDICT: {"PASS" if ok else "FAIL"}'
      f' (weights-identical:{w_ok} scores-equal:{score_ok} RF=0:{rf==0})')
sys.exit(0 if ok else 1)
