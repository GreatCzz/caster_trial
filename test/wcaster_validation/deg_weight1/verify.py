#!/usr/bin/env python3
"""Test1 verify: alignment_wcaster (weights->1) should match CASTER.

Checks:
  1. per-species weights ~= 1
  2. per-line Score: values close (search order identical: same seed, no taxon prioritizing)
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

cs = parse_scores(RES / 'caster.log')
ws = parse_scores(RES / 'alignment_wcaster.log')
wts = parse_weights(RES / 'alignment_wcaster.log')

# 1. weights ~= 1
print('== weights (should be ~=1) ==')
for k, v in wts.items():
    print(f'  {k}: {v}')
w_ok = all(abs(v - 1) < 0.01 for k, v in wts.items() if k != 'ref')

# 2. per-line score closeness (same search trajectory)
print(f'\n== scores == caster: {len(cs)} lines, wcaster: {len(ws)} lines')
n = min(len(cs), len(ws))
close = []
for i in range(n):
    if abs(cs[i]) < 1e-9:
        rel = abs(cs[i] - ws[i])
    else:
        rel = abs(cs[i] - ws[i]) / abs(cs[i])
    close.append(rel)
# ignore first trivial (tree-building) steps maybe; report max rel diff over all
max_rel = max(close) if close else 0
print(f'  max |relative diff| over {n} aligned Score lines: {max_rel:.6g}')
# Compare also scale-adjusted: score scales ~ product of weights; expect within a few %.
score_ok = max_rel < 0.05

# 3. RF final tree
rf = rf_distance((RES / 'caster.tre').read_text().strip(),
                 (RES / 'alignment_wcaster.tre').read_text().strip())
print(f'\nRF distance: {rf}')

ok = w_ok and score_ok and rf == 0
print(f'\nVERDICT: {"PASS" if ok else "FAIL"}'
      f' (weights~1:{w_ok} scores-close:{score_ok} RF=0:{rf==0})')
sys.exit(0 if ok else 1)
