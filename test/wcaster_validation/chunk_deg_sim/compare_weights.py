#!/usr/bin/env python3
"""compare_weights.py — Compare per-(chunk,species) weights between chunk_wtrial and chunk_wcaster.

Both logs contain blocks:
  Chunk weights for alignment file: <path>
  chunk <i> <species> <weight>
  ...
Parse into {(chunk, species): weight}; pass iff identical species/chunk sets and |Δ|<1e-12.
"""
import re, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RES = HERE / 'results'

def parse(log):
    w = {}
    in_block = False
    for line in open(log):
        s = re.sub(r'^\s*-\s*', '', line).strip()  # strip LogInfo " - " prefix
        if s.startswith('Chunk weights for alignment file:'):
            in_block = True
            continue
        if not in_block:
            continue
        m = re.match(r'^chunk (\d+) (\S+) ([-0-9.eE+]+)$', s)
        if m:
            w[(int(m.group(1)), m.group(2))] = float(m.group(3))
        else:
            in_block = False  # block ended (next non-chunk line)
    return w

wt = parse(RES / 'chunk_wtrial.log')
wc = parse(RES / 'chunk_wcaster.log')

print(f'wtrial entries: {len(wt)}, wcaster entries: {len(wc)}')
if not wt or not wc:
    print('ERROR: no chunk-weight block found (did you pass --dump-chunk-weights?)')
    sys.exit(1)

common = sorted(set(wt) & set(wc))
only_wt = sorted(set(wt) - set(wc))
only_wc = sorted(set(wc) - set(wt))

print(f'common (chunk,species): {len(common)}, only-wtrial: {len(only_wt)}, only-wcaster: {len(only_wc)}')
for k in only_wt[:10]:
    print(f'  only wtrial: chunk={k[0]} {k[1]}')
for k in only_wc[:10]:
    print(f'  only wcaster: chunk={k[0]} {k[1]}')

maxdiff = 0.0
for k in common:
    maxdiff = max(maxdiff, abs(wt[k] - wc[k]))
print(f'max |diff| over common: {maxdiff:.3g}')

# sample of the first few rows for eyeballing
print('\nsample (first 6 common entries):')
for k in common[:6]:
    print(f'  chunk {k[0]:>4} {k[1]:<4} wtrial={wt[k]:.10f} wcaster={wc[k]:.10f}')

ok = not only_wt and not only_wc and maxdiff < 1e-12
print(f'\nVERDICT: {"PASS" if ok else "FAIL"}  (sets equal: {not only_wt and not only_wc}, max|diff|={maxdiff:.3g})')
sys.exit(0 if ok else 1)
