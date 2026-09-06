#!/usr/bin/env python3
"""compare_weights.py — Compare per-species weights between alignment_wtrial and alignment_wcaster.

Both tools print a "Species weights for alignment file:" block in their logs (debug dump).
We align by species name and report per-species difference. Pass criterion: |Δ| < 1e-9.
"""
import re, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
RES = HERE / 'results'

def parse_weights(log):
    """Return dict: {species_name: weight}. Reuse weight lines '  <name> = <value>'."""
    w = {}
    for line in open(log):
        m = re.search(r'^\s*-\s*(\S+)\s*=\s*([0-9.eE+\-]+)\s*$', line)
        if m:
            name, val = m.group(1), float(m.group(2))
            if name in ('weights',):  # no-op guard
                continue
            w[name] = val
    return w

def write_out(wt, wc):
    out = HERE / 'results' / 'weight_compare.txt'
    names = sorted(set(wt) | set(wc))
    maxdiff = 0.0
    with open(out, 'w') as f:
        f.write(f"{'species':<28} {'alignment_wtrial':>18} {'alignment_wcaster':>20} {'diff':>16}\n")
        f.write('-' * 88 + '\n')
        for n in names:
            a = wt.get(n, float('nan'))
            b = wc.get(n, float('nan'))
            d = abs(a - b) if not (a != a or b != b) else float('nan')
            if d == d:
                maxdiff = max(maxdiff, d)
            f.write(f'{n:<28} {a:>18.10f} {b:>20.10f} {d:>16.3g}\n')
        f.write('-' * 88 + '\n')
        f.write(f'max |diff| = {maxdiff:.3g}\n')
    print(open(out).read())
    return maxdiff

wt = parse_weights(RES / 'alignment_wtrial.log')
wc = parse_weights(RES / 'alignment_wcaster.log')

if not wt or not wc:
    print('ERROR: no weight lines found in logs (need --log with weight dump)')
    sys.exit(1)

maxdiff = write_out(wt, wc)
ok = set(wt) == set(wc) and maxdiff < 1e-9
print(f'VERDICT: {"PASS" if ok else "FAIL"}  (species sets equal: {set(wt)==set(wc)}, max|diff|={maxdiff:.3g})')
sys.exit(0 if ok else 1)
