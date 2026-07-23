#!/usr/bin/env python3
"""Phase 1 validation: verify wtrial == CASTER_TRI when all weights are 1."""

import sys, re
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4] / 'lib'))
from tree_utils import rf_distance

def parse_log(logfile: str) -> dict:
    scores, nnistarted = [], False
    placement, nni = [], []
    with open(logfile) as f:
        for line in f:
            if 'Performing NNI search' in line:
                nnistarted = True
            m = re.search(r'^\s*-?\s*Score:\s+([0-9.eE+\-]+)', line)
            if m:
                val = float(m.group(1))
                scores.append(val)
                (nni if nnistarted else placement).append(val)
    return {'placement': placement, 'nni': nni}

def parse_bootstrap(treefile: str) -> list:
    with open(treefile) as f:
        newick = f.read().strip()
    return [float(m) for m in re.findall(r'\)([0-9]+(?:\.[0-9]+)?)[,:\)]', newick)]

RESULTS = Path(__file__).resolve().parent / 'results'
ctri_log = RESULTS / 'caster_tri.log'
wt_log   = RESULTS / 'wtrial.log'
ctri_tre = RESULTS / 'caster_tri.tre'
wt_tre   = RESULTS / 'wtrial.tre'

ctri = parse_log(ctri_log)
wt   = parse_log(wt_log)
ct_bs = parse_bootstrap(ctri_tre)
wt_bs = parse_bootstrap(wt_tre)

# Compare trees
rf = rf_distance(open(ctri_tre).read(), open(wt_tre).read())

# Compare scores (check all placement scores and all NNI scores match)
score_ok = True
# Compare at least the final scores of each round
for i, (p_c, p_w) in enumerate(zip(ctri['placement'], wt['placement'])):
    if abs(p_c - p_w) > 1e-6:
        print(f"  FAIL placement[{i}]: CASTER_TRI={p_c:.6f} wtrial={p_w:.6f}")
        score_ok = False
for i, (n_c, n_w) in enumerate(zip(ctri['nni'], wt['nni'])):
    if abs(n_c - n_w) > 1e-6:
        print(f"  FAIL NNI[{i}]: CASTER_TRI={n_c:.6f} wtrial={n_w:.6f}")
        score_ok = False

# Compare bootstrap
bs_ok = ct_bs == wt_bs

print(f"RF distance:     {rf}  {'PASS' if rf == 0 else 'FAIL'}")
print(f"Placement score: {'PASS' if score_ok else 'FAIL'} (scorePos validation)")
print(f"Boostrap values: {'PASS' if bs_ok else 'FAIL'} (quadPos validation)")
print(f"NNI scores:      {'PASS' if all(abs(a-b)<1e-6 for a,b in zip(ctri['nni'], wt['nni'])) else 'FAIL'}")
print(f"\nCASTER_TRI placement: {ctri['placement'][:3]}")
print(f"wtrial placement:     {wt['placement'][:3]}")
print(f"CASTER_TRI bootstrap: {ct_bs}")
print(f"wtrial bootstrap:     {wt_bs}")
