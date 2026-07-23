#!/usr/bin/env python3
"""Extract placement scores, NNI scores, and bootstrap values from CASTER_TRI / wtrial logs and trees."""

import sys, re
from pathlib import Path

def parse_log(logfile: str) -> dict:
    """Extract all Score: lines from the log file. Returns {stage: score}."""
    scores = []
    with open(logfile) as f:
        for line in f:
            m = re.search(r'^\s*-?\s*Score:\s+([0-9.eE+\-]+)', line)
            if m:
                scores.append(float(m.group(1)))
    return scores

def parse_bootstrap(treefile: str) -> list:
    """Extract all bootstrap values from a Newick string (numbers after ')')."""
    with open(treefile) as f:
        newick = f.read().strip()
    # Match )<number>: or )<number>,
    return [float(m) for m in re.findall(r'\)([0-9]+(?:\.[0-9]+)?)[,:\)]', newick)]

def parse_all_scores(logfile: str) -> dict:
    """Separate placement scores (scorePos) from NNI scores (quadPos).
    In the log, placement scores appear first (scorePos via elementScore),
    then NNI scores appear after "Performing NNI search" (quadPos via elementQuadripartitionScores).
    """
    scores = []
    nnistarted = False
    placement = []
    nni = []
    with open(logfile) as f:
        for line in f:
            if 'Performing NNI search' in line:
                nnistarted = True
            m = re.search(r'^\s*-?\s*Score:\s+([0-9.eE+\-]+)', line)
            if m:
                val = float(m.group(1))
                scores.append(val)
                if nnistarted:
                    nni.append(val)
                else:
                    placement.append(val)
    return {'all': scores, 'placement': placement, 'nni': nni}

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <c_tri_log> <wtrial_log>")
        sys.exit(1)

    ctri = parse_all_scores(sys.argv[1])
    wt = parse_all_scores(sys.argv[2])
    ctree_bs = parse_bootstrap(sys.argv[1].replace('.log', '.tre'))
    wtree_bs = parse_bootstrap(sys.argv[2].replace('.log', '.tre'))

    print("=== CASTER_TRI ===")
    print(f"  placement scores: {len(ctri['placement'])}  [{', '.join(f'{s:.2f}' for s in ctri['placement'][:5])}...]")
    print(f"  NNI scores:       {len(ctri['nni'])}  [{', '.join(f'{s:.2f}' for s in ctri['nni'][:5])}...]")
    print(f"  bootstrap:        {ctree_bs}")

    print("\n=== wtrial ===")
    print(f"  placement scores: {len(wt['placement'])}  [{', '.join(f'{s:.2f}' for s in wt['placement'][:5])}...]")
    print(f"  NNI scores:       {len(wt['nni'])}  [{', '.join(f'{s:.2f}' for s in wt['nni'][:5])}...]")
    print(f"  bootstrap:        {wtree_bs}")
