#!/usr/bin/env python3
"""
rf_verify.py — Verify topological identity of CASTER / CASTER_TRI / wtrial trees.
Uses dendropy's symmetric_difference via lib/tree_utils.

Usage:  python3 rf_verify.py <caster.tre> <caster_tri.tre> <wtrial.tre>
Output: pairwise Robinson-Foulds distances (RF=0 = identical topologies)
"""

import sys
from pathlib import Path

# Add lib/ to path for tree_utils
SYS_LIB = Path(__file__).resolve().parents[4] / 'lib'
sys.path.insert(0, str(SYS_LIB))
from tree_utils import rf_distance

def pairwise(files: list[str]) -> int:
    trees = {}
    for f in files:
        p = Path(f)
        trees[p.stem] = p.read_text().strip()

    keys = list(trees.keys())
    pairs = [(keys[0], keys[1]), (keys[0], keys[2]), (keys[1], keys[2])]
    all_zero = True
    for a, b in pairs:
        d = rf_distance(trees[a], trees[b])
        print(f'{a} vs {b}: RF={d}')
        if d != 0:
            all_zero = False

    print()
    if all_zero:
        print("VERDICT: All three trees are topologically identical (RF=0).")
    else:
        print("VERDICT: Trees DIFFER — non-zero RF distance detected.")
    return 0 if all_zero else 1

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print(f"Usage: {sys.argv[0]} caster.tre caster_tri.tre wtrial.tre")
        sys.exit(1)
    sys.exit(pairwise(sys.argv[1:4]))
