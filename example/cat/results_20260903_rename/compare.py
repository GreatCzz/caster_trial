#!/usr/bin/env python3
"""compare.py — Verify topological identity of renamed binaries vs old reference trees."""
import sys
from pathlib import Path

SYS_LIB = Path(__file__).resolve().parents[4] / 'lib'
sys.path.insert(0, str(SYS_LIB))
from tree_utils import rf_distance

HERE = Path(__file__).resolve().parent
REF = HERE.parent / 'results_20260722_1'

pairs = [
    ('trial (new)', HERE / 'trial.tre', 'caster_tri (old ref)', REF / 'caster_tri.tre'),
    ('alignment_wtrial (new)', HERE / 'alignment_wtrial.tre', 'wtrial (old ref)', REF / 'wtrial.tre'),
    ('chunk_wtrial (new)', HERE / 'chunk_wtrial.tre', 'caster (old ref)', REF / 'caster.tre'),
    ('chunk_wtrial (new)', HERE / 'chunk_wtrial.tre', 'caster_tri (old ref)', REF / 'caster_tri.tre'),
]

def read(p):
    return p.read_text().strip()

all_zero = True
for an, af, bn, bf in pairs:
    d = rf_distance(read(af), read(bf))
    print(f'{an} vs {bn}: RF={d}')
    if d != 0:
        all_zero = False

# three new trees pairwise
new = ['trial.tre', 'alignment_wtrial.tre', 'chunk_wtrial.tre']
print()
for i in range(len(new)):
    for j in range(i + 1, len(new)):
        d = rf_distance(read(HERE / new[i]), read(HERE / new[j]))
        print(f'{new[i]} vs {new[j]} (new): RF={d}')
        if d != 0:
            all_zero = False

print()
if all_zero:
    print("VERDICT: All renamed binaries match references (RF=0).")
else:
    print("VERDICT: MISMATCH — non-zero RF detected.")
sys.exit(0 if all_zero else 1)
