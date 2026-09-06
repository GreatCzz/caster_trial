#!/usr/bin/env python3
"""Generate Test2 quartet data: 4 taxa (ref + 3), ~2000bp genuinely diverse.

Goal: with only 4 taxa the sole quartet necessarily contains ref, so wcaster's global
scoring must EQUAL alignment_wtrial's ref-restricted scoring on the same topology.
Weights are real (not 1) to exercise the weighted XXYY/scorePos math.
"""
import random

random.seed(20260904)

N = 2000
TAXA = ['ref', 'a', 'b', 'c']

def rand_seq(n):
    return ''.join(random.choice('ACGT') for _ in range(n))

def mutate(seq, rate):
    return ''.join(c if random.random() >= rate else random.choice('ACGT') for c in seq)

# ref: random; a/b: close to ref; c: more divergent (all share informative block, no A tail)
ref = rand_seq(N)
a = mutate(ref, 0.06)   # ~94% identity to ref
b = mutate(ref, 0.07)
c = mutate(ref, 0.18)   # ~82% identity to ref

# ensure a and b are closer to each other (share some derived state) for a clean (a,b)|(ref,c) split
# (topology is not asserted; both tools must simply agree)
seqs = {'ref': ref, 'a': a, 'b': b, 'c': c}

with open('test.fasta', 'w') as f:
    for t in TAXA:
        f.write(f'>{t}\n{seqs[t]}\n')

with open('fasta2ref.txt', 'w') as f:
    f.write('test.fasta    ref\n')

print('expected similarity/weight (nonzero, moderate):')
for t in TAXA[1:]:
    h = sum(1 for i in range(N) if seqs[t][i] != ref[i])
    sim = 1.0 - h / N
    w = 0.0 if sim < 0.25 else (sim - 0.25) / 0.75
    print(f'  {t}: hamming={h} sim={sim:.4f} weight={w:.4f}')
