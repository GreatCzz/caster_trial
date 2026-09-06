#!/usr/bin/env python3
"""Generate Test1 degenerate data: 7 taxa, ~800bp diverse block + 1,000,000bp conserved A tail.

Goal: per-species similarity vs ref ~ 1 - (mismatches in block)/(huge tail) -> weight -> 1,
so alignment_wcaster scoring must approach integer CASTER scoring.
Chunk = 800 so the diverse block is its own element with balanced eqFreqs.
"""
import random

random.seed(20260904)

N_BLOCK = 800
N_TAIL = 1_000_000
TAXA = ['ref', 's1', 's2', 's3', 's4', 's5', 's6']

def rand_seq(n):
    return ''.join(random.choice('ACGT') for _ in range(n))

# ref: random diverse block (informative) + long conserved A tail
ref = rand_seq(N_BLOCK) + 'A' * N_TAIL

seqs = {'ref': ref}
for t in TAXA[1:]:
    # each species: independent random diverse block + identical A tail (tail fully conserved vs ref)
    seqs[t] = rand_seq(N_BLOCK) + 'A' * N_TAIL

with open('test.fasta', 'w') as f:
    for t in TAXA:
        f.write(f'>{t}\n{seqs[t]}\n')

with open('fasta2ref.txt', 'w') as f:
    f.write('test.fasta    ref\n')

# report expected weights
print('expected similarity/weight (approx, random block ~ 600/800 mismatches vs ref):')
for t in TAXA[1:]:
    h = sum(1 for i in range(N_BLOCK) if seqs[t][i] != ref[i])
    sim = 1.0 - h / (N_BLOCK + N_TAIL)
    w = 0.0 if sim < 0.25 else (sim - 0.25) / 0.75
    print(f'  {t}: hamming={h} sim={sim:.6f} weight={w:.6f}')
print('total length:', len(ref))
