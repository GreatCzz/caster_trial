#!/usr/bin/env python3
"""Generate Test1 data: 8 species, 100,000 sites, per-chunk weight comparison data.

Structure: concatenate regions of differing divergence so that different chunks get
different per-chunk weights (not all = 1), then run chunk_wcaster / chunk_wtrial with
the SAME --chunk and compare per-(chunk,species) weights + final topologies.

Species tree (rooted at ref; unrooted topology ((a,b),(c,d),(e,f),g,ref) style):
We just simulate divergence from ref directly with block-varying rates; each non-ref
species i gets divergence profile d_i(block) that differs per block so chunk weights differ.
"""
import random

random.seed(20260906)

LEN = 100_000
BLOCK = 1000  # must match --chunk used in run.sh
N_BLOCK = LEN // BLOCK  # 100 chunks

TAXA = ['ref', 'a', 'b', 'c', 'd', 'e', 'f', 'g']

# divergence per species per 1000-site chunk varies sinusoidally in [lo, hi]
# so per-chunk similarity/weight differs across chunks and across species
def rate(i_species, i_chunk):
    base = 0.03 + 0.05 * i_species          # species baseline 0.03..0.38
    wave = 0.5 * (1 + (i_species * 0.3 + i_chunk * 0.12)) % 1.0
    return max(0.01, min(0.6, base * (0.5 + wave)))  # per-site substitution prob

ref = ''.join(random.choice('ACGT') for _ in range(LEN))

seqs = {'ref': ref}
for i, t in enumerate(TAXA[1:], start=1):
    seq = []
    for c in range(N_BLOCK):
        r = rate(i, c)
        seg = ref[c * BLOCK:(c + 1) * BLOCK]
        seq.append(''.join(ch if random.random() >= r else random.choice('ACGT') for ch in seg))
    seqs[t] = ''.join(seq)

with open('test.fasta', 'w') as f:
    for t in TAXA:
        f.write(f'>{t}\n{seqs[t]}\n')
with open('fasta2ref.txt', 'w') as f:
    f.write('test.fasta    ref\n')

print(f'species: {TAXA}')
print(f'length: {LEN}, chunks of {BLOCK}: {N_BLOCK}')
print('per-species mean weight over chunks (approx):')
for t in TAXA[1:]:
    sims = []
    for c in range(N_BLOCK):
        s, t_ = c * BLOCK, (c + 1) * BLOCK
        seg_r, seg_t = ref[s:t_], seqs[t][s:t_]
        h = sum(1 for a, b in zip(seg_r, seg_t) if a != b)
        sim = 1 - h / BLOCK
        sims.append(sim)
    w = sum(max(0.0, (s - 0.25) / 0.75) for s in sims) / len(sims)
    print(f'  {t}: mean weight={w:.3f}, range=[{max(0,(min(sims)-0.25)/0.75):.3f}, {(max(sims)-0.25)/0.75:.3f}]')
