# Cat Dataset — CASTER / CASTER_TRI / wtrial Comparison

**Date**: 2026-07-22
**Data**: `example/cat/test_full.fasta` (570 MB, 10 species, Felidae family)
**Reference species**: Felis_catus
**Root**: Canis_lupus_familiaris

## Methods

| Tool | Binary | Command |
|------|--------|---------|
| CASTER | `bin/caster` | `-i test_full.fasta -t 4 --initial-round 4 --subsequent-round 2 --root Canis_lupus_familiaris` |
| CASTER_TRI | `bin/caster_tri` | Same params, input via fasta2ref.txt |
| wtrial | `bin/wtrial` | Same params, input via fasta2ref.txt |

## Performance

| Tool | Time | Peak Memory | Binary Size |
|------|------|------------|-------------|
| CASTER | 14s | ~545 MB | 1.4 MB |
| CASTER_TRI | 28s | ~545 MB | 1.2 MB |
| wtrial | 21s | ~601 MB | 544 KB |

wtrial is faster than CASTER_TRI (21s vs 28s) due to single DataClass variant (no type-fallback retries). Both are slower than CASTER (14s) due to R-subtraction scoring overhead. wtrial uses ~56 MB more memory than CASTER_TRI because `cnts` stores `double` (8 bytes) instead of `unsigned char` (1 byte).

## Results

### Topologies — verified with Robinson-Foulds distance

| Pair | RF distance | Conclusion |
|------|------------|------------|
| CASTER vs CASTER_TRI | 0 | Identical |
| CASTER vs wtrial | 0 | Identical |
| CASTER_TRI vs wtrial | 0 | Identical |

RF distance computed via `lib/tree_utils.rf_distance()` using dendropy's symmetric_difference. RF=0 means all three trees share the same unrooted bipartition set — **topologically identical**.
验证拓扑结构命令：cd /home/great_czz/project/caster_trial/example/cat/results_20260722_1 && python3 rf_verify.py caster.tre caster_tri.tre wtrial.tre

| Tool | Tree | Bootstrap (Leopardus/Caracal) |
|------|------|------|
| CASTER | `(((Prionailurus,Otocolobus)100,Felis)100,Lynx)100,...`  | 98.9 |
| CASTER_TRI | `(((Prionailurus,Otocolobus)100,Felis)100,Lynx)100,...`  | 94.2 |
| wtrial | `((Otocolobus,Prionailurus)100,Felis)100,Lynx)100,...`    | 95.8 |

### Bootstrap comparison

| Clade | CASTER | CASTER_TRI | wtrial |
|------|--------|-----------|--------|
| (Prionailurus,Otocolobus) | 100 | 100 | 100 |
| (...Felis) | 100 | 100 | 100 |
| (...Lynx) | 100 | 100 | 100 |
| (...Acinonyx) | 100 | 100 | 100 |
| (...Leopardus) | 98.9 | 94.2 | 95.8 |
| (...Caracal) | 100 | 100 | 100 |
| (...Panthera) | 100 | 100 | 100 |

### Weighted trial diagnostic

| Non-ref species | Hamming distance | nonGapLen | Similarity | Weight |
|------|------|------|------|------|
| (computed per position during read, not logged) | — | — | — | — |

### Observations

1. **All three methods agree** on the tree topology.
2. **wtrial bootstrap on Leopardus/Caracal** (95.8) falls between CASTER (98.9) and CASTER_TRI (94.2), suggesting weighting attenuates the signal difference between methods.
3. CASTER and CASTER_TRI produce identical unrooted topology; minor sibling-order differences are random.
