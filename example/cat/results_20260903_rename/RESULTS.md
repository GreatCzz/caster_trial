# Cat Dataset — Rename Verification (2026-09-03)

**Purpose**: Verify the Phase-1 rename (caster_tri→trial, weighted_trial→alignment_wtrial, chunk_weighted_trial→chunk_wtrial) introduced no functional regression.

**Data**: `example/cat/test_full.fasta` (570 MB, 10 species, Felidae family)
**Reference species**: Felis_catus
**Root**: Canis_lupus_familiaris

## Methods

| Tool (renamed) | Command |
|----------------|---------|
| `bin/trial` | `--input fasta2ref.txt --output trial.tre -t 8 --initial-round 4 --subsequent-round 2 --root Canis_lupus_familiaris` |
| `bin/alignment_wtrial` | same params |
| `bin/chunk_wtrial` | same params |

Run from `example/cat/` (fasta2ref points to `test_full.fasta`). Old reference trees from `results_20260722_1/`.

## Results — Robinson-Foulds verification (`compare.py`)

| Pair | RF |
|------|----|
| trial (new) vs caster_tri (old ref) | **0** |
| alignment_wtrial (new) vs wtrial (old ref) | **0** |
| chunk_wtrial (new) vs caster (old ref) | **0** |
| chunk_wtrial (new) vs caster_tri (old ref) | **0** |
| trial vs alignment_wtrial (new) | 0 |
| trial vs chunk_wtrial (new) | 0 |
| alignment_wtrial vs chunk_wtrial (new) | 0 |

**VERDICT: All renamed binaries match references (RF=0).**

## Bootstrap comparison (Leopardus/Caracal clade)

| Tool | Old reference | Renamed (2026-09-03) |
|------|--------------|----------------------|
| trial (=caster_tri) | 94.2 | 94.2 |
| alignment_wtrial (=wtrial) | 95.8 | 95.8 |
| chunk_wtrial | — | 95.0 |

Topology and support values are reproduced exactly, confirming the rename is behavior-preserving.
