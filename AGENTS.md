# caster_trial Project

Phylogenetic species tree inference tool series. Contains ASTRAL, CASTER, SISTER, CASTER_TRI, and wtrial.

## Build Commands

```bash
make              # Build all (astral, caster, sister, caster_tri, wtrial)
make caster_tri   # Build CASTER-TRI only
make wtrial       # Build Weighted Trial only
```

- Compiler: g++ >= 13
- Standard: C++20 with `-std=c++20 -march=native -Ofast`
- All tools share `src/driver.cpp`; selected via preprocessor define (`-D ASTRAL`, `-D CASTER`, `-D SISTER`, `-D CASTER_TRI`, `-D WEIGHTED_TRIAL`)

No test framework exists. Validate by comparing output trees against known results.

## Project Structure

```
src/
├── driver.cpp              # Main entry point, dispatches to tool via #ifdef
├── driver.hpp              # DRIVER concept for tool conformance
├── common.hpp              # Shared utilities: LogInfo, InputParser, AnnotatedBinaryTree, Random
├── stepwise_colorable.hpp  # Interface concepts: STEPWISE_COLORABLE, TAXON_ORDER_PRIORITIZING, etc.
├── optimization_algorithm.hpp  # Heuristic search, TaxonOrderGenerator, subsample, placement, NNI, DP
├── placement_algorithm.hpp     # Stepwise color placement algorithm
├── nni_algorithm.hpp           # NNI optimization
├── constrained_dp_algorithm.hpp# Constrained DP tree assembly
├── quadripartition_support.hpp # Quadripartition support annotation
├── alignment_utilities.hpp     # FASTA/Phylip alignment parser
├── threadpool.hpp              # Multi-threaded pool
├── astral.hpp             # ASTRAL tool
├── caster.hpp             # CASTER-site tool
├── caster_tri.hpp         # CASTER_TRI tool (reference-triangulated)
├── weighted_trial.hpp     # wtrial tool (weighted CASTER_TRI)
├── sister.hpp             # SISTER tool
└── documentation.hpp      # DocumentationBase
```

## CASTER_TRI Architecture

### Scoring
- Only scores quartets containing the reference species R
- Uses exact R-subtraction: `score(all_quartets) - score(quartets_without_R)`
- Multi-ref: each element carries its own `iReferenceTaxonId`; `posRefIndex` maps global positions to ref-indexes; per-ref `rColorCnts` stored for R-subtraction.

### Priority Ordering
- Uses `TAXON_ORDER_PRIORITIZING` concept (defined in `stepwise_colorable.hpp`)
- `Color::taxonOrderPrioritizing(order)` is a **static method** that modifies `order` in-place
- All reference species are moved to the front (preserving their shuffled order), followed by non-reference species
- `TaxonOrderGenerator` in `optimization_algorithm.hpp` calls this automatically — no other file changes needed

### External File Changes
- `driver.cpp`: Added `#ifdef CASTER_TRI` block for compilation dispatch
- `makefile`: Added `caster_tri` build target
- No changes to `optimization_algorithm.hpp`, `common.hpp`, or `stepwise_colorable.hpp` — they provide generic hooks used by caster_tri

### Input
- `-i fasta2ref.txt` — fasta2ref file: `<fasta_path> <reference_species>` per line
- Multiple fasta files supported, each with potentially different reference species

## Weighted Trial (wtrial) Architecture

### Fork
Forked from CASTER_TRI v2.6.1 (2026-07-06). Located in `src/weighted_trial.hpp`.

### Scoring
- Identical quartet logic to CASTER_TRI (only quartets containing the reference species R).
- **Weighting**: each non-ref species gets a weight based on Hamming distance to the reference
  in its alignment (gaps excluded):
  ```
  similarity = 1 - hammingDist / nonGapLen
  weight = 0 (if similarity < 0.25), else (similarity - 0.25) / 0.75
  ```
- Weighted nucleotide counts (`cnts += weight` instead of `cnts += 1`) flow through
  unchanged XXYY/quadXXYY formulas.
- Reference species weight is fixed at 1.0.

### DataClass
- **Single variant, double-only**: `cnt_taxon_t = cnt_t = cnt4_t = double`.
- No type-fallback loop needed (5 try-catch → 1 direct call).
- Tradeoff: ~56 MB extra memory vs CASTER_TRI (8-byte `double` vs 1-byte `unsigned char`
  per count), but ~25% faster (no retries).

### Input / Output
- Same fasta2ref format as CASTER_TRI.
- Binary: `bin/wtrial`. Usage: `bin/wtrial -i fasta2ref.txt -o out.tre`.

## Utility Library

- `../lib/tree_utils.py` — dendropy-powered RF distance, tree cleaning, pruning.
  Use `rf_distance(t1, t2)` to verify topological identity (RF=0 = same tree).
- Script template: `rf_verify.py` in each results directory.

## Test Data
```
example/
├── cat/                       # Full-genome cat test (570 MB, 10 Felidae)
│   ├── test_full.fasta
│   ├── fasta2ref.txt (Felis_catus)
│   └── results_*/             # CASTER / CASTER_TRI / wtrial comparison runs
├── test_small/
│   └── test_small.fasta + fasta2ref.txt
└── test_multi_ref/
    ├── test_small_felis.fasta + test_small_otocolobus.fasta
    └── fasta2ref.txt (Felis_catus + Otocolobus_manul)
```
