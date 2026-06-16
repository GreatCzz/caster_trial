# caster_trial Project

Phylogenetic species tree inference tool series. Contains ASTRAL, CASTER, SISTER, and CASTER_TRI.

## Build Commands

```bash
make              # Build all (astral, caster, sister, caster_tri)
make caster_tri   # Build CASTER-TRI only
```

- Compiler: g++ >= 13
- Standard: C++20 with `-std=c++20 -march=native -Ofast`
- All tools share `src/driver.cpp`; selected via preprocessor define (`-D ASTRAL`, `-D CASTER`, `-D SISTER`, `-D CASTER_TRI`)

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

## Test Data
```
example/
├── test_small.fasta           # Single-ref test (Felis_catus, 500K bp × 10 species)
├── fasta2ref.txt              # Single-ref fasta2ref
└── test_caster_tri/
    ├── test_small/
    │   └── test_small.fasta + fasta2ref.txt
    └── test_multi_ref/
        ├── test_small_felis.fasta + test_small_otocolobus.fasta
        └── fasta2ref.txt (Felis_catus + Otocolobus_manul)
```
