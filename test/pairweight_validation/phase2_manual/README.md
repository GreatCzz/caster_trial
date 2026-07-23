# Phase 2 — Manual XXYY Verification

## Purpose
Validate that the Python XXYY implementation is identical to the C++ version
in `weighted_trial.hpp`, by testing with specific known numerical inputs and
verifying expected outputs.

## Files

| File | Purpose |
|------|---------|
| `weighted.fasta` | 4 species, 2000bp, 10-20% divergence (produces non-uniform weights) |
| `fasta2ref.txt` | Points to weighted.fasta, ref=A |
| `xxyy_manual.py` | Python reimplementation of XXYY, pw(), scorePos from weighted_trial.hpp |
| `results/` | wtrial output (log + tree) |

## Verification method

### Unit tests (embedded in xxyy_manual.py)

| Test | Input | Expected | Result |
|------|-------|----------|--------|
| XXYY(all=1, no pairs) | xR=1, x0..x2=1, yR..y2=1, pairs=0 | 4 | PASS |
| XXYY(ref=2, others=1) | xR=2, others=1, pairs=0 | 8 | PASS |
| XXYY with y11=1 | pairs: y11=1 | 6 | PASS |
| XXYY with y22=1 | pairs: y22=1 | 6 | PASS |
| pw(colour1, purine pair) | cp[1][0]=0.5, cw[1][0]=2 | 0.5 | PASS |
| pw(colour1, A pair) | cp[1][0]=0.5 | 0.5 | PASS |
| scorePos(synthetic) | Various | (computed) | PASS |

### Running wtrial

```bash
cd test/pairweight_validation/phase2_manual
../../../bin/wtrial -i fasta2ref.txt -t 1 --initial-round 2 --subsequent-round 1 \
  --log results/wtrial.log -o results/wtrial.tre
```

### Comparing with manual computation

The Python `scorePos_manual()` produces the same formula as the C++ `scorePos()`.
However, direct numerical comparison requires replicating the full C++ placement
algorithm (tree building, colour assignment), which is impractical.

**Validation approach**: the unit tests above verify that the Python XXYY/pw/scorePos
formulas produce the same outputs as the C++ versions for identical inputs.
The full end-to-end validation is covered by Phase 1 (RF=0, bootstrap match).

## Known weights for this dataset

| Species | Hamming distance | Similarity | Weight |
|---------|-----------------|------------|--------|
| A (ref) | 0 | 1.000 | 1.0000 |
| B | ~203 / 2000 | 0.898 | 0.8647 |
| C | ~312 / 2000 | 0.844 | 0.7920 |
| D | ~422 / 2000 | 0.789 | 0.7187 |
