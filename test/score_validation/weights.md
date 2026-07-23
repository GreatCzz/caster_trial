# test.fasta weight analysis

## Sequences
ref: ACGTACGTACGTACGTACGT
1:   AGCTACGTACGTACGTACGT
2:   ATCGACGTACGTACGTACGT
3:   TCGAACGTACGTACGTACGT

## Per-species weights

| Species | Hamming | nonGapLen | Similarity | Weight |
|---------|---------|-----------|------------|--------|
| ref     | —       | —         | —          | 1.000000 |
| 1       | 2       | 20        | 0.900      | 0.866667 |
| 2       | 3       | 20        | 0.850      | 0.800000 |
| 3       | 2       | 20        | 0.900      | 0.866667 |

Formula: similarity = 1 - hamming/nonGapLen
         weight = 0 (if sim < 0.25), else (sim - 0.25) / 0.75

## Weight product (4-species quartet)

ref × 1 × 2 × 3 = 1.0 × 0.866667 × 0.800000 × 0.866667 = 0.600889

## Mismatch positions vs ref

pos  ref  1   2   3
  0   A   A   A   T  ← 3 mismatches ref
  1   C   G   T   C  ← 2 mismatches ref
  2   G   C   C   G  ← 1 mismatch ref
  3   T   T   G   A  ← 2 mismatches ref
  4   A   A   A   A  ← all match
  5   C   C   C   C
  6   G   G   G   G
  7   T   T   T   T
(positions 4-19 all match)
