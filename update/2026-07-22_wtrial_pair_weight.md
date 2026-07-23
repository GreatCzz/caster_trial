# wtrial Pair-Weight Refactoring (2026-07-22)

## Mathematical problem

The original wtrial replaced integer nucleotide counts with weighted doubles
(`cnts += weight`) and passed them directly into XXYY.  The XXYY formula
contains terms like `y2 * (y2 - 1)` which are combinatorial: "choose 2
distinct individuals from colour 2, both of type Y".  This is correct only
for **integer counts** (where y2 = number of individuals).

With weighted fractional counts, `y2 * (y2 - 1)` ≠ true weighted-pair sum.
Example: 3 species in colour 2, all purine, weights 0.5, 0.7, 0.3.
  y2 = 1.5  →  y2*(y2-1) = 0.75
  True weighted pairs: 0.5*0.7 + 0.5*0.3 + 0.7*0.3 = 0.71  ≠ 0.75

## Solution

Track weighted-pair sums explicitly in a new `colorPairWeight` array,
updated incrementally during `elementSetOrClearTaxonColor`.

**Add species (weight w, nucleotide count c)**:
```
colorPairWeight[pos][color][nuc] += colorWeight[pos][color][nuc] * (c * w)
colorWeight[pos][color][nuc]     += c * w
```

**Remove species**:
```
colorWeight[pos][color][nuc]     -= c * w
colorPairWeight[pos][color][nuc] -= colorWeight[pos][color][nuc] * (c * w)
```

## Changed functions

| Function | Change |
|----------|--------|
| `elementSetOrClearTaxonColor` | Maintains colorWeight + colorPairWeight instead of colorCnts |
| `XXYY` | Added 4 pair-weight params (`x11,x22,y11,y22`); `*(...-1)` → `*2* pairWeight` |
| `scorePos` | New signature `(colorWeight, colorPairWeight, refCnt, pi)`; computes pair-weights via `pw()` lambda |
| `quadPos` (weight_t overload) | New, for elementQuadripartitionScores |
| `quadPosSingle` (weight_t overload) | New, internal types = weight_t |
| `elementScore` | Reads colorWeight/colorPairWeight instead of colorCnts |
| `elementQuadripartitionScores` | Uses weight_t quadPos |

## The `pw()` lambda (in scorePos)

```cpp
auto pw = [&](int col, int n1, int n2) -> weight_t {
    if (n1 == n2) return cp[col][n1];
    return cp[col][n1] + cp[col][n2] + cw[col][n1] * cw[col][n2];
};
```

Computes the weighted-pair sum for nucleotides `n1` and `n2` in colour group `col`:
- Same nucleotide (e.g. A+A): `cp[col][nuc]` — already tracked pairs.
- Different nucleotides (e.g. A+G for purines): pairs where BOTH are A (from cp), BOTH are G (from cp), plus cross-nucleotide products `cw[A] * cw[G]`.

Used to derive `r11, r22, y11, y22, a11, a22, g11, g22, c11, c22, t11, t22`
which replace the 12 pair-weight arguments to the 9 XXYY calls.

## New data structures

| Name | Type | Size | Purpose |
|------|------|------|---------|
| `colorWeight` | `vector<array<array<weight_t,4>,4>>` | nGenomePos × 4 × 4 | Weighted nucleotide sums per position/colour |
| `colorPairWeight` | `vector<array<array<weight_t,4>,4>>` | nGenomePos × 4 × 4 | Weighted-pair sums per position/colour/nuc |
| `Element::speciesWeights` | `vector<weight_t>` | nTaxa | Per-species-row sequence-similarity weight |
| `weight_t` | `= double` | — | Global type alias for pair-weight arithmetic |

## Removed

- `colorCnts` (old weighted-count array) — replaced by colorWeight + colorPairWeight
- Old `XXYY` signature without pair-weight params
- Old `scorePos` signature without pair-weight params
- Old wtrial's `speciesWeights` applied as `cnts += weight` (integer count replaced by double)

## Performance

| Version | Time | Memory | Math correct? |
|---------|------|--------|---------------|
| Old wtrial | 21s | 615 MB | ❌ |
| New wtrial (pair-weight) | 27s | 615 MB | ✅ |
| CASTER_TRI | 28s | 545 MB | ✅ |
