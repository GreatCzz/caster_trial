#!/bin/bash
# Test1: chunk_wcaster vs chunk_wtrial — 8 species, 100kb, per-chunk weight + topology comparison
# Both tools run with the SAME --chunk so chunk boundaries (and per-chunk weights) are directly comparable.
cd "$(dirname "$0")" || exit 1

WTRIAL=../../../bin/chunk_wtrial
WCASTER=../../../bin/chunk_wcaster
CHUNK=1000
RESULTS=results
mkdir -p "$RESULTS"

echo "=== chunk_wtrial ==="
$WTRIAL -i fasta2ref.txt -t 1 --initial-round 4 --subsequent-round 2 --chunk "$CHUNK" \
  --dump-chunk-weights \
  --log "$RESULTS/chunk_wtrial.log" -o "$RESULTS/chunk_wtrial.tre" 2>&1 | tail -2
echo ""
echo "=== chunk_wcaster ==="
$WCASTER -i fasta2ref.txt -t 1 --initial-round 4 --subsequent-round 2 --chunk "$CHUNK" \
  --dump-chunk-weights \
  --log "$RESULTS/chunk_wcaster.log" -o "$RESULTS/chunk_wcaster.tre" 2>&1 | tail -2
echo ""
echo "=== weight comparison ==="
python3 compare_weights.py
echo ""
echo "=== topology comparison ==="
python3 compare_topology.py
