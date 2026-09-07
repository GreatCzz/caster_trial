#!/bin/bash
# Test2: alignment_wcaster (global) vs alignment_wtrial (ref-restricted) — 4 taxa, sole quartet contains ref
cd "$(dirname "$0")" || exit 1

WCASTER=../../../bin/alignment_wcaster
WTRIAL=../../../bin/alignment_wtrial
RESULTS=results
mkdir -p "$RESULTS"

echo "=== alignment_wcaster ==="
$WCASTER -i fasta2ref.txt -t 1 --initial-round 4 --subsequent-round 2 --chunk 10000 \
  --dump-chunk-weights \
  --log "$RESULTS/alignment_wcaster.log" -o "$RESULTS/alignment_wcaster.tre" 2>&1 | tail -3
echo ""
echo "=== alignment_wtrial ==="
$WTRIAL -i fasta2ref.txt -t 1 --initial-round 4 --subsequent-round 2 --chunk 10000 \
  --dump-chunk-weights \
  --log "$RESULTS/alignment_wtrial.log" -o "$RESULTS/alignment_wtrial.tre" 2>&1 | tail -3
echo ""
echo "=== trees ==="
cat "$RESULTS/alignment_wcaster.tre"; echo
cat "$RESULTS/alignment_wtrial.tre"; echo
