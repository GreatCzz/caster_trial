#!/bin/bash
# Test1: alignment_wcaster (weight->1) vs CASTER — degenerate scoring validation
cd "$(dirname "$0")" || exit 1

CASTER=../../../bin/caster
WCASTER=../../../bin/alignment_wcaster
RESULTS=results
mkdir -p "$RESULTS"

echo "=== CASTER ==="
$CASTER -i test.fasta -t 1 --initial-round 4 --subsequent-round 2 --chunk 800 \
  --log "$RESULTS/caster.log" -o "$RESULTS/caster.tre" 2>&1 | grep -E "Score:|Final tree" | tail -8
echo ""
echo "=== alignment_wcaster ==="
$WCASTER -i fasta2ref.txt -t 1 --initial-round 4 --subsequent-round 2 --chunk 800 \
  --dump-chunk-weights \
  --log "$RESULTS/alignment_wcaster.log" -o "$RESULTS/alignment_wcaster.tre" 2>&1 | grep -E "Species weights|Score:|Final tree" | tail -20
echo ""
echo "=== trees ==="
cat "$RESULTS/caster.tre"; echo
cat "$RESULTS/alignment_wcaster.tre"; echo
