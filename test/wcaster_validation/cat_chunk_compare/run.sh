#!/bin/bash
# Test2: chunk_wcaster vs chunk_wtrial on real cat data, --chunk 1000, NO weight dump.
# Compare both to each other and to the reference true species tree.
cd "$(dirname "$0")" || exit 1

WTRIAL=../../../bin/chunk_wtrial
WCASTER=../../../bin/chunk_wcaster
FASTA2REF=../../../example/cat/fasta2ref.txt   # relative: test_full.fasta  Felis_catus
ROOT=Canis_lupus_familiaris
CHUNK=1000
RESULTS=results
mkdir -p "$RESULTS"

echo "=== chunk_wtrial ==="
$WTRIAL -i "$FASTA2REF" -t 8 --initial-round 4 --subsequent-round 2 --root "$ROOT" --chunk "$CHUNK" \
  --log "$RESULTS/chunk_wtrial.log" -o "$RESULTS/chunk_wtrial.tre" 2>&1 | tail -2
echo ""
echo "=== chunk_wcaster ==="
$WCASTER -i "$FASTA2REF" -t 8 --initial-round 4 --subsequent-round 2 --root "$ROOT" --chunk "$CHUNK" \
  --log "$RESULTS/chunk_wcaster.log" -o "$RESULTS/chunk_wcaster.tre" 2>&1 | tail -2
echo ""
echo "=== topology comparison ==="
python3 compare_topology.py
