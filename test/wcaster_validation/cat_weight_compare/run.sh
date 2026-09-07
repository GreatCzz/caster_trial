#!/bin/bash
# cat_weight_compare — compare per-species weights between alignment_wtrial and alignment_wcaster
# Uses the real cat dataset at example/cat (ref = Felis_catus).
cd "$(dirname "$0")" || exit 1

WTRIAL=../../../bin/alignment_wtrial
WCASTER=../../../bin/alignment_wcaster
FASTA2REF=../../../example/cat/fasta2ref.txt   # relative: test_full.fasta  Felis_catus
ROOT=Canis_lupus_familiaris
RESULTS=results
mkdir -p "$RESULTS"

echo "=== alignment_wtrial ==="
$WTRIAL -i "$FASTA2REF" -t 8 --initial-round 4 --subsequent-round 2 --root "$ROOT" \
  --dump-chunk-weights \
  --log "$RESULTS/alignment_wtrial.log" --log-verbose 5 -o "$RESULTS/alignment_wtrial.tre" 2>&1 | tail -2
echo ""
echo "=== alignment_wcaster ==="
$WCASTER -i "$FASTA2REF" -t 8 --initial-round 4 --subsequent-round 2 --root "$ROOT" \
  --dump-chunk-weights \
  --log "$RESULTS/alignment_wcaster.log" --log-verbose 5 -o "$RESULTS/alignment_wcaster.tre" 2>&1 | tail -2
echo ""
echo "=== weight comparison ==="
python3 compare_weights.py
