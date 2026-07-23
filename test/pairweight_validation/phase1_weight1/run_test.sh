#!/bin/bash
# Phase 1: weight=1 validation of wtrial pair-weight scoring
cd "$(dirname "$0")" || exit 1

CASTER_TRI=../../../bin/caster_tri
WTRIAL=../../../bin/wtrial
RESULTS=results

mkdir -p "$RESULTS"

echo "=== CASTER_TRI ==="
$CASTER_TRI -i fasta2ref.txt -t 1 --initial-round 4 --subsequent-round 2 \
  --log "$RESULTS/caster_tri.log" -o "$RESULTS/caster_tri.tre" 2>&1 | grep -E "(Score:|Final tree|score:)" | head -20

echo ""
echo "=== wtrial ==="
$WTRIAL -i fasta2ref.txt -t 1 --initial-round 4 --subsequent-round 2 \
  --log "$RESULTS/wtrial.log" -o "$RESULTS/wtrial.tre" 2>&1 | grep -E "(Score:|Final tree|score:)" | head -20

echo ""
echo "=== Results ==="
echo "CASTER_TRI tree:"
cat "$RESULTS/caster_tri.tre"
echo ""
echo "wtrial tree:"
cat "$RESULTS/wtrial.tre"
