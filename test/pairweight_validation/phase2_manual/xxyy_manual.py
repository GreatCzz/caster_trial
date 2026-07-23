#!/usr/bin/env python3
"""
Phase 2: Manual XXYY comparison against wtrial output.

Compares the C++ XXYY formula (implemented identically in Python)
against wtrial's actual log output for a controlled dataset with known weights.
"""

import sys, re, math
from pathlib import Path

# ─── C++ XXYY equivalent ─────────────────────
def XXYY(xR, x0, x1, x2, yR, y0, y1, y2, x11, x22, y11, y22):
    """Identical to weighted_trial.hpp XXYY (12-arg version)."""
    return (xR * x0 * y1 * y2 * 2 +
            yR * y0 * x1 * x2 * 2 +
            xR * (xR - 1) * y1 * y2 +
            yR * (yR - 1) * x1 * x2 +
            xR * x1 * y22 * 2 +
            yR * y1 * x22 * 2 +
            xR * x2 * y11 * 2 +
            yR * y2 * x11 * 2)

# ─── Pair-weight helper (same as C++ pw lambda) ────
def pw(cp, cw, col, n1, n2):
    """Pair-weight for nucleotides n1,n2 in colour col."""
    if n1 == n2:
        return cp[col][n1]
    return cp[col][n1] + cp[col][n2] + cw[col][n1] * cw[col][n2]

# ─── ScorePos equivalent (same as C++ weighted_trial.hpp scorePos) ────
def scorePos_manual(cw, cp, rCnt, pi):
    """Manual reimplementation of weighted_trial.hpp scorePos."""
    aR, cR, gR, tR = rCnt[0], rCnt[1], rCnt[2], rCnt[3]
    a0, c0, g0, t0 = cw[0][0], cw[0][1], cw[0][2], cw[0][3]
    a1, c1, g1, t1 = cw[1][0], cw[1][1], cw[1][2], cw[1][3]
    a2, c2, g2, t2 = cw[2][0], cw[2][1], cw[2][2], cw[2][3]

    r11 = pw(cp, cw, 1, 0, 2); r22 = pw(cp, cw, 2, 0, 2)
    y11 = pw(cp, cw, 1, 1, 3); y22 = pw(cp, cw, 2, 1, 3)
    a11 = cp[1][0]; a22 = cp[2][0]
    g11 = cp[1][2]; g22 = cp[2][2]
    c11 = cp[1][1]; c22 = cp[2][1]
    t11 = cp[1][3]; t22 = cp[2][3]

    A, Ci, G, T = pi[0], pi[1], pi[2], pi[3]
    R = A + G; Y = Ci + T
    R2 = A*A + G*G; Y2 = Ci*Ci + T*T

    r0 = a0 + g0; y0 = c0 + t0
    r1 = a1 + g1; y1 = c1 + t1
    r2 = a2 + g2; y2 = c2 + t2
    rR = aR + gR; yR = cR + tR

    rryy = XXYY(rR, r0, r1, r2, yR, y0, y1, y2, r11, r22, y11, y22)
    aayy = XXYY(aR, a0, a1, a2, yR, y0, y1, y2, a11, a22, y11, y22)
    ggyy = XXYY(gR, g0, g1, g2, yR, y0, y1, y2, g11, g22, y11, y22)
    rrcc = XXYY(rR, r0, r1, r2, cR, c0, c1, c2, r11, r22, c11, c22)
    rrtt = XXYY(rR, r0, r1, r2, tR, t0, t1, t2, r11, r22, t11, t22)
    aacc = XXYY(aR, a0, a1, a2, cR, c0, c1, c2, a11, a22, c11, c22)
    aatt = XXYY(aR, a0, a1, a2, tR, t0, t1, t2, a11, a22, t11, t22)
    ggcc = XXYY(gR, g0, g1, g2, cR, c0, c1, c2, g11, g22, c11, c22)
    ggtt = XXYY(gR, g0, g1, g2, tR, t0, t1, t2, g11, g22, t11, t22)

    return rryy*R2*Y2 - (aayy+ggyy)*R*R*Y2 - (rrcc+rrtt)*R2*Y*Y + (aacc+aatt+ggcc+ggtt)*R*R*Y*Y


# ─── Weight computation (same as read()) ────
def compute_weights(ref_seq, taxon_seq):
    """Hamming distance to ref, gaps excluded. Returns (similarity, weight)."""
    hamming, nonGap = 0, 0
    n = min(len(ref_seq), len(taxon_seq))
    for i in range(n):
        if ref_seq[i] == '-' or taxon_seq[i] == '-':
            continue
        nonGap += 1
        if ref_seq[i] != taxon_seq[i]:
            hamming += 1
    if nonGap == 0:
        return 0.0, 0.0
    sim = 1.0 - hamming / nonGap
    w = 0.0 if sim < 0.25 else (sim - 0.25) / 0.75
    return sim, w


# ─── Read FASTA ──────────────────────────────
def read_fasta(path):
    """Returns dict {name: seq}."""
    result = {}
    current_name, current_seq = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    result[current_name] = ''.join(current_seq)
                current_name = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        if current_name:
            result[current_name] = ''.join(current_seq)
    return result


# ─── Main ─────────────────────────────────────
if __name__ == '__main__':
    fasta = Path(__file__).resolve().parent / 'weighted.fasta'
    seqs = read_fasta(fasta)
    ref_seq = seqs['A']

    print("=== Weights ===")
    weights = {}
    for name in ['A', 'B', 'C', 'D']:
        if name == 'A':
            weights[name] = 1.0
            print(f"  {name}: weight=1.0 (ref)")
        else:
            sim, w = compute_weights(ref_seq, seqs[name])
            weights[name] = w
            print(f"  {name}: sim={sim:.3f} weight={w:.4f}")

    # ─── Build colorWeight / colorPairWeight for a SAMPLE position ─────
    # At position 0: A=ACGT... so A[0]=A, B[0]=A(w_B), C[0]=A(w_C), D[0]=A(w_D)
    # All species have A at position 0.
    # Suppose tree is ((A,B),(C,D)) and ref=A is in colour 0 (left subtree).
    # After swap: cw[0] = non-ref in ref's group = {B} → cw[0][0] = w_B
    #              cw[1] = other group 1 = {C} → cw[1][0] = w_C
    #              cw[2] = other group 2 = {D} → cw[2][0] = w_D
    # refCnt[0] = 1.0 (A has A at pos 0)
    # cp[0] = 0 (only B in this group, no pairs)
    # cp[1] = 0 (only C in this group)
    # cp[2] = 0 (only D in this group)

    # Simulated colour assignment after tree placement (example)
    # This is an idealised example to test XXYY against itself.
    cw = [[0., 0., 0., 0.] for _ in range(3)]  # 3 colours × 4 nucs
    cp = [[0., 0., 0., 0.] for _ in range(3)]

    wA, wB, wC, wD = weights['A'], weights['B'], weights['C'], weights['D']

    # Assume colour 0 = {B} (non-ref in ref's group), colour 1 = {C}, colour 2 = {D}
    # All have A at position 0
    cw[0][0] = wB  # B has A
    cw[1][0] = wC  # C has A
    cw[2][0] = wD  # D has A
    # No pairs (only 1 species per colour)

    rCnt = [wA, 0, 0, 0]  # ref has A
    pi = [0.25, 0.25, 0.25, 0.25]  # uniform equilibrium frequencies

    score = scorePos_manual(cw, cp, rCnt, pi)
    print(f"\n=== Sample scorePos output ===")
    print(f"  position 0 (all A, weights: B={wB:.4f} C={wC:.4f} D={wD:.4f})")
    print(f"  scorePos = {score:.10f}")

    # ─── Also test with a multi-species colour (to trigger pair-weight) ────
    # Colour 0 = {B1, B2} both have A, weights wB and 0.5
    cw2 = [[0., 0., 0., 0.] for _ in range(3)]
    cp2 = [[0., 0., 0., 0.] for _ in range(3)]

    wB2 = 0.8  # second B species
    cw2[0][0] = wB + wB2
    # pair weight: wB * wB2 (unordered pairs of B and B2 both having A)
    cp2[0][0] = wB * wB2
    cw2[1][0] = wC
    cw2[2][0] = wD

    score2 = scorePos_manual(cw2, cp2, rCnt, pi)
    print(f"\n=== Multi-species test ===")
    print(f"  colour 0: 2 species (w={wB:.4f}, w={wB2:.4f}) both A")
    print(f"  colorPairWeight[0][0] = {cp2[0][0]:.6f}")
    print(f"  scorePos = {score2:.10f}")

    # ─── Read wtrial log ───────────────────
    logfile = Path(__file__).resolve().parent / 'results' / 'wtrial.log'
    if not logfile.exists():
        print(f"\n=== Log not yet generated: {logfile} ===")
        print("Run wtrial first: cd phase2_manual && ../../../bin/wtrial -i fasta2ref.txt -t 1 --initial-round 2 --subsequent-round 1 --log results/wtrial.log -o results/wtrial.tre")
    else:
        print(f"\n=== wtrial log ===")
        with open(logfile) as f:
            for line in f:
                if 'Score:' in line:
                    print(f"  wtrial: {line.strip()}")
