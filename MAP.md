# The map

## 0. The object

Weil's criterion: RH holds iff the Weil functional is nonnegative on all test functions g = f * f~. Restricting f to [-L/2, L/2]
(L = log x) gives the **window form** W_x (E4's convention):

    W(f) = (1/2pi) int |fhat(t)|^2 M_inf(t) dt  -  sum_{n<x} (Lambda(n)/sqrt n) (g(log n) + g(-log n))  +  ghat(i/2) + ghat(-i/2),
    M_inf(t) = Re psi(1/4 + it/2) - log pi.

By the explicit formula W(f) = sum over zeros of ghat(gamma). Window positivity is monotone in x [proved, E1], so
**RH <=> W_x >= 0 for every x**, and certifying a window x certifies every smaller one.

## 1. Components and status

| # | Component | Statement (short) | Status | Page |
|---|---|---|---|---|
| A1 | Window certificates | W_x >= 0 for x <= 13 | [proved, computer-assisted], awaiting verification | results/01 |
| A2 | Invisible-harmonic method | replace each prime's comb by -C_p + positive atoms (Caratheodory-Toeplitz); harmonics above the window are invisible | [proved] | results/01 |
| A3 | Krein optimality | best certificate constant = compressed atom norm Lambda_x ~ 1.1 sqrt x; near the off-line-zero detection limit | [proved] / [checked] | results/01 |
| B1 | Positivity sieve | positivity + Gamma_R + pole re-derive prime powers 2..49 with von Mangoldt weights | [checked] | results/02 |
| B2 | Maximality | positivity forbids prime weights above von Mangoldt | [proved] (E2) | results/02 |
| B3 | Lemma U' | no positive prime measure at finite total variation from von Mangoldt other than von Mangoldt | [proved] (E2, H1) | results/02 |
| B4 | Coarse deficit lemma | positivity forbids weights below von Mangoldt after a delay | [proved mod admissibility] (E2) | results/02 |
| C1 | Arrival Lemma | a new prime power enters only as a boundary coupling B_n = 2 R_b on a thin edge strip | [proved] (H1) | results/03 |
| C2 | Universal arrival family | the directions each arrival must repair form one universal family in a layer coordinate z | [checked] out of sample | results/03 |
| D1 | Near-null ladder | near-null space = Connes-Consani prolate images; parity rule; floor rate 4 pi x | [checked]; parity rule is in Connes-Consani | results/04 |
| D2 | Boundary-layer law (L3) | universal arrival law from the prolate edge profile, no fitted profile | [proved] parts + [checked]; sharpness conditional | results/04 |
| D3 | K ~ 1 | Hadamard/Laguerre rate of the ground state = Slepian's identity | [checked], RH-neutral | results/04 |
| E1 | Growth = RH content | per window: an angle gamma(x) < 1 plus a positive bulk | [heuristic] equivalence, [checked] data | results/05 |
| E2 | Pair-locking | gamma(x) is eaten by tight zero pairs just inside the horizon 4 pi x | [checked], out of sample | results/05 |
| E3 | Causality | window data up to x never contains what window x + dx needs | [heuristic], sharp (long-range births) | results/05 |
| E4 | Phase conspiracies | the prime phases align pointwise thousands of times (in zero gaps) but never at window resolution | [checked] | results/05 |
| F1 | Existence + rigidity | RH <= ARCHPOS + Rig | [proved] chain (E2) | results/06 |
| F2 | ARCHPOS | one inequality on Gamma_R and the pole; RH => ARCHPOS | open; tight [checked] | results/06 |
| F3 | Rig | every positive pair has prime measure = von Mangoldt | [proved] at finite TV; first steps by certificate | results/06 |
| G1 | Gluing lemma | exact splice of zero side (verified zeros) and prime side | [proved] | results/01 + open-problems/P3 |
| H  | Earlier phases | common-signal class (quasi-RH 7/8), its wall, four requirements | [proved]/[checked] | results/07 |

## 2. The shape of the difficulty (one paragraph)

Everything that also holds for fake L-functions (Davenport-Heilbronn, zeta(s-a)zeta(s+a), Epstein, Gonek's zeta_X) is explained by one
archimedean object, the Connes-Consani prolate structure: the sieve, the arrival family, the near-null ladder, K ~ 1, the certificates.
Everything RH-sensitive sits at one place per window: the height where the zero density equals the window's resolution (~2 pi x), where
tight zero pairs eat the margin. Finite-window certificates cannot reach infinity (Zhu's barrier, and our Krein optimality shows the
improved certificates are near the limit). The step must use a global input built from the Euler product, the functional equation and
unitary local factors together; each fake keeps two of the three and fails. The cleanest open targets are in [open-problems/](open-problems/README.md).

## 3. Fake worlds used as kill tests

| Fake | Keeps | Lacks | Behaviour |
|---|---|---|---|
| Davenport-Heilbronn | functional equation, Gamma factor | Euler product | windows go negative from x ~ 31 (births at x = 30.989 even, 31.616 odd) |
| L_a = zeta(s-a) zeta(s+a) | Euler product, functional equation, positive weights | unitary local factors | all zeros off the line; births predicted to 0.5% |
| Epstein zeta | functional equation | Euler product | births tracked to x = 80 |
| Gonek's zeta_X | Euler product (truncated) | functional equation | used as a control |
