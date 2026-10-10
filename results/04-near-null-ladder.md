# 04. The near-null ladder and the boundary-layer law

## The near-null space is the Connes-Consani prolate structure [checked]
- Radical vectors E(h)(u) = e^{u/2} sum_m h(m e^u) (h even, band-limited to [-sqrt x, sqrt x]) vanish exactly left of the window
  (Poisson) and are killed by the Weil form; truncated to the window, **W(chi E(h)) = W(t_h)**, the Weil energy of the exterior tail
  [proved modulo admissibility, E1]. The near-null energies are time-leakage (Slepian) energies at c = 2 pi x / q.
- Parity rule: E(psi_n) is ~even for n = 0 mod 4, ~odd for n = 2 mod 4; one effective constraint per sector; even and odd levels
  interleave exactly (x = 6: 8.6e-23 < 1.7e-19 < 1.6e-16 < ...). **Credit**: this bookkeeping is in Connes-Consani, arXiv:2106.01715,
  Sec. 3 (rediscovered by agent E1); the quantitative ladder below is new.
- Floor: mu_1 ~ x^b e^{-4 pi x / q} with b ~ n0 + 1/2 + ~0.4 (n0 = first admissible prolate index); measured across families:
  chi_5 (n0 = 0): 0.77; chi_-4 (n0 = 1): 2.06; zeta (pole kills psi_0, n0 = 4): 5.36 even [preregistered, passed]. The pole of zeta raises
  the floor by exactly one constraint.
- Certified window margins follow the same rate: log10 lowest level -7.3 / -16.7 / -27.2 / -48.2 at x = 3 / 5 / 7 / 11 (about 4 pi per
  unit x).

## K ~ 1 is Slepian's identity [checked, RH-neutral]
The Hadamard/Laguerre rate of the ground state (Lemma G of E3: L mu_1' = P(f_1)) stays ~1 for zeta; it equals Slepian's derivative
identity d lambda_n/dc = (2/c) psi_n(1)^2 lambda_n for the prolate image. It is a sharp birth detector in fake worlds (DH: ~1 until just
before its birth, then 1.8 / 3.8 / 15 / 125) but gives no early warning.

## Boundary-layer law L3 (agent E1)
Edge profile of the k-th near-null state: f(L/2 - s)/f(L/2) = e^{-s/2} I_0(beta_k sqrt(1 - e^{-s})), beta_k = sqrt(2 (c^2 - chi_{4k})),
c = 2 pi x / q (the Frobenius solution at the prolate's regular singular endpoint). Fed into the whitened arrival block, it predicts the
full-deficit horizons of the next prime to 0.86-0.90 of the measured ones at n = 5, 7, 11, 13, 17, and the universal arrival law
delta = z^2 / (2 c^2) with the (n/q)^{-2} dependence exact.
Status: [proved] for the block lemmas, [checked] to machine precision for the prolate approximation, [heuristic] for one constant;
the sharp version needs a lower bound off the trial space (i.e. positivity), so it is RH-neutral but conditional as stated.
