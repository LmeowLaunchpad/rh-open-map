# 01. Window certificates and how we beat Zhu's barrier

## Background
- X. Zhu, arXiv:2608.24827 (Sep 2026): certified Weil positivity Q(f) >= 8.9e-18 ||f||^2 for f supported in [-0.8, 0.8]
  (autocorrelation support 1.6, i.e. window x = e^{1.6} = 4.953), and a **barrier theorem** (Thm 1.4): any certificate that bounds the
  prime comb pointwise needs frequencies up to T1 = 2 pi e^{A}, A = 2 sum_{n<x} Lambda(n)/sqrt(n) ~ 4 sqrt(x), so cost grows doubly
  exponentially; certification past support 3.2 is "computationally void" by such methods.
- Frozen-symbol engine (agent E4): W >= W_low, the archimedean symbol frozen at its value beyond T0 (it is increasing in |t|), the
  remainder an entire-kernel compact operator plus a rank-2 pole; Legendre-Galerkin in Arb balls with rigorous tail and quadrature
  bounds; certified by ball LDL^T. Certified x = 3 (T0 = 120), x = 5 (T0 = 600), x = 7 (T0 = 900) [proved, computer-assisted, independent
  reruns]. x = 8 would need T0 ~ 3600 (Zhu's wall).

## The invisible-harmonic certificate [proved]
For a prime p with powers p, ..., p^K below x, its part of the comb is a trigonometric polynomial of degree K in theta = t log p:
F_p(theta) = -2 sum_{k<=K} c_{p^k} cos(k theta). Harmonics of order > K have t-frequencies k log p >= log p^{K+1} > log x: they integrate
to zero against |fhat|^2 (whose inverse transform lives in [-log x, log x]). So the comb may be replaced by -C_p + mu_p with mu_p >= 0
atomic, where by Caratheodory-Toeplitz the best constant is C_p = -lambda_min of the (K+1)x(K+1) Toeplitz matrix [hat F_p(i-j)] (the
extremal measure sits on the Pisarenko support). Then exactly

    W(f) = (1/2pi) int |fhat|^2 (M_inf - C_x) dt + (1/2pi) sum_p int |fhat|^2 dmu_p + pole,   C_x = sum_p C_p,

and the frozen-symbol skeleton applies with alpha' = M_inf(T0) - C_x: everything except alpha' I stays smoothing (finitely many band atoms
are rank-one terms). For sqrt(x) < p < x (K = 1) the constant is exactly half the pointwise value. This is outside the class covered by
Zhu's Thm 1.4 (it bounds comb + invisible harmonics, not the comb), so it does not contradict it.

| x | A_x (pointwise) | C_x (this method) | T0 needed, pointwise (alpha = 0.583) | T0 needed, this method |
|---|---|---|---|---|
| 5 | 2.94 | 1.52 | 213 | 52 |
| 7 | 4.38 | 2.24 | 900 | 106 |
| 11 | 7.08 | 3.77 | 13,305 | 487 |
| 13 | 8.52 | 4.49 | 56,493 | 1,003 |
| 17 | 10.29 | 5.47 | 331,446 | 2,673 |
| 31 | 16.64 | 8.93 | 189,772,404 | 85,084 |

In practice T0 is also limited by the freezing loss on the deepest ladder level (alpha' ~ 1 needed).

## Certified runs (this method)
| Window | T0 | N (Legendre modes) | bits | certified lowest levels (even / odd) | status |
|---|---|---|---|---|---|
| x = 7 | 220 | 400 | 256 | 4.4e-28 / 1.1e-24 | [proved, computer-assisted], awaiting verification |
| x = 11 | 800 | 1280 | 384 | 6.7e-49 / 4.2e-45 | [proved, computer-assisted], awaiting verification |
| x = 13 | 1700 | 2700 | 448 | 2.1e-59 / 1.9e-55 | [proved, computer-assisted], awaiting verification |

x = 11: quadrature error unit 1.0e-75, tail eps 1.4e-73, alpha' = 1.0797, ladder K = 8 per sector. The certified levels sit just
below the true ones measured independently (true even ground 1.86e-48, odd 8.22e-45 at x = 11). At T0 = 550 the even sector certified
but the odd sector failed at -3.7e-40 (freezing loss), fixed by T0 = 800.
Reproduce: `NOCERT=1 NSMALL=80 python wcert_h.py 11 1280 800 0 0 384 15 5 200 3` then
`python cert2.py GH_x11_N1280_T800.pkl tailH_x11_N1280_T800.pkl 1e-7 16` (code/certificates/).

## Krein optimality [proved identification, checked values]
Positivity on a window is equivalent (Krein) to the existence of a positive measure whose Fourier transform agrees with the Weil
distribution on the window. The best constant obtainable by invisible harmonics of the prime comb alone (jointly over all primes) is
therefore the top of the compressed atom operator, Lambda_x (E4): 1.77 / 5.56 / 11.09 / 34.95 at x = 7 / 30 / 100 / 1000, ~1.1 sqrt(x).
An off-line zero far from the line can break window x from heights up to ~e^{sqrt x - log x} (E3, long-range births), so window
certificates of any kind need frequencies ~e^{c sqrt x}: the improved method is near the intrinsic limit. Windows alone cannot reach RH.

## Related
- A. Groskin, arXiv:2607.02828: certifies (Arb, 9000 bits) that the N = 200 Galerkin compression at x = 100 has no negative eigenvalue
  (a finite-subspace statement, not full-window positivity). Should be cited alongside.
- Gluing with verified zeros (open-problems/P3): an exact splice of the zero side and the prime side [proved]; the missing piece is an
  explicit sampling inequality at the zeta zeros.
