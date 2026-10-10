# 05. Where RH lives: the growth step, the horizon, and the prime phases

## The growth step carries all RH content
Split of the window space: S = span of the prolate images chi E(psi_n) (near-null, energies e^{-4 pi x}), B = its complement (bulk).
Schur: W >= 0 <=> W|B >= 0 and the angle gamma(x) between S and B (strengthened Cauchy-Schwarz) is < 1 [proved, given a bulk gap].
In prolate-whitened coordinates the condition is an O(1) matrix M(x) > 0: no e^{-4 pi x} anywhere [heuristic scaling, checked in the
scalar case]. By "every off-line zero is eventually resolved" [heuristic], RH <=> gamma(x) < 1 for every x.

## Data [checked]
- zeta: gamma^2 = 0.13 / 0.32-0.33 / 0.41 / 0.39 / 0.51 at x = 5 / 10 / 15 / 20 / 25 (finite-N lower bounds; N-converged ~0.414 at
  x = 15). Not saturated.
- Davenport-Heilbronn climbs 0.44 -> 0.59 -> 0.91 -> >1 into its birth.
- **Pair-locking** (E2, out of sample): the coupling sits on the zeros flanking the tightest zero pair just inside the horizon 4 pi x
  (Lehmer pair #34/#35 at x = 10-12; #135/#136 at 0.94 of the horizon for x = 25).
- **Slide test** (E2 and E3, two independent implementations agreeing to 3 digits): sliding one real pair together at x = 10 drives
  gamma^2 to 0.998 while the window stays positive; with one extra "midpoint" trial vector per tight pair the angle stays ~0.25-0.29
  even for an exact double zero, but then the RH content moves into the sign of a small pair block (blind angle).
- Edge margin shrinks slowly: double-zero floor 0.104 -> 0.070 for x = 8 -> 12.

## Causality wall [heuristic, sharp]
Each window resolves zeros up to ~its own horizon; knowing positivity up to x never contains what window x + dx needs. The
long-range exception is sharp in the other direction: an off-line zero far from the line can break window x from heights up to
~e^{sqrt x - log x} (E3; checked in DH, whose off-line zeros near 1e5 make windows negative from x = 125). Verified zeros up to
3e12 therefore give heuristic window positivity only up to x ~ 1e3, not 3e12/4pi.

## The prime phases conspire pointwise, never at resolution [checked]
D_x(t) = 2 sum_{n<x} Lambda(n) n^{-1/2} cos(t log n) exceeds M_inf(t) on 61 / 1518 / 3762 / 8224 intervals for x = 11 / 30.5 / 50.5 /
100.5 (t <= 1e5-2e5), excess up to +9.5, every peak inside an unusually large gap between zeros (1.8-2.4 x the mean spacing).
Reason [proved, standard]: on the critical line Re(zeta'/zeta) = -theta' + pi sum over real zeros of delta (off-line pairs cancel in
the real part), and theta' = M_inf/2. Smoothing at the window's own resolution (Fejer weights) removes every exceedance; the margin is
tightest at the lowest heights (0.017-0.07 at t ~ 20-50). So "the prime phases cannot align" must be read at window resolution, where it
is exactly window positivity.

## Fakes and the necessary ingredients
Each of Euler product, functional equation and unitary local factors is necessary: DH has the functional equation but no Euler
product; zeta(s-a)zeta(s+a) has both but non-unitary local factors; both fail. A proof must use all three jointly.
