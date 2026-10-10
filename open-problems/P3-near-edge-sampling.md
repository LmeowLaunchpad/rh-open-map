# P3. Sampling at the zeta zeros next to a heavy neighbour

**Setting.** Gluing lemma [proved]: for x not a prime power and 0 < delta < log(X+/x) (X+ the next prime power), any phi with
0 <= phi <= 1 and phi-check supported in [-delta, delta] gives exactly
W(f) = sum_rho phi(gamma) ghat(gamma) + <prime-side certificate, (1 - phi) ghat>. With zeros verified real below H0 = 3e12
(Platt-Trudgian), the errors (unknown zeros above H0, low-band leak) cost only frequency room, never precision.

**Missing statement (Sub-lemma S).** For F of exponential type L/2: sum over verified zeros of phi(gamma)|F(gamma)|^2 + (prime-side
positive part above R) >= eps_x ||F||^2, with log(1/eps_x) allowed as large as ~1e9. The obstruction: functions whose bulk sits just
below the oversampled band (zero density between 1x and ~5x Nyquist), whose band energy is exponentially small and sits next to the
bulk. Max-gap theorems (Groechenig) and disc interpolation need the zeros ~3-10x denser than Nyquist; |F|-based estimates cannot see the
cancellation; interpolation through zeta itself picks up the trivial zeros.

**Payoff.** Rigorous window positivity for every x up to ~100 (LMFDB zero data, ~3e10) or ~150 (3e12 data); ~400+ with P4.
**Wanted.** An explicit Beurling-type sampling theorem for sets with density only slightly above Nyquist, or a quantitative version of
the boundary-layer law (results/04). Six failed attempts are listed in the graveyard.
