# 06. Existence and rigidity: RH <= ARCHPOS + Rig

## Definitions (agent E2)
S := 2 cosh(u/2) + W_inf (the archimedean + pole distribution). A **positive pair** is (pi, nu), pi >= 0 on (0, inf), nu >= 0 on R,
with  int hhat dnu = <S, h> - 2 int h dpi  for all even h in C_c^inf. Zeta's pair is (von Mangoldt, zeros); it is positive iff RH.

- **ARCHPOS**: A(lambda) := 2 int lambda cosh(u/2) du + (1/2pi) int M_inf lambdahat dt >= 0 for every even lambda >= 0 with
  lambdahat >= 0. It involves no primes and no zeros. RH => ARCHPOS [proved]. ARCHPOS holds unconditionally on the Poitou cone
  (Poitou/Odlyzko mechanism) [proved]; on the full doubly positive cone it is open. A single violating lambda would disprove RH.
- **Rig**: every positive pair has pi = von Mangoldt. [proved] at finite total variation (Lemma U').
- **Chain** [proved, modulo tightness bookkeeping]: ARCHPOS + Rig => RH.

## Results of this round
- **Reformulation** [proved]: Rig <=> the only positive measure nu with nu-check <= S on (0, inf) (plus the matching at u = 0) is the
  zero measure of zeta.
- **Local rigidity** [proved, elementary]: no positive pair has nu = zeros + (nonzero smooth density >= 0): it would force the prime
  measure negative on (0, log 2), where von Mangoldt has no mass.
- **No simple smooth pair** [checked]: the family pi = s(u) chi(u) (all positive mass of S beyond a smooth cutoff) leaves the zero side
  negative, best min -0.054 near t ~ 0-2.7. Consistent with E2's LP: near-feasible pairs pile prime mass onto log n with von Mangoldt
  masses (0.491 / 0.629 / 0.358 / 0.697 near log 2/3/4/5 vs 0.490 / 0.634 / 0.347 / 0.720).
- **ARCHPOS is tight** [checked]: minimising A over doubly positive lambda (exact positive-definiteness, lambda(0) = 1) stays positive
  at every support tested but decreases with resolution (x = e: 0.0029 -> 0.00096 from step 0.02 to 0.01), independently recomputed on
  the frequency side. The near-minimiser balances pole +1.115 against archimedean -1.114 with lambdahat parked below and between the first
  zeros. With slack, Rig would be false (a second pair with a uniform zero density); tightness is what Rig predicts. Consequence:
  any proof of ARCHPOS has to locate the zeros, so ARCHPOS (with Rig) is RH-strength.
- **Rig, first step by certificate** [checked, not yet rigorous]: for doubly positive h >= 1 on (0, u1], every positive pair has
  pi((0, u1]) <= A(h)/2. LP optima give <= 0.00093 on (0, 0.3] and <= 0.039 on (0, 0.6] (still decreasing with resolution).

## What it means
If Rig is proved unconditionally (an RH-neutral rigidity theorem: Gamma factor + pole + positivity force the primes), RH becomes
equivalent to ARCHPOS, a single inequality about the Gamma function and the pole. That reduction would be a significant theorem by
itself. ARCHPOS itself is tight and RH-strength.
