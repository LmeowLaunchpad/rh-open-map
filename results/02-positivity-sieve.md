# 02. The positivity sieve: positivity regenerates the primes

## Statement [checked]
Treat the prime-power weights as unknowns. Keep only the archimedean data (Gamma_R and the pole) and demand that the window form be
positive on every window. Then:

- **One step at a time, true past fixed** (agent H1, fast double-precision engine, basis N >= 2 x log x): positivity alone places every
  prime power 5..49 at its exact location with weight Lambda(n)/sqrt(n) to 6-10 digits (prime powers 8, 9, 16, 25, 27, 32, 49 with
  Lambda(n) = log p, not log n), and rejects all 25 composites in between. Examples: 43 found at 42.9999992 with weight 0.573577
  (true 0.573578); 49 = 7^2 at 48.9999997 with weight 0.2779871 = log 7 / 7. The first steps (2, 3, 4) were found with weights
  log 2/sqrt 2, log 3/sqrt 3, Lambda(4)/2.
- **Jointly, nothing fixed** (all weights free, windows 2.9 -> 31.99): atoms within ~2e-3 of 2, 3, 4, 5, 7 with weights within ~1%;
  optimiser-limited (lambda_min ~ -2.7e-4), so this is the crude version.
- **Full measure class** (agent E2, log-barrier SDP, x = 2.9): any nonnegative prime measure on (0, L] that keeps the window positive
  must put a lump of mass ~0.49 within ~0.03 of log 2 (RMS spread <= 0.026), mass <= 0.0027 on (0, 0.6): positivity forces atoms.
- **The squeeze**: the admissible set of weights shrinks onto von Mangoldt like e^{-3x}; von Mangoldt sits at the TOP (excess corner)
  of every admissible interval (e.g. at x = 3 the log 2 weight is admissible in [0.489211, 0.490131]; the true 0.4901291 is 2e-6 below
  the cap). The cap tightens faster than the floor (two-rate squeeze, E1).

## Proved pieces
- **Maximality** [proved, E2]: no admissible prime measure exceeds von Mangoldt (no excess).
- **Lemma U'** [proved, E2 and H1 independently]: a positive prime measure at finite total variation from von Mangoldt equals it.
- **Coarse deficit lemma** [proved modulo admissibility, E2]: weights below von Mangoldt become inadmissible after a delay.
  Spectral reason: window-truncated radical vectors chi E(h) (Connes-Consani) probe the weight at n through
  B_n(E(h)) = (1/pi) int |zeta(1/2+it) htilde(1/2+it)|^2 cos(t log n) dt, which takes both signs.

## What it means
Uniqueness is essentially settled and **RH-neutral** (fake worlds obey it with their own data). It does not by itself prove RH: the
missing half is existence (that von Mangoldt itself is admissible on every window), which is RH. See results/06 for the existence +
rigidity architecture. Novelty: no prior source for the sieve was found by the team's credit search; a proper literature check is
still owed before claiming novelty.
