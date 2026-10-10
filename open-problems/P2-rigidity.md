# P2. Rigidity at infinite total variation

**Statement.** Every positive pair (pi, nu) for zeta's archimedean data S = 2 cosh(u/2) + W_inf has pi = von Mangoldt.
Equivalently: the only positive measure nu with nu-check <= S on (0, inf) (and the matching at u = 0) is zeta's zero measure.

**Known.** [proved] at finite total variation from von Mangoldt (Lemma U'). [proved] no pair of the form zeros + smooth positive
density. [checked] near-feasible pairs concentrate their prime mass at log n with von Mangoldt masses (E2's LP); positivity forces a lump
of mass ~0.49 near log 2 (E2, window 2.9). First steps by certificate: pi((0, 0.3]) <= 0.00093, pi((0, 0.6]) <= 0.039 for every positive
pair (LP, not yet rigorous).

**Why it matters.** RH-neutral, and a striking theorem on its own ("the Gamma factor, the pole and positivity force the primes").
Together with P1 it gives RH; alone it makes RH equivalent to the single inequality ARCHPOS.

**Entry points.** Make the step-1 certificates rigorous (repair positive-definiteness, evaluate A(h) in Arb), then push the
localisation to log 2, log 3, ... with window (not only doubly positive) tests; the infinite-TV case needs a uniform version of the
squeeze. code/experiments/rig_mass.py, archpos_ac.py.
