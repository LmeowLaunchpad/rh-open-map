# P1. ARCHPOS

**Statement.** For every even, compactly supported lambda with lambda >= 0 and lambdahat >= 0,

    A(lambda) = 2 int lambda(u) cosh(u/2) du + (1/2pi) int M_inf(t) lambdahat(t) dt >= 0,   M_inf(t) = Re psi(1/4 + it/2) - log pi.

**Known.** RH => ARCHPOS [proved]. Holds on the Poitou cone (lambda(u) cosh(bu) positive definite for 0 <= b <= 1/2) [proved,
classical mechanism]. Holds for all supports up to log 11 (from the certified windows via Boas-Kac) [proved, computer-assisted].
Tight: the infimum of A/lambda(0) appears to be 0, approached by lambda whose transform nearly vanishes at the zeta zeros [checked].

**Why it matters.** With P2 (Rig), ARCHPOS => RH. Without P2 it is still a clean necessary condition: a single explicit lambda with
A(lambda) < 0 would disprove RH by one computation (Odlyzko-type explicit formula argument).

**Entry points.** code/experiments/archpos_lp2.py (dual LP with exact positive-definiteness), archpos_check.py (frequency-side check).
An analytic proof must "see" the zeros through the Gamma function alone; known partial: Connes-Consani's archimedean positivity on a
finite-codimension subspace (support log 2).
