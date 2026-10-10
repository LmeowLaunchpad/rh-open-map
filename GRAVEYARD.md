# Graveyard: routes that fail, and why

One line each. "Dead" means the route cannot work as stated; the reason is the useful part.

## Window certificates
- **Pointwise envelope of the prime comb** (frozen symbol with 2 sum c_n): cost ~e^{4 sqrt x} by Zhu's Thm 1.4; x = 8 already needs T0 ~ 3600.
- **Compressed-atom norm Lambda_x used in u-space** (E4): atoms are translations, non-smoothing, so the band cannot be truncated. Fixed by moving the same bound to the Fourier side (invisible harmonics).
- **Any finite-window certificate reaching RH**: impossible; long-range births force frequencies ~e^{sqrt x} at window x.

## Base lemma (gluing verified zeros with the prime side), Sub-lemma S
- **Disc interpolation at zeros** (global bound, no pollution): needs zero density >= ~3.25x Nyquist; reaches x ~ 13 only.
- **Groechenig max-gap on a localised F**: window suppression is sub-exponential, Slepian leakage forces bandwidth ~(1 + 2e) sigma; reaches tiny x only.
- **Explicit-formula localisation** (test w|F|^2): exact but re-derives the prime side; circular.
- **Gram-lattice Poisson** (exact for a lattice): zeros deviate 0.5-3 spacings; perturbation needs ~10x Nyquist.
- **Lagrange interpolation through zeta** ((z + i/2) zeta(1/2 + iz) e^{i tau z}): trivial zeros enter with weight ~cosh(2 pi e^{sigma + tau}).
- **Near-null split + max-gap sampling on the complement**: needs an explicit Beurling theorem in the 1-5x Nyquist zone (open, P3).
- **|F|-based estimates** (Poisson majorants, Bernstein, Remez/Turan, Logvinenko-Sereda): blind to the cancellation that makes leakage small.

## The growth step / induction
- **Non-shrinking margin gamma0^2 ~ 0.41**: superseded at x = 25 (gamma^2 = 0.51), eaten by tight zero pairs at the horizon.
- **Arrival dip of the angle**: artefact of the basis (withdrawn by H1).
- **Edge-layer explanation of the angle** (E3 prediction 0.03): refuted (measured 0.13-0.36).
- **Look-ahead induction** (H1): no look-ahead above the horizon 4 pi x for the angle or positivity.
- **Plain time/band separation for a uniform angle**: by Landau-Pollak the cosine of the angle tends to 1 as c = 2 pi x grows.
- **Spectral split by the archimedean eigenbasis** (E2): the atoms couple the blocks with O(1) non-smoothing elements.
- **V_K splitting** (E4): fails its click test at x = 3 (rho = 1.39 / 2.30).
- **Convex combinations of prime-side and zero-side bounds**: fail on ladder vectors (e^{-4 pi x} vs O(1)).
- **Hadamard rate K as an early-warning detector**: sharp at births but no early warning (DH: ~1 until just before).

## Phases and positivity mechanisms
- **"The prime phases never align" pointwise**: false; thousands of exceedances in zero gaps. Only the window-resolution version holds, and it is RH.
- **Damped windows e^{-eta|u|} W as a global flow**: positivity <=> eta >= sup Re(rho) - 1/2; just the zero-free half-plane.
- **Splitting the archimedean form among the primes** (each bracket w_n (P + A) - 2 c_n R_n PSD): impossible once P + A has a negative direction and there are >= 2 primes.
- **Delsarte / magic-function dual of existence**: its value equals the maximin margin m(x), exponentially small (2.7e-9 at x = 5); no robust margin.
- **Reflection-group frustration alone**: gives Lambda_x ~ 0.3 of the pointwise worst case; not enough without the archimedean term.
- **Positivity of the weights alone**: zeta(s-a) zeta(s+a) has positive weights and all zeros off the line.
- **Torus-uniform positivity**: fails already for S = {inf, 2} away from the exact phase.

## Existence / rigidity
- **ARCHPOS with slack**: false numerically (tight); slack would also refute Rig.
- **Smooth (absolutely continuous) positive pair, family pi = s chi**: zero side negative (min -0.054).
- **Existence (L4) as an easier route**: cannot be easier than the direct target (E2).

## Earlier phase (common-signal class, results/07)
- ~50 attempts at the 7/8 wall, each with a named failing step; outside-in bridges A1-A13; integrality/divisibility engines (closed at
  Mertens scale); Landau-Deligne positivity class (void over Q); Jacobi-sum Moebius source (fails at self-reference); trace functions on
  A^1 ((M) and (S) incompatible); Steinberg/Kottwitz sign routes (die at (F)); t-averaged Jacobi probe (coherence cancelled by shifts
  sharing a prime with c).
