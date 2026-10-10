# 03. Arrival of a prime power: one universal object

## Arrival Lemma and reflection identity [proved, H1]
When the window passes a prime power n, the new atom enters only as a boundary coupling on a thin edge strip
J = [log n - L/2, L/2]: its form is B_n = 2 R_b, twice the reflection about b = log(n)/2 restricted to J. Odd-about-b functions on the
strip are lifted, even ones lowered.

## First step anatomy, S = {inf, 2} below log 3 [checked]
The archimedean form alone has exactly one negative level on windows x in (2.271, 3]; it is odd about log sqrt 2 and is lifted by the
atom at 2; the even second level caps the weight from above. Converged margin of the first semilocal step: 5.54e-8 (gap to the second
level 1.46e-3). Positivity holds only at the exact phase (E3): shifting the phase of the prime 2 on its torus breaks it
(lowest eigenvalue -0.27 at the opposite phase).

## The universal family (hypothesis of R. Malcolm, tested by H1) [checked]
Count of negative directions each arrival must repair (higher precision, dps 30-66):

| n | 2 | 3 | 4 | 5 | 7 | 8 | 9 | 11 | 13 | 16 | 19 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| directions | 1 | 2 | 2 | 4 | 3 | 3 | 5 | 5 | >= 6 | 4 (consistent) | 10 resolvable |

- The count is a monotone function of one layer coordinate z (one new direction per ~7.5 in z), prime powers included.
- The j-th direction is the j-th odd strip mode with 2j - 1 nodes at universal positions; shapes overlap 0.84-1.00 across primes.
- Directions are inherited step to step (principal cosines 0.93-1.00; 0.52 across the long gap 5 -> 7): about one genuinely new
  direction per step.
- Out of sample: n = 16 (prediction 4, observed 4 to dps 64, a fifth possible below precision); n = 19 at z = 81 (beyond the fitted
  range): ten resolvable modes of the same Sturm-ordered family.
- Without any primes (Gamma_R + pole only) the number of negative levels is only 1, 2, 3 for x up to 60 (thresholds 2.27, 8.8, 24.7).
- E1's block model with no fitting reproduces the counts exactly at n = 2, 4, 7, 8 and within one at 3, 5, 9.

## What it means
Arrivals are completely understood as structure, and **RH-neutral**. A single lemma in the layer coordinate (the boundary-layer law,
results/04) covers every arrival.
