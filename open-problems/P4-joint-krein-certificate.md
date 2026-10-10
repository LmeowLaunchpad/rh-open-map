# P4. An explicit joint (Krein-optimal) invisible-harmonic certificate

The separable certificate (results/01) treats each prime's comb on its own circle, with constant C_x ~ 2 sqrt(x). The best possible
constant using invisible harmonics of all primes jointly is the compressed atom norm Lambda_x ~ 1.1 sqrt(x) (Krein duality), but no
explicit extension realising it is known. Build one (e.g. an LP/SDP over mixed harmonics sum k_p log p > log x, made exact as in
wcert_h.py's Pisarenko construction).
**Payoff.** The frequency cutoff drops from ~2 pi e^{C_x} to ~2 pi e^{Lambda_x}: roughly a factor 4 at x ~ 13, putting windows 17-20
within reach of a laptop, and improving the reach of P3.
