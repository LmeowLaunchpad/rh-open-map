# Code

Requirements: Python 3.12, numpy, scipy, mpmath, python-flint 0.9 (Arb). Run from the directory containing the scripts.

## certificates/ (rigorous, Arb ball arithmetic)
- `t2_cert.py`, `wcert3.py`, `cert2.py`, `wlow.py`: the frozen-symbol engine of agent E4 (pointwise envelope of the comb). Certified
  x = 3, 5, 7. `wlow.py` is a double-precision preview.
- `wcert_h.py`: the invisible-harmonic certificate (wcert3.py with the comb replaced by -C_p + positive atoms, Caratheodory-Toeplitz /
  Pisarenko, exact moment solve in Arb). `patch.py` shows the exact diff from wcert3.py (it expects the original as `wcert3_orig.py`).
- `wlow_h.py`: double-precision preview of the invisible-harmonic certificate (`mode = h` or `zhu`).
- `cx.py`, `base_table.py`: the constant C_x vs Zhu's A_x and the frequency cutoffs.

Window 11, two stages (about 12 minutes on a laptop, ~2.4 GB peak):

    NOCERT=1 NSMALL=80 python wcert_h.py 11 1280 800 0 0 384 15 5 200 3
    python cert2.py GH_x11_N1280_T800.pkl tailH_x11_N1280_T800.pkl 1e-7 16

## experiments/ (numerical, not rigorous)
- `conspire.py x Tmax`: pointwise vs window-smoothed prime-phase conspiracies.
- `archpos_lp2.py X ...`: dual LP for ARCHPOS (exact positive-definiteness via the hat-basis cosine polynomial); `archpos_check.py`
  recomputes the optimum on the frequency side.
- `archpos_ac.py`: test for a smooth (absolutely continuous) positive pair.
- `rig_mass.py`: first-step rigidity bounds on the prime mass below u1 for every positive pair.

Larger engines (H1's sieve, E1's prolate/CBS engines, E2's bootstrap LP and slide tests, E3's fake-world builders) exist in the original
run and can be added on request.
