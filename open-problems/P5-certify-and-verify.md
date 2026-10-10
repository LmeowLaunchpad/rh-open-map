# P5. Verify the certificates; certify further windows

- Independently re-run x = 7 and x = 11 with different parameters (cutoff, basis size, precision, quadrature panels), and audit
  wcert_h.py (the invisible-harmonic patch of the frozen-symbol engine): the invisibility condition (K_p + 1) log p > log x, the
  normalisation of the atoms against the (1/pi) int_0^T0 convention, the exact moment solve and positivity of the atom weights, the
  added tail term.
- Push to x = 13 (running at the time of writing: T0 = 1700, N = 2700, 448 bits) and beyond with P4.
- Cross-check against Zhu (support <= 1.6) and Groskin (Galerkin compressions to x = 100).
