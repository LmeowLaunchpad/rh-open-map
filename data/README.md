# Data

- `E4_frozen_symbol_x7_cert.txt`: certification log of the frozen-symbol engine at x = 7 (agent E4; T0 = 900, N = 1000).
- Invisible-harmonic certificates (re-create with code/certificates, commands in code/README.md):

| Window | T0 | N | bits | quadrature unit | tail eps | alpha' | certified levels, even (lowest 4) | certified levels, odd (lowest 4) |
|---|---|---|---|---|---|---|---|---|
| 7 | 220 | 400 | 256 | 2.9e-53 | 7.4e-68 | 1.314 | 4.58e-11, 5.55e-16, 1.13e-21, 4.41e-28 | 5.56e-9, 2.26e-13, 9.28e-19, 1.10e-24 |
| 11 | 800 | 1280 | 384 | 1.0e-75 | 1.4e-73 | 1.080 | 5.53e-29, 5.82e-35, 1.29e-41, 6.73e-49 | 2.14e-26, 7.64e-32, 2.91e-38, 4.25e-45 |
| 13 | 1700 | 2700 | 448 | 2.6e-81 | 1.8e-59 | 1.110 | (running) | (running) |

Separable constants C_x (Caratheodory-Toeplitz) vs pointwise A_x: x = 5: 1.52 / 2.94; 7: 2.24 / 4.38; 11: 3.77 / 7.08; 13: 4.49 / 8.52.
