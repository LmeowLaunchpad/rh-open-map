# P6. Is the per-window margin uniform?

Per window, RH-content = (bulk positive) and gamma(x) < 1, gamma the angle between the prolate near-null span and the bulk. In
prolate-whitened coordinates this is an O(1) matrix condition M(x) > 0. Data: gamma^2 = 0.13 / 0.33 / 0.41 / 0.39 / 0.51 at
x = 5 / 10 / 15 / 20 / 25, driven by tight zero pairs just inside the horizon 4 pi x (pair-locking). Questions:
1. Does lambda_min(M(x)) stay above a fixed m0 > 0? (Selberg's CLT and Conrey-Li suggest a strictly uniform bound may fail at
   astronomically large x.)
2. Is there a global identity (built from the Euler product, the functional equation and unitary local factors jointly) that forces
   gamma(x) < 1? This is the induction step; RH-strength.
Entry points: H1's whitened block (results/05), E1's chunked CBS engine, the certified matrices for x = 3, 5, 7, 11, 13.
