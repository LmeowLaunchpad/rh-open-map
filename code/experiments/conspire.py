# Do the prime phases "conspire" pointwise along the line, while the certified window form stays positive?
# D_x(t) = 2 sum_{n<x} Lambda(n) n^{-1/2} cos(t log n);  M_inf(t) = Re psi(1/4 + it/2) - log pi.
# Fejer-smoothed symbol (modulated f = indicator on [-L/2, L/2]): Dfej(t) = 2 sum c_n (1 - log n / L) cos(t log n).
import numpy as np, sys
from scipy.special import psi
x = float(sys.argv[1]); Tmax = float(sys.argv[2]); L = np.log(x)
def lam(n):
    for p in range(2, n+1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return np.log(p) if m == 1 else 0.0
ns = [n for n in range(2, int(x)+1) if n < x and lam(n) > 0]
c = np.array([lam(n)/np.sqrt(n) for n in ns]); ln = np.log(ns)
A = 2*c.sum()
best = (-1e9, 0); bestf = (-1e9, 0); cnt = 0; first = None; widths = []
step = 0.004; t0 = 20.0
while t0 < Tmax:
    t = np.arange(t0, min(t0 + 2e5*step, Tmax), step)
    M = np.real(psi(0.25 + 0.5j*t)) - np.log(np.pi)
    C = np.cos(np.outer(t, ln))
    D = 2*(C @ c); Df = 2*(C @ (c*(1 - ln/L)))
    e = D - M; ef = Df - M
    i = np.argmax(e)
    if e[i] > best[0]: best = (e[i], t[i])
    j = np.argmax(ef)
    if ef[j] > bestf[0]: bestf = (ef[j], t[j])
    pos = e > 0
    if pos.any():
        cnt += np.count_nonzero(np.diff(pos.astype(int)) == 1) + (1 if pos[0] else 0)
        if first is None: first = t[np.argmax(pos)]
        # widths of exceedance intervals
        d = np.diff(np.concatenate([[0], pos.astype(int), [0]]))
        s, f_ = np.where(d == 1)[0], np.where(d == -1)[0]
        widths += list((f_ - s)*step)
    t0 = t[-1] + step
print(f"x = {x}: A_x = {A:.3f}, scanned t in [20, {Tmax:g}]")
print(f"  pointwise: max (D_x - M_inf) = {best[0]:+.3f} at t = {best[1]:.3f}; exceedance intervals: {cnt}; first at t = {first}")
if widths: print(f"  exceedance widths: max {max(widths):.3f}, median {np.median(widths):.3f}  (1/L = {1/L:.3f}, window resolution ~2/L = {2/L:.3f})")
print(f"  Fejer-smoothed: max (Dfej - M_inf) = {bestf[0]:+.3f} at t = {bestf[1]:.3f}")
