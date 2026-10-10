# Click test (E2's ARCHPOS / Rig architecture): does zeta's archimedean+pole distribution S = 2cosh(u/2) + W_inf admit an
# EXACT absolutely-continuous positive pair?  pi = s(u) chi(u) on u > 0 (s = density of S off 0, >= 0 for |u| > ~0.28),
# nu = Rhat/(2 pi) with R = S (1 - chi).  If Rhat(t) >= 0 for all t: positive pair exists (ARCHPOS holds), and it is not zeta's.
import numpy as np
from math import log, pi
EG = 0.5772156649015329
def s_dens(u): return 2*np.cosh(u/2) - np.exp(-u/2)/(1 - np.exp(-2*u))
def chi(u, a, b):   # smooth step 0 (u<=a) -> 1 (u>=b)
    x = np.clip((u - a)/(b - a), 0, 1)
    return np.where(x <= 0, 0.0, np.where(x >= 1, 1.0, 1/(1 + np.exp(1/np.maximum(1 - x, 1e-300) - 1/np.maximum(x, 1e-300)))))
def Rhat(t, a, b, n=200001):
    u = np.linspace(1e-7, b, n); w = np.full(n, (b - 1e-7)/(n - 1)); w[0] = w[-1] = w[0]/2
    one_m = 1 - chi(u, a, b)
    c = np.cos(np.outer(t, u))
    pole = 4*(c*(np.cosh(u/2)*one_m)) @ w
    arch_in = ((np.exp(-2*u)[None, :] - (np.exp(-u/2)*one_m)[None, :]*c)/(1 - np.exp(-2*u))[None, :]) @ w
    tail = -0.5*np.log(1 - np.exp(-2*b))          # int_b^inf e^{-2u}/(1-e^{-2u}) du  (chi = 1 there)
    return pole - (log(pi) + EG) + 2*(arch_in + tail)
print("s(u) sign change near", [round(x, 3) for x in np.linspace(0.2, 0.35, 16) if s_dens(x) > 0][:1])
t = np.concatenate([np.linspace(0, 20, 801), np.linspace(20, 400, 1521)])
best = None
for a in [0.28, 0.32, 0.36, 0.40, 0.45, 0.50, 0.55, 0.60]:
    for b in [a + d for d in (0.05, 0.1, 0.2, 0.3, 0.5, 0.8)]:
        r = Rhat(t, a, b, n=40001)
        m = r.min()
        if best is None or m > best[0]: best = (m, a, b, t[np.argmin(r)])
print("best min Rhat over t in [0,400]: %.4f at a=%.2f b=%.2f (argmin t=%.2f)" % best)
m, a, b, _ = best
r = Rhat(t, a, b, n=200001)
print("refined: min %.5f at t=%.3f; Rhat(0)=%.4f, Rhat(100)=%.4f, Rhat(400)=%.4f (expect ~log t growth)" % (r.min(), t[np.argmin(r)], r[0], r[np.searchsorted(t,100)], r[-1]))
