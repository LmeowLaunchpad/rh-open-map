# Rig, first step, as a rigorous-by-certificate bound valid for EVERY positive pair (pi, nu) with zeta's Gamma_R + pole data:
# for even h >= 0 with hhat >= 0:  2 int h dpi = A(h) - int hhat dnu <= A(h).  If h >= 1 on [0, u1] then pi((0, u1]) <= A(h)/2.
# LP: minimise A(h) over doubly positive piecewise-linear h (support [-X, X]) with h >= 1 on [0, u1].
import numpy as np, sys
from scipy.optimize import linprog
import importlib.util
spec = importlib.util.spec_from_file_location("lp", "archpos_lp2.py"); sys_argv = sys.argv; sys.argv = ["x"]
lp = importlib.util.module_from_spec(spec); spec.loader.exec_module(lp); sys.argv = sys_argv
def solve(X, u1, h=0.01, over=8):
    # reuse the objective construction from archpos_lp2 by calling solve() internals: rebuild cost via a trick
    u = np.arange(0, X + 1e-12, h); n = len(u)
    res0, _ = lp.solve(X, h=h, over=over)          # to get the same cost vector we recompute it below
    # recompute cost exactly as in archpos_lp2.solve
    from math import log, pi
    EG = 0.5772156649015329; psi14 = -EG - pi/2 - 3*log(2)
    uf = np.linspace(0, X, 40*(n - 1) + 1)
    hat = np.maximum(0, 1 - np.abs(uf[None, :] - u[:, None])/h)
    pole = 4*np.trapezoid(hat*np.cosh(uf/2)[None, :], uf, axis=1)
    ufs = uf.copy(); ufs[0] = 1e-12
    a_f = np.where(uf > 0, 2*np.exp(-uf/2)/(1 - np.exp(-2*ufs)), 0.0)
    A0 = np.trapezoid(a_f*(1 - hat[0]), uf)
    coef_sing = -np.trapezoid(hat[1:]*a_f[None, :], uf, axis=1)
    ug = np.linspace(X, X + 80, 80001); tailX = 2*np.trapezoid(np.exp(-ug/2)/(1 - np.exp(-2*ug)), ug)
    cost = pole.copy(); cost[0] += -(log(pi) + EG) + (psi14 + EG) + A0 + tailX; cost[1:] += coef_sing
    t = np.linspace(0, pi/h, over*n + 1)
    Pm = np.cos(np.outer(t, u)); Pm[:, 1:] *= 2
    lb = np.where(u <= u1 + 1e-12, 1.0, 0.0)
    res = linprog(cost, A_ub=-Pm, b_ub=np.zeros(len(t)), bounds=list(zip(lb, [None]*n)), method='highs')
    return res, u
for u1 in [0.3, 0.6]:
    for X in [1.0, 2.0, 3.0]:
        for hh in [0.02, 0.01]:
            res, u = solve(X, u1, h=hh)
            print(f"u1 = {u1}, X = {X}, h = {hh}: min A(h) = {res.fun:.5f}  => pi((0,{u1}]) <= {res.fun/2:.5f} for every positive pair", flush=True)
