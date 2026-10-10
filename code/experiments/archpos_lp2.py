# ARCHPOS dual LP, exact-PD version: lambda piecewise linear on grid u_j = j h (hat basis), even, supp [-X, X].
# lambdahat(t) = h sinc^2(t h / 2) * P(t),  P(t) = lam_0 + 2 sum_{j>=1} lam_j cos(t u_j)  (period 2pi/h, even):
# PD <=> P >= 0 on [0, pi/h]  (enforced on a dense grid).  Objective A(lambda) = <S, lambda>, normalisation lambda(0) = 1.
import numpy as np, sys
from scipy.optimize import linprog
from math import log, pi
EG = 0.5772156649015329
psi14 = -EG - pi/2 - 3*log(2)
def solve(X, h=0.02, over=8):
    u = np.arange(0, X + 1e-12, h); n = len(u)
    # integrals of hat functions against a weight wfun on [0, X]: use fine quadrature
    uf = np.linspace(0, X, 40*(n - 1) + 1)
    hat = np.maximum(0, 1 - np.abs(uf[None, :] - u[:, None])/h)            # (n, nf)
    def integ(wf):  return np.trapezoid(hat*wf[None, :], uf, axis=1)
    pole = 4*integ(np.cosh(uf/2))
    ufs = uf.copy(); ufs[0] = 1e-12
    a_f = np.where(uf > 0, 2*np.exp(-uf/2)/(1 - np.exp(-2*ufs)), 0.0)
    # 2 int_0^X e^{-u/2}(lam0 - lam(u))/(1-e^{-2u}) du : lam0 - lam(u) = sum_j lam_j (delta_{j0} - hat_j(u))
    # coefficient for lam_j : -int a hat_j  (j >= 0), plus lam_0 * int_0^X a  [the hat_0 part cancels the 1/u singularity]
    A0 = np.trapezoid(a_f*(1 - hat[0]), uf)                               # lam0*int a (1 - hat_0): finite
    coef_sing = -np.trapezoid(hat[1:]*a_f[None, :], uf, axis=1)            # j >= 1
    ug = np.linspace(X, X + 80, 80001); tailX = 2*np.trapezoid(np.exp(-ug/2)/(1 - np.exp(-2*ug)), ug)
    cost = pole.copy()
    cost[0] += -(log(pi) + EG) + (psi14 + EG) + A0 + tailX
    cost[1:] += coef_sing
    t = np.linspace(0, pi/h, over*n + 1)
    Pm = np.cos(np.outer(t, u)); Pm[:, 1:] *= 2
    res = linprog(cost, A_ub=-Pm, b_ub=np.zeros(len(t)), A_eq=np.eye(1, n), b_eq=[1.0], bounds=(0, None), method='highs')
    return res, u
for X in [float(s) for s in sys.argv[1:]]:
    res, u = solve(X)
    if res.status != 0: print(X, res.message); continue
    lam = res.x
    print(f"X = {X:.3f} (x = {np.exp(X):6.2f}): min A/lambda(0) = {res.fun:+.5f}; int lambda = {2*np.trapezoid(lam, u):.4f}; max lambda = {lam.max():.3f}", flush=True)
