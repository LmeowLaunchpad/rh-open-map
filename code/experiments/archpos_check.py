# Independent check of the ARCHPOS LP optimum: recompute A(lambda) on the FREQUENCY side,
# A = 2 int lambda cosh(u/2) du + (1/2pi) int M_inf(t) lambdahat(t) dt, for the LP minimiser; and refine h.
import numpy as np
from math import pi
from scipy.special import psi
import importlib.util, sys
spec = importlib.util.spec_from_file_location("lp", "archpos_lp2.py"); sys.argv = ["x"]; lp = importlib.util.module_from_spec(spec); spec.loader.exec_module(lp)
for X, h in [(1.0, 0.02), (1.0, 0.01), (2.0, 0.01), (0.693, 0.01)]:
    res, u = lp.solve(X, h=h)
    lam = res.x
    pole = 4*np.trapezoid(lam*np.cosh(u/2), u)          # piecewise linear: trapezoid exact for linear*... approx
    t = np.linspace(0, 4000, 400001)
    P = lam[0] + 2*np.cos(np.outer(t[:, None].ravel(), u[1:])) @ lam[1:] if False else None
    # lambdahat(t) = h sinc^2(t h/2) P(t)
    Pt = np.empty_like(t); B = 20000
    for i in range(0, len(t), B):
        tt = t[i:i+B]; Pt[i:i+B] = lam[0] + 2*np.cos(np.outer(tt, u[1:])) @ lam[1:]
    sinc2 = np.where(t == 0, 1.0, (np.sin(t*h/2)/(t*h/2 + (t == 0)))**2)
    lh = h*sinc2*Pt
    M = np.real(psi(0.25 + 0.5j*t)) - np.log(np.pi)
    arch = (1/pi)*np.trapezoid(M*lh, t)                   # (1/2pi) int over R = (1/pi) int_0^inf (even)
    # tail beyond 4000: lh ~ h sinc^2 ~ 4/(h t^2) * avg(P) ; bound with M ~ log t
    tailest = (1/pi)*np.trapezoid((np.log(t[-1]/2/pi)+0)*(4/(h*t[-1]**2))*abs(Pt[-2000:]).mean()*np.ones(2), [t[-1], 1e9]) if False else 0.0
    print(f"X={X} h={h}: LP value {res.fun:+.5f} | frequency-side A = {pole + arch:+.5f} (pole {pole:.4f}, arch {arch:+.4f}); min P on grid {Pt.min():.2e}")
