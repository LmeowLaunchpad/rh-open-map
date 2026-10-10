# W_low with the invisible-harmonic (Caratheodory-Toeplitz) certificate: per prime p, the comb part
# F_p(theta) = -2 sum_{k<=K_p} c_{p^k} cos(k theta) is replaced by -C_p + mu_p (mu_p >= 0 atomic, Fourier coeffs matching
# C_p + F_p up to order K_p); the difference has only harmonics k > K_p, i.e. frequencies k log p > log x: invisible to the window.
# Usage: python wlow_h.py x N T0 [mode]   mode = 'h' (new) or 'zhu' (E4/Zhu pointwise freeze)   (double precision preview)
import sys, numpy as np
from scipy.special import psi, spherical_jn, spherical_in
from numpy.polynomial.legendre import leggauss
x = float(sys.argv[1]); N = int(sys.argv[2]); T0 = float(sys.argv[3]); mode = sys.argv[4] if len(sys.argv) > 4 else 'h'
L = np.log(x); a = L/2
def isprime(p): return p > 1 and all(p % q for q in range(2, int(p**.5)+1))
Minf = lambda t: np.real(psi(0.25 + 0.5j*t)) - np.log(np.pi)
atoms_t, atoms_w, C, A = [], [], 0.0, 0.0
combs = []
for p in range(2, int(x)+1):
    if not isprime(p): continue
    K = int(np.floor(np.log(x)/np.log(p) - 1e-12))
    if K < 1: continue
    c = np.array([np.log(p)*p**(-k/2) for k in range(1, K+1)]); A += 2*c.sum()
    combs.append((np.log(p), c))
    T = np.zeros((K+1, K+1))
    for i in range(K+1):
        for j in range(K+1):
            if i != j: T[i, j] = -c[abs(i-j)-1]
    ev, V = np.linalg.eigh(T); lam, v = ev[0], V[:, 0]; Cp = -lam; C += Cp
    # extremal (Pisarenko) measure: support = roots of q(z) = sum v_j z^j (on |z| = 1), weights from moments k = 0..K
    r = np.roots(v[::-1]); th = np.angle(r)
    Mk = np.array([[np.exp(-1j*k*tj) for tj in th] for k in range(K+1)])
    rhs = np.array([Cp] + list(-c), dtype=complex)            # hat(Cp + F_p)_k : k=0 -> Cp, k>=1 -> -c_k
    wts = np.linalg.lstsq(Mk, rhs, rcond=None)[0]
    res = np.abs(Mk @ wts - rhs).max()
    print(f"p={p} K={K} pointwise={2*c.sum():.4f} toeplitz={Cp:.4f} |z|-1={np.abs(np.abs(r)-1).max():.1e} "
          f"weights={np.round(wts.real,4)} imag={np.abs(wts.imag).max():.1e} moment_res={res:.1e}")
    lp = np.log(p)
    for tj, wj in zip(th, wts.real):
        for m in range(-int(T0*lp/(2*np.pi))-2, int(T0*lp/(2*np.pi))+3):
            t = (tj + 2*np.pi*m)/lp
            if 0 < t <= T0: atoms_t.append(t); atoms_w.append(2*np.pi*wj/lp)
            elif t == 0: atoms_t.append(0.0); atoms_w.append(np.pi*wj/lp)
atoms_t = np.array(atoms_t); atoms_w = np.array(atoms_w)
h = 0.5; P = int(np.ceil(T0/h)); gx, gw = leggauss(24)
edges = np.linspace(0, T0, P+1)
t = ((edges[:-1, None]+edges[1:, None])/2 + (edges[1:, None]-edges[:-1, None])/2*gx[None, :]).ravel()
w = ((edges[1:, None]-edges[:-1, None])/2*gw[None, :]).ravel()
if mode == 'zhu':
    m = Minf(t) - sum(2*(cc*np.cos(np.outer(t, lp*np.arange(1, len(cc)+1)))).sum(1) for lp, cc in combs)
    alpha = Minf(T0) - A
else:
    m = Minf(t) - C
    alpha = Minf(T0) - C
cn = np.sqrt((2*np.arange(N+1)+1)*a/2)*2
F = np.array([spherical_jn(n, t*a) for n in range(N+1)])*cn[:, None]
G = (F*((m - alpha)*w)[None, :]) @ F.T/np.pi
if mode != 'zhu' and len(atoms_t):
    Fa = np.array([spherical_jn(n, atoms_t*a) for n in range(N+1)])*cn[:, None]
    G += (Fa*atoms_w[None, :]) @ Fa.T/np.pi
n = np.arange(N+1); ph = np.real(1j**((n[None, :]-n[:, None]) % 4)); same = ((n[None, :]-n[:, None]) % 2 == 0)
G = np.where(same, G*ph, 0.0)
W = G + alpha*np.eye(N+1)
Av = np.array([np.sqrt((2*k+1)/(2*a))*2*a*spherical_in(k, a/2) for k in range(N+1)]); Bv = Av*(-1.0)**np.arange(N+1)
W = W + np.outer(Av, Bv) + np.outer(Bv, Av)
print(f"x={x} N={N} T0={T0} mode={mode} A_x={A:.4f} C_x={C:.4f} alpha={alpha:.4f} atoms_in_band={len(atoms_t)}")
for name, sel in (('even', n % 2 == 0), ('odd', n % 2 == 1)):
    ev = np.linalg.eigvalsh(W[np.ix_(sel, sel)])
    print(f"  {name}: lowest " + ' '.join('%.4e' % e for e in ev[:4]))
