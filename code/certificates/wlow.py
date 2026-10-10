# E4 (T2) click test: lower-symbol form W_low (symbol M_eff on |t|<=T0, constant alpha beyond) in an orthonormal Legendre basis.
# W(f) = (1/2pi) int |fhat|^2 M_eff dt + P(f),  M_eff(t) = Re psi(1/4+it/2) - log pi - 2 sum_{n<=x} Lambda(n)/sqrt(n) cos(t log n)
# Usage: python wlow.py x N T0 [alpha|auto] [par]
import sys, numpy as np
from scipy.special import psi, spherical_jn
from numpy.polynomial.legendre import leggauss, legval
x = float(sys.argv[1]); N = int(sys.argv[2]); T0 = float(sys.argv[3])
alpha_arg = sys.argv[4] if len(sys.argv) > 4 else 'auto'
L = np.log(x); a = L/2
def lam(n):
    for p in range(2, n+1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return np.log(p) if m == 1 else 0.0
    return 0.0
atoms = [(np.log(n), lam(n)/np.sqrt(n)) for n in range(2, int(np.floor(x))+1) if lam(n) > 0 and np.log(n) < L]
def Meff(t):
    v = np.real(psi(0.25 + 0.5j*t)) - np.log(np.pi)
    for ln, c in atoms: v = v - 2*c*np.cos(t*ln)
    return v
# quadrature in t on [0, T0]: panels of length h
h = 0.5; P = int(np.ceil(T0/h)); gx, gw = leggauss(24)
edges = np.linspace(0, T0, P+1)
t = ((edges[:-1, None] + edges[1:, None])/2 + (edges[1:, None]-edges[:-1, None])/2*gx[None, :]).ravel()
w = ((edges[1:, None]-edges[:-1, None])/2*gw[None, :]).ravel()
m = Meff(t)
if alpha_arg == 'auto':
    tt = np.linspace(T0, 50*T0, 200001); alpha = Meff(tt).min()
elif alpha_arg == 'rig':
    alpha = np.real(psi(0.25+0.5j*T0)) - np.log(np.pi) - 2*sum(c for _, c in atoms)
else:
    alpha = float(alpha_arg)
# basis e_n(u) = sqrt((2n+1)/(2a)) P_n(u/a); ehat_n(t) = sqrt((2n+1) a/2) * 2 i^n j_n(t a)
J = np.array([spherical_jn(n, t*a) for n in range(N+1)])          # (N+1, nt)
cn = np.sqrt((2*np.arange(N+1)+1)*a/2)*2
# M_jk = (1/2pi) int_{-T0}^{T0} (m - alpha) conj(ehat_j) ehat_k dt = (1/pi) int_0^T0 (m-alpha) c_j c_k i^{k-j} j_j j_k dt (same parity)
F = J*cn[:, None]
G = (F*((m - alpha)*w)[None, :]) @ F.T/np.pi
n = np.arange(N+1); ph = np.real(1j**((n[None, :]-n[:, None]) % 4)); same = ((n[None, :]-n[:, None]) % 2 == 0)
G = np.where(same, G*ph, 0.0)
W = G + alpha*np.eye(N+1)
# pole: P(f) = 2 (int f e^{u/2})(int f e^{-u/2}); A_j = int e_j e^{u/2} du via Gauss in u
from scipy.special import spherical_in
A = np.array([np.sqrt((2*k+1)/(2*a))*2*a*spherical_in(k, a/2) for k in range(N+1)]); B = A*(-1.0)**np.arange(N+1)
W = W + np.outer(A, B) + np.outer(B, A)
even = n % 2 == 0; odd = ~even
for name, sel in (('even', even), ('odd', odd)):
    ev = np.linalg.eigvalsh(W[np.ix_(sel, sel)])
    print('x=%g N=%d T0=%g alpha=%.4f %s: lowest %s' % (x, N, T0, alpha, name, ' '.join('%.4e' % e for e in ev[:4])), flush=True)
