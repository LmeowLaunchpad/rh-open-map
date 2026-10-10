# E4 (T2/T3) rigorous certificate engine, python-flint (Arb ball arithmetic). 2026-10-09.
#
# Claim certified: for x = e^L and every f in L^2[-L/2, L/2] (real or complex), the window Weil form
#     W(f) = (1/2pi) int |fhat(t)|^2 M_eff(t) dt + P(f),
#     M_eff(t) = Re psi(1/4 + it/2) - log pi - 2 sum_{n <= x, Lambda(n) > 0} Lambda(n)/sqrt(n) cos(t log n),
#     P(f) = 2 (int f e^{u/2})(int f e^{-u/2}),
# is >= 0 (this is Weil's quadratic functional on test functions g = f * f~ supported in [-L, L]; atoms are all prime powers < x).
#
# Proof skeleton (each step rigorous):
#  (1) M_inf(t) = Re psi(1/4+it/2) - log pi is increasing in |t| (termwise in psi(w) = -gamma + sum (1/(n+1) - 1/(n+w))).
#      Hence for |t| >= T0: M_eff(t) >= alpha := M_inf(T0) - 2 sum c_n  (lower endpoint of an Arb ball).
#  (2) W >= W_low := (1/2pi) int |fhat|^2 m + P with m = M_eff on |t| <= T0 and alpha beyond. W_low = alpha I + K + P,
#      K the self-adjoint operator with kernel (1/2pi) int_{-T0}^{T0} (M_eff - alpha)(t) e^{it(u-v)} dt.
#  (3) Galerkin in the orthonormal Legendre basis e_n(u) = sqrt((2n+1)/(2a)) P_n(u/a), a = L/2, n <= N:
#      ehat_n(t) = c_n i^n j_n(t a), c_n = sqrt(2(2n+1)a);  G_jk = alpha delta_jk + (1/pi) (-1)^{(k-j)/2} c_j c_k
#      int_0^T0 (M_eff - alpha) j_j(ta) j_k(ta) dt (same parity; 0 otherwise) + P_jk, P_jk = A_j A_k ((-1)^j + (-1)^k),
#      A_n = int e_n e^{u/2} = sqrt((2n+1)/(2a)) 2a i_n(a/2).
#      Integrals: composite Gauss-Legendre (Arb nodes/weights) + Trefethen's Bernstein-ellipse error bound with the
#      integrand bounded on each ellipse by Arb evaluation over a ball covering of the ellipse boundary and |j_n(z)| <= e^{|Im z|}.
#  (4) Tail: eps >= ||Pi_perp (K + P)|| via |j_n(x)| <= x^n/(2n+1)!!, i_n(z) <= z^n e^z/(2n+1)!! (x, z >= 0).
#      Then W_low(v + w) >= lam ||v||^2 - 2 eps ||v|| ||w|| + (alpha - eps) ||w||^2 >= 0 if lam (alpha - eps) >= eps^2.
#  (5) lam: G_sector - lam I is certified positive definite by an LDL^T factorisation in ball arithmetic (all pivots > 0).
# Usage: python t2_cert.py x N T0 lam_even lam_odd [prec]
import sys, time
import flint
from flint import arb, acb, arb_mat

x_s, N, T0_s = sys.argv[1], int(sys.argv[2]), sys.argv[3]
lam_even, lam_odd = arb(sys.argv[4]), arb(sys.argv[5])
prec = int(sys.argv[6]) if len(sys.argv) > 6 else 256
flint.ctx.prec = prec
t_start = time.time()
X = arb(x_s); L = X.log(); a = L/2; T0 = arb(T0_s)
pi = arb.pi(); logpi = pi.log()
def lam_vm(n):
    for p in range(2, n+1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return p if m == 1 else 0
    return 0
xf = float(x_s)
atoms = []
for n in range(2, int(xf) + 1):
    p = lam_vm(n)
    if p and n < xf:
        atoms.append((arb(n).log(), arb(p).log()/arb(n).sqrt()))
csum = sum((c for _, c in atoms), arb(0))
def Minf_t(t):          # real t (arb)
    return acb(arb(1)/4, t/2).digamma().real - logpi
def Meff_z(z):          # analytic continuation, complex ball z
    w1 = acb(arb(1)/4) + acb(0, 1)*z/2
    w2 = acb(arb(1)/4) - acb(0, 1)*z/2
    v = (w1.digamma() + w2.digamma())/2 - acb(logpi)
    for ln, c in atoms:
        v = v - 2*c*(z*ln).cos()
    return v
alpha_ball = Minf_t(T0) - 2*csum
alpha = arb(alpha_ball.lower())            # exact number <= M_eff(t) for all |t| >= T0
print('x = %s  L = %s  N = %d  T0 = %s  atoms %d  alpha = %s' % (x_s, L.str(10), N, T0_s, len(atoms), alpha.str(12)), flush=True)

# ---- quadrature setup
hh = arb('0.25'); npan = int(round(float(T0)/0.5)); ng = 30
rho = arb('3.48')
gl = [arb.legendre_p_root(ng, k, weight=True) for k in range(ng)]
# Bernstein-ellipse bound per panel
K_arc = 96
two_a = 2*a
def ellipse_bound(cen):
    best = arb(0)
    semi = hh*(rho + 1/rho)/2
    r_arc = semi*pi/K_arc*arb('1.02')
    for k in range(K_arc):
        th = 2*pi*(arb(k) + arb(1)/2)/K_arc
        e = acb(th.cos(), th.sin())
        zc = acb(cen) + acb(hh)*(acb(rho)*e + 1/(acb(rho)*e))/2
        z = acb(arb(zc.real.mid(), r_arc), arb(zc.imag.mid(), r_arc))
        val = (Meff_z(z) - acb(alpha)).abs_upper() * (two_a*abs(z.imag).upper()).exp()
        ub = arb(val.upper())
        if ub > best: best = ub
    return best
t1 = time.time()
Msum = arb(0)
for p in range(npan):
    cen = hh*(2*p + 1)
    Msum += hh*ellipse_bound(cen)
qfac = arb(64)/15/(rho*rho - 1)*rho**(-2*ng)
quad_err_unit = qfac*Msum/pi           # error of entry (j,k) <= c_j c_k * quad_err_unit
print('ellipse bounds done in %.1f s; sum_p hh*M_p = %s; quad error unit = %s' % (time.time()-t1, Msum.str(6), quad_err_unit.str(6)), flush=True)

# ---- Galerkin matrix
cn = [(2*(2*n+1)*a).sqrt() for n in range(N+1)]
G = arb_mat(N+1, N+1)
t1 = time.time()
chunk = 20
for p0 in range(0, npan, chunk):
    nodes = []
    for p in range(p0, min(npan, p0+chunk)):
        cen = hh*(2*p + 1)
        for (s, wgt) in gl:
            nodes.append((cen + hh*s, hh*wgt))
    F = arb_mat(N+1, len(nodes)); FW = arb_mat(len(nodes), N+1)
    for i, (t, wgt) in enumerate(nodes):
        hval = Minf_t(t) - alpha
        for ln, c in atoms:
            hval = hval - 2*c*(t*ln).cos()
        ta = t*a
        pref = (pi/(2*ta)).sqrt()
        wh = wgt*hval/pi
        for n in range(N+1):
            jn = pref*ta.bessel_j(arb(n) + arb(1)/2)
            if jn.rad() > arb('1e-30'):
                flint.ctx.prec = 2*prec
                jn = (pi/(2*ta)).sqrt()*ta.bessel_j(arb(n) + arb(1)/2)
                flint.ctx.prec = prec
            v = cn[n]*jn
            F[n, i] = v
            FW[i, n] = v*wh
    G = G + F*FW
    print('  panels %d-%d done (%.1f s)' % (p0, min(npan, p0+chunk)-1, time.time()-t1), flush=True)
# signs, parity, alpha, quadrature radius, pole
A = [((2*n+1)/(2*a)).sqrt()*2*a*((pi/a).sqrt()*(a/2).bessel_i(arb(n) + arb(1)/2)) for n in range(N+1)]
# note: i_n(z) = sqrt(pi/(2z)) I_{n+1/2}(z); with z = a/2, sqrt(pi/(2z)) = sqrt(pi/a)
Gs = {}
for par in (0, 1):
    idx = [n for n in range(N+1) if n % 2 == par]
    M = arb_mat(len(idx), len(idx))
    for r, j in enumerate(idx):
        for c, k in enumerate(idx):
            sgn = 1 if ((k - j)//2) % 2 == 0 else -1
            v = sgn*G[j, k]
            v = v + arb(0, (cn[j]*cn[k]*quad_err_unit).upper())
            if j == k: v = v + alpha
            v = v + A[j]*A[k]*((-1)**j + (-1)**k)
            M[r, c] = v
    Gs[par] = M
print('matrix assembled (%.1f s total)' % (time.time() - t_start), flush=True)

# ---- tail bound eps >= ||Pi_perp (K + P)||
# ||Pi_perp e^{it.}||^2 = sum_{n>N} c_n^2 j_n(ta)^2 <= sum_{n>N} 2(2n+1)a (T0 a)^{2n}/((2n+1)!!)^2  (|t| <= T0)
def dfact(m):
    r = arb(1)
    for k in range(m, 0, -2): r *= k
    return r
xa = T0*a
tau2 = arb(0); n = N + 1
while True:
    term = 2*(2*n+1)*a*xa**(2*n)/dfact(2*n+1)**2
    tau2 += term
    if n > N + 5 and term < tau2*arb('1e-30'): break
    n += 1
tau2 = tau2*2       # generous factor for the truncated remainder of the series
hmax = max(abs(Minf_t(arb(0)) - 2*csum - alpha).upper(), abs(Minf_t(T0) + 2*csum - alpha).upper())
epsK = (2*a).sqrt()*(1/(2*pi))*2*T0*arb(hmax)*tau2.sqrt()
# pole tail: ||Pi_perp e^{+-u/2}||^2 = sum_{n>N} (2n+1)/(2a) (2a i_n(a/2))^2,  i_n(z) <= z^n e^z/(2n+1)!!
pt2 = arb(0); n = N + 1
z = a/2
while True:
    term = (2*n+1)/(2*a)*(2*a*z**n*z.exp()/dfact(2*n+1))**2
    pt2 += term
    if n > N + 5 and term < pt2*arb('1e-30'): break
    n += 1
pt2 = pt2*2
normE = (((a).exp() - (-a).exp())).sqrt()        # ||e^{+-u/2}||_{L^2[-a,a]} = sqrt(2 sinh a)
epsP = 2*2*pt2.sqrt()*normE
eps = arb((epsK + epsP).upper())
print('tail: tau_N(T0)^2 <= %s, eps_K <= %s, eps_P <= %s' % (tau2.str(5), epsK.str(5), epsP.str(5)), flush=True)

# ---- certify G - lam I > 0 by ball LDL^T
def ldl_ok(M, lam):
    n = M.nrows()
    Aa = [[M[i, j] for j in range(n)] for i in range(n)]
    for i in range(n): Aa[i][i] = Aa[i][i] - lam
    D = [None]*n; Lm = [[arb(0)]*n for _ in range(n)]
    minpiv = None
    for j in range(n):
        s = Aa[j][j]
        for k in range(j):
            s = s - Lm[j][k]*Lm[j][k]*D[k]
        D[j] = s
        if not (s > 0):
            return False, j, s
        if minpiv is None or s.lower() < minpiv.lower(): minpiv = s
        for i in range(j+1, n):
            t = Aa[i][j]
            for k in range(j):
                t = t - Lm[i][k]*Lm[j][k]*D[k]
            Lm[i][j] = t/s
    return True, n, minpiv
ok_all = True
for par, lam in ((0, lam_even), (1, lam_odd)):
    ok, j, piv = ldl_ok(Gs[par], lam)
    tries = 0
    while not ok and tries < 8:
        print('  %s sector: not certified at lam = %s (failed at pivot %d: %s); lowering lam by 10' % ('even' if par == 0 else 'odd', lam.str(4), j, piv.str(5)), flush=True)
        lam = lam/10; tries += 1
        ok, j, piv = ldl_ok(Gs[par], lam)
    cond = (lam*(alpha - eps) >= eps*eps)
    print('%s sector: G - %s I positive definite: %s (pivots checked %d, min pivot %s); tail condition lam(alpha-eps) >= eps^2: %s' %
          ('even' if par == 0 else 'odd', lam.str(4), ok, j, piv.str(5) if piv is not None else '-', cond), flush=True)
    ok_all = ok_all and ok and cond
# floating check of the midpoint spectrum
import numpy as np
for par in (0, 1):
    M = Gs[par]; n = M.nrows()
    Mf = np.array([[float(M[i, j].mid()) for j in range(n)] for i in range(n)])
    print('  midpoint lowest eigenvalues (%s): %s' % ('even' if par == 0 else 'odd', ' '.join('%.6e' % e for e in np.linalg.eigvalsh(Mf)[:3])))
print('CERTIFIED: W >= 0 on L^2[-L/2, L/2] for x = %s' % x_s if ok_all else 'NOT certified', '(%.1f s)' % (time.time()-t_start), flush=True)
