# E4 rigorous window-positivity certificate, v2 (python-flint / Arb). 2026-10-09.
# Same theorem and proof skeleton as t2_cert.py (see T2_THEOREM.md); changes for larger windows:
#   * quadrature: panels of length 1/2 (30 nodes, rho = 3.48) on [0, ts]; panels of length 2*hb (nb nodes) on [ts, T0] with
#     Bernstein ellipses of semi-minor axis bb (singularities of psi(1/4 +- iz/2) lie on Re z = 0 only);
#   * j_n(ta), n <= N, by downward recurrence j_{n-1} = (2n+1)/x j_n - j_{n+1} from two Arb values (ball arithmetic, rigorous);
#   * tail: ||Pi_perp e^{it.}||^2 = sum_{n>N} c_n^2 j_n(ta)^2 is increasing in |t| on [0, T0] when N >= T0 a (a positive solution
#     of the spherical Bessel ODE has no local maximum where x^2 < n(n+1)); evaluated at t = T0 by direct Arb sums plus a crude
#     x^n/(2n+1)!! remainder.
# Usage: python wcert.py x N T0 lam_even lam_odd [prec] [ts] [hb] [nb] [bb]
import sys, time
import flint
from flint import arb, acb, arb_mat

x_s, N, T0_s = sys.argv[1], int(sys.argv[2]), sys.argv[3]
lam0 = {0: arb(sys.argv[4]), 1: arb(sys.argv[5])}
prec = int(sys.argv[6]) if len(sys.argv) > 6 else 256
ts = float(sys.argv[7]) if len(sys.argv) > 7 else 15.0
hb = float(sys.argv[8]) if len(sys.argv) > 8 else 5.0
nb = int(sys.argv[9]) if len(sys.argv) > 9 else 90
bb = float(sys.argv[10]) if len(sys.argv) > 10 else 3.0
flint.ctx.prec = prec
t_start = time.time()
X = arb(x_s); L = X.log(); a = L/2; T0 = arb(T0_s)
pi = arb.pi(); logpi = pi.log()
def vm_prime(n):
    for p in range(2, n+1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return p if m == 1 else 0
    return 0
xf = float(x_s)
atoms = []
for n in range(2, int(xf) + 1):
    p = vm_prime(n)
    if p and n < xf:
        atoms.append((arb(n).log(), arb(p).log()/arb(n).sqrt()))
csum = sum((c for _, c in atoms), arb(0))

# ---- (H) invisible-harmonic certificate (orchestrator, 2026-10-09): per prime p with powers p..p^K < x, replace the comb part
# F_p(theta) = -2 sum_k c_{p^k} cos(k theta) (theta = t log p) by -C_p + mu_p, mu_p >= 0 atomic with hat(mu_p)(0) = C_p and
# hat(mu_p)(k) = -c_{p^k}, 1 <= k <= K. The difference mu_p - C_p - F_p has only harmonics |k| > K, i.e. t-frequencies
# k log p >= (K+1) log p > log x = L, so it integrates to 0 against |fhat|^2 = ghat (supp g in [-L, L]): W is unchanged.
import os as _os
import numpy as _np
from scipy.optimize import linprog as _linprog
GRID = int(_os.environ.get('HGRID', '4000'))
byp = {}
for n in range(2, int(xf) + 1):
    p = vm_prime(n)
    if p and n < xf: byp.setdefault(p, []).append(n)
Cr = arb(0); hatoms = []
for p, ns in sorted(byp.items()):
    K = len(ns); lp = arb(p).log()
    ck = [lp/arb(p)**(arb(k)/2) for k in range(1, K+1)]
    ckf = [float(c) for c in ck]
    phis = _np.pi*_np.arange(GRID+1)/GRID; rho = _np.where((_np.arange(GRID+1) == 0) | (_np.arange(GRID+1) == GRID), 1.0, 2.0)
    Aeq = _np.array([rho*_np.cos(k*phis) for k in range(1, K+1)]); beq = -_np.array(ckf)
    lpres = _linprog(rho, A_eq=Aeq, b_eq=beq, bounds=(0, None), method='highs')
    assert lpres.status == 0, lpres.message
    sup = [int(g) for g in _np.argsort(-lpres.x)[:K]]
    sup = [g for g in sup if lpres.x[g] > 0] if K > 1 else [GRID]
    assert len(sup) == K, (p, K, sup)
    Msys = arb_mat(K, K); rhs = arb_mat(K, 1)
    for k in range(1, K+1):
        for i, g in enumerate(sup):
            Msys[k-1, i] = (1 if g in (0, GRID) else 2)*(pi*arb(k*g)/GRID).cos()
        rhs[k-1, 0] = -ck[k-1]
    vsol = Msys.solve(rhs)
    vs = [vsol[i, 0] for i in range(K)]
    assert all(v > 0 for v in vs), ('nonpositive weight', p, [v.str(5) for v in vs])
    Cp = sum(((1 if g in (0, GRID) else 2)*v for g, v in zip(sup, vs)), arb(0))
    Cr += Cp
    print('  (H) p = %d K = %d support phi/pi = %s weights %s C_p = %s (pointwise 2 sum c = %s)' % (p, K, [g/GRID for g in sup],
          [v.str(6) for v in vs], Cp.str(8), (2*sum(ck, arb(0))).str(8)), flush=True)
    for g, v in zip(sup, vs):
        phi = pi*arb(g)/GRID
        signs = (1,) if g in (0, GRID) else (1, -1)
        mass = 2*pi/lp*v
        mmax = int(float(T0)*float(lp)/(2*_np.pi)) + 2
        for s in signs:
            for m in range(-mmax, mmax + 1):
                tt = (s*phi + 2*pi*m)/lp
                if float(tt) > 1e-9 and float(tt) <= float(T0):
                    hatoms.append((tt, mass))
                elif abs(float(tt)) <= 1e-9:
                    raise SystemExit('atom at t = 0 not handled')
atom_mass_total = 2*sum((ms for _, ms in hatoms), arb(0))
print('(H) C_x = %s vs pointwise 2 sum c = %s; %d band atoms (t in (0, T0]), total band atom mass %s' % (Cr.str(10), (2*csum).str(10), len(hatoms), atom_mass_total.str(6)), flush=True)
def Minf_t(t):
    return acb(arb(1)/4, t/2).digamma().real - logpi
def Meff_z(z):
    w1 = acb(arb(1)/4) + acb(0, 1)*z/2
    w2 = acb(arb(1)/4) - acb(0, 1)*z/2
    v = (w1.digamma() + w2.digamma())/2 - acb(logpi) - acb(Cr)
    return v
alpha = arb((Minf_t(T0) - Cr).lower())
assert alpha > 0, 'alpha <= 0: increase T0'
print('x = %s  L = %s  N = %d  T0 = %s  atoms %d (sum c = %s)  alpha = %s  prec %d' % (x_s, L.str(10), N, T0_s, len(atoms), csum.str(8), alpha.str(12), prec), flush=True)
xmax = float(T0)*float(a)
assert N >= xmax + 1, 'need N >= T0*a for the tail monotonicity'

# ---- panels: (center, half-width, nodes, rho)
panels = []
nsmall = int(round(ts/0.5))
rho_s = arb('3.48')
for p in range(nsmall):
    panels.append((arb(0.25)*(2*p + 1), arb('0.25'), int(__import__('os').environ.get('NSMALL', '40')), rho_s))
tcur = ts
Tf = float(T0)
while tcur < Tf - 1e-9:
    h = min(hb, (Tf - tcur)/2)
    rr = bb/h
    rho = arb(rr + (rr*rr + 1)**0.5)        # rho - 1/rho = 2 bb/h  (semi-minor = bb)
    panels.append((arb(tcur) + arb(h), arb(h), nb, rho))
    tcur += 2*h
GL = {}
for (_, _, ng, _) in panels:
    if ng not in GL:
        GL[ng] = [arb.legendre_p_root(ng, k, weight=True) for k in range(ng)]
K_arc = 128
two_a = 2*a
def ellipse_bound(cen, hh, rho):
    best = arb(0)
    semi = hh*(rho + 1/rho)/2
    assert float(cen - semi) > 0 or float(hh) <= 0.25 + 1e-12
    r_arc = semi*pi/K_arc*arb('1.02')
    for k in range(K_arc):
        th = 2*pi*(arb(k) + arb(1)/2)/K_arc
        e = acb(th.cos(), th.sin())
        zc = acb(cen) + acb(hh)*(acb(rho)*e + 1/(acb(rho)*e))/2
        z = acb(arb(zc.real.mid(), r_arc), arb(zc.imag.mid(), r_arc))
        val = (Meff_z(z) - acb(alpha)).abs_upper()*(two_a*abs(z.imag).upper()).exp()
        ub = arb(val.upper())
        if ub > best: best = ub
    return best
t1 = time.time()
qerr = arb(0)
for (cen, hh, ng, rho) in panels:
    Mp = ellipse_bound(cen, hh, rho)
    qerr += hh*arb(64)/15/(rho*rho - 1)*rho**(-2*ng)*Mp
qerr = qerr/pi                        # entry (j,k) quadrature error <= c_j c_k qerr
print('%d panels, %d nodes; quadrature error unit %s (%.1f s)' % (len(panels), sum(p[2] for p in panels), qerr.str(5), time.time()-t1), flush=True)

cn = [(2*(2*n+1)*a).sqrt() for n in range(N+1)]
half = arb(1)/2
def jdirect(n, xx):
    # spherical j_n(xx) by Arb with adaptive working precision (the hypergeometric evaluation cancels for n >> xx)
    # evaluate at the exact midpoint; |j_n'(x)| <= 1 on the real line (j_n' = (n j_{n-1} - (n+1) j_{n+1})/(2n+1), |j_m| <= 1),
    # so the input radius adds at most rad(xx)
    m = arb(xx.mid()); r = xx.rad()
    p = prec
    while True:
        flint.ctx.prec = p
        v = (arb.pi()/(2*m)).sqrt()*m.bessel_j(arb(n) + arb(1)/2)
        flint.ctx.prec = prec
        if v.rad() < abs(v).upper()*arb(2)**(-prec+16) or p > 16*prec:
            break
        p *= 2
    # input-radius term. For 0 < y < sqrt(n(n+1)): j_n(y) > 0, j_n(y)/y^n is decreasing (derivative -y^{-n} j_{n+1} <= 0 before
    # the first zero of j_{n+1}), hence 0 <= j_n'(y) <= (n/y) j_n(y) and j_n(m+r) <= j_n(m) ((m+r)/m)^n.
    if float(m + r) < (n*(n+1))**0.5 and n >= 1:
        lo = m - r
        bound = r*n/lo*abs(v).upper()*((m + r)/m)**n
        return v + arb(0, arb(bound.upper()))
    return v + arb(0, r)
PREC_HI = 4*prec
def jtable(xx):
    # j_n(xx), n = 0..N, at the exact midpoint m = mid(xx):
    #  n <= ns (~ x + 12, oscillatory): j_0 = sin m/m, j_1 = sin m/m^2 - cos m/m, UPWARD recurrence at 4*prec bits (the ball
    #     wrapping growth <= ~2.42^n is absorbed by the extra precision);
    #  n > ns (evanescent, x^2 < n(n+1)): DOWNWARD recurrence from direct Arb values at N, N-1 (stable).
    # Input radius r: |j_n'| <= 1 on the real line; in the evanescent range the relative bound r n/(m-r) ((m+r)/m)^n (see jdirect).
    m = arb(xx.mid()); r = xx.rad()
    ns = min(N - 2, int(float(m)) + 12)
    j = [None]*(N+2)
    flint.ctx.prec = prec + int(1.4*float(m)) + 96     # wrapping growth of the upward recurrence <= ~2.42^n
    sn, cs = m.sin(), m.cos()
    j[0] = sn/m
    if ns >= 1:
        j[1] = sn/(m*m) - cs/m
    for n in range(1, ns):
        j[n+1] = (2*n+1)/m*j[n] - j[n-1]
    flint.ctx.prec = prec
    for n in range(0, ns + 1):
        j[n] = j[n] + arb(0, r)
    mm = arb(m.mid())
    _p = prec
    globals()['prec'] = prec + 400       # starting values need ~e^{120} extra relative accuracy (downward wrapping)
    top = [jdirect(N, mm), jdirect(N-1, mm)]
    globals()['prec'] = _p
    flint.ctx.prec = _p + 400          # ball wrapping in the downward recurrence grows like prod (1 + 2/c_n^2) <= e^{~120}
    j[N], j[N-1] = top
    for n in range(N-1, ns + 1, -1):
        j[n-1] = (2*n+1)/mm*j[n] - j[n+1]
    flint.ctx.prec = prec
    for n in range(ns + 1, N + 1):
        if r > 0:
            j[n] = j[n] + arb(0, arb((r*n/(m - r)*abs(j[n]).upper()*((m + r)/m)**n).upper()))
    return [j[n] + 0 for n in range(N+1)]
idx = {0: [n for n in range(N+1) if n % 2 == 0], 1: [n for n in range(N+1) if n % 2 == 1]}
G = {par: arb_mat(len(idx[par]), len(idx[par])) for par in (0, 1)}
t1 = time.time()
nodes = []
for (cen, hh, ng, rho) in panels:
    for (s, wgt) in GL[ng]:
        nodes.append((cen + hh*s, hh*wgt))
nb_q = len(nodes)
nodes = nodes + [(tt, None) for (tt, _) in hatoms]
hatoms_map = [(nb_q + j, ms/pi) for j, (tt, ms) in enumerate(hatoms)]
hval_atom = {}
chunk = 400
maxrad = arb(0)
for c0 in range(0, len(nodes), chunk):
    part = nodes[c0:c0+chunk]
    hval_atom = {i: ms/pi for i, (tt, ms) in enumerate([(None, None)]*0)}
    for i, (tt, wgt) in enumerate(part):
        if wgt is None:
            hval_atom[i] = dict(hatoms_map)[c0 + i]
    F = {par: arb_mat(len(idx[par]), len(part)) for par in (0, 1)}
    FW = {par: arb_mat(len(part), len(idx[par])) for par in (0, 1)}
    for i, (t, wgt) in enumerate(part):
        if wgt is None:
            wh = hval_atom[i]
        else:
            hval = Minf_t(t) - Cr - alpha
            wh = wgt*hval/pi
        jt = jtable(t*a)
        for par in (0, 1):
            for r, n in enumerate(idx[par]):
                v = cn[n]*jt[n]
                F[par][r, i] = v
                FW[par][i, r] = v*wh
        rr = max((jt[n].rad() for n in range(N+1)), key=lambda q: float(q))
        if rr > maxrad: maxrad = rr
    for par in (0, 1):
        G[par] = G[par] + F[par]*FW[par]
    print('  nodes %d-%d (%.1f s; max Bessel radius so far %s)' % (c0, c0+len(part)-1, time.time()-t1, arb(maxrad).str(3)), flush=True)
A = [((2*n+1)/(2*a)).sqrt()*2*a*((pi/a).sqrt()*(a/2).bessel_i(arb(n) + half)) for n in range(N+1)]
Gs = {}
for par in (0, 1):
    m = len(idx[par]); M = arb_mat(m, m)
    for r, j in enumerate(idx[par]):
        for c, k in enumerate(idx[par]):
            sgn = 1 if ((k - j)//2) % 2 == 0 else -1
            v = sgn*G[par][r, c] + arb(0, (cn[j]*cn[k]*qerr).upper())
            if j == k: v = v + alpha
            v = v + A[j]*A[k]*((-1)**j + (-1)**k)
            M[r, c] = v
    Gs[par] = M
print('matrix assembled (%.1f s total)' % (time.time() - t_start), flush=True)
import pickle
def ser(M):
    out = []
    for i in range(M.nrows()):
        row = []
        for j in range(M.ncols()):
            mid = M[i, j].mid(); rad = M[i, j].rad()
            mm_, me_ = mid.man_exp() if not mid.is_zero() else (0, 0)
            rm_, re_ = rad.man_exp() if not rad.is_zero() else (0, 0)
            row.append((int(mm_), int(me_), int(rm_), int(re_)))
        out.append(row)
    return out
pickle.dump({'x': x_s, 'N': N, 'T0': T0_s, 'prec': prec, 'alpha': alpha.mid().man_exp() if True else None,
             'G0': ser(Gs[0]), 'G1': ser(Gs[1])}, open('GH_x%s_N%d_T%s.pkl' % (x_s, N, T0_s), 'wb'))

# ---- tail
def dfact(m):
    r = arb(1)
    for k in range(m, 0, -2): r *= k
    return r
xa = T0*a
tau2 = arb(0)
M2 = N + 400
for n in range(N+1, M2+1):
    jn = jdirect(n, xa)
    tau2 += 2*(2*n+1)*a*jn*jn
rem = arb(0); n = M2 + 1
while True:
    term = 2*(2*n+1)*a*xa**(2*n)/dfact(2*n+1)**2
    rem += term
    if term < rem*arb('1e-30'): break
    n += 1
tau2 = tau2 + 2*rem
hmax = max(abs(Minf_t(arb(0)) - Cr - alpha).upper(), abs(Minf_t(T0) - Cr - alpha).upper())
epsK = (2*a).sqrt()*(1/(2*pi))*(2*T0*arb(hmax) + atom_mass_total)*tau2.sqrt()
pt2 = arb(0); n = N + 1; z = a/2
while True:
    term = (2*n+1)/(2*a)*(2*a*z**n*z.exp()/dfact(2*n+1))**2
    pt2 += term
    if n > N + 5 and term < pt2*arb('1e-30'): break
    n += 1
epsP = 4*(2*pt2).sqrt()*((a.exp() - (-a).exp())).sqrt()
eps = arb((epsK + epsP).upper())
print('tail: tau_N(T0)^2 <= %s (remainder %s), eps_K <= %s, eps_P <= %s' % (tau2.str(5), rem.str(3), epsK.str(5), epsP.str(5)), flush=True)

def ldl_ok(M, lam):
    n = M.nrows()
    Aa = [[M[i, j] for j in range(n)] for i in range(n)]
    for i in range(n): Aa[i][i] = Aa[i][i] - lam
    D = [None]*n; Lm = [[None]*n for _ in range(n)]
    minpiv = None
    for j in range(n):
        s = Aa[j][j]
        Lj = Lm[j]
        for k in range(j):
            s = s - Lj[k]*Lj[k]*D[k]
        D[j] = s
        if not (s > 0):
            return False, j, s
        if minpiv is None or s.lower() < minpiv.lower(): minpiv = s
        for i in range(j+1, n):
            t = Aa[i][j]; Li = Lm[i]
            for k in range(j):
                t = t - Li[k]*Lj[k]*D[k]
            Li[j] = t/s
    return True, n, minpiv
if __import__('os').environ.get('NOCERT'):
    print('NOCERT: matrix pickled; tail eps = %s, alpha = %s' % (eps.str(5), alpha.str(10)), flush=True)
    pickle.dump({'eps': (eps.mid().man_exp(), eps.rad().man_exp() if not eps.rad().is_zero() else (0,0)), 'alpha_str': alpha.str(40)}, open('tailH_x%s_N%d_T%s.pkl' % (x_s, N, T0_s), 'wb'))
    raise SystemExit(0)
import numpy as np
def precond(M):
    n = M.nrows()
    Mf = np.array([[float(M[i, j].mid()) for j in range(n)] for i in range(n)])
    w, V = np.linalg.eigh((Mf + Mf.T)/2)
    order = np.argsort(-w)                     # largest first, near-null directions factorised last
    V = V[:, order]
    Va = arb_mat([[arb(float(V[i, j])) for j in range(n)] for i in range(n)])
    H = Va.transpose()*M*Va
    S = Va.transpose()*Va
    return H, S
def ldl_ok2(H, S, lam):
    n = H.nrows()
    Aa = [[H[i, j] - lam*S[i, j] for j in range(n)] for i in range(n)]
    D = [None]*n; Lm = [[None]*n for _ in range(n)]
    for j in range(n):
        s_ = Aa[j][j]; Lj = Lm[j]
        for k in range(j):
            s_ = s_ - Lj[k]*Lj[k]*D[k]
        D[j] = s_
        if not (s_ > 0):
            return False, j, s_, D
        for i in range(j+1, n):
            t = Aa[i][j]; Li = Lm[i]
            for k in range(j):
                t = t - Li[k]*Lj[k]*D[k]
            Li[j] = t/s_
    return True, n, D[n-1], D
ok_all = True
for par in (0, 1):
    H, S = precond(Gs[par])
    lam = lam0[par]
    if float(lam.mid()) <= 0:
        lam = arb((2*eps*eps/(alpha - eps) + arb('1e-70')).upper())   # smallest lam the tail argument needs
    ok, j, piv, D = ldl_ok2(H, S, lam)
    tries = 0
    while not ok and tries < 6:
        print('  %s: not certified at lam = %s (pivot %d of %d: %s); lam /= 10' % ('even' if par == 0 else 'odd', lam.str(4), j, H.nrows(), piv.str(5)), flush=True)
        lam = lam/10; tries += 1
        ok, j, piv, D = ldl_ok2(H, S, lam)
    cond = (lam*(alpha - eps) >= eps*eps) and (alpha > eps)
    print('%s sector: V^T(G - %s I)V > 0: %s; last pivots %s; tail condition: %s' % ('even' if par == 0 else 'odd', lam.str(4), ok,
          ', '.join(D[k].str(5) for k in range(max(0, H.nrows()-3), H.nrows()) if D[k] is not None), cond), flush=True)
    ok_all = ok_all and ok and cond
print(('CERTIFIED: Weil form >= 0 on L^2[-L/2, L/2], x = %s' % x_s) if ok_all else 'NOT certified', '(%.1f s)' % (time.time()-t_start), flush=True)
