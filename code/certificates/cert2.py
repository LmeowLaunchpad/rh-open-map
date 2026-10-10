# E4 certification stage v2 (2026-10-09): certify G_sector - lam I > 0 for the Galerkin matrices pickled by wcert3.py (NOCERT=1).
# Basis change V = [V_bulk | X_ladder]: V_bulk = double-precision eigenvectors with eigenvalue >= theta, X_ladder = near-null subspace
# refined by block inverse iteration at full precision (arb solve on the exact midpoint matrix). Then, rigorously in Arb:
#   H = V^T G V, S = V^T V (G with its radii); B = H_bb - lam S_bb certified PD by Weyl (min diagonal - ||offdiag||_F > 0);
#   Schur complement Sch = (H_ll - lam S_ll) - C B^{-1} C^T (C = H_lb - lam S_lb, B^{-1} C^T by Arb solve) certified PD by ball LDL^T.
#   (B > 0 and Sch > 0  <=>  H - lam S > 0  <=>  G - lam I > 0, V invertible by congruence.)
# Usage: python cert2.py G_pickle tail_pickle [theta] [iters]
import sys, pickle, time
import numpy as np
import flint
from flint import arb, arb_mat, fmpz

t_start = time.time()
D = pickle.load(open(sys.argv[1], 'rb'))
T = pickle.load(open(sys.argv[2], 'rb'))
theta = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-5
iters = int(sys.argv[4]) if len(sys.argv) > 4 else 12
prec = D['prec']; flint.ctx.prec = prec
def mk(m, e):
    return arb(fmpz(int(m)))*arb(2)**int(e) if m != 0 else arb(0)
def deser(rows, with_rad=True):
    n = len(rows)
    M = arb_mat(n, n); Mm = arb_mat(n, n)
    for i in range(n):
        for j in range(n):
            mm, me, rm, re = rows[i][j]
            mid = mk(mm, me); rad = mk(rm, re)
            Mm[i, j] = mid
            M[i, j] = arb(mid, rad) if with_rad else mid
    return M, Mm
alpha = arb(T['alpha_str'])
(em, ee), (rm, re) = T['eps']
eps = arb((mk(em, ee) + mk(rm, re)).upper())
lam = arb((2*eps*eps/(alpha - eps) + arb('1e-70')).upper())
print('x = %s N = %s T0 = %s prec %d: alpha = %s, eps <= %s, lam = %s' % (D['x'], D['N'], D['T0'], prec, alpha.str(8), eps.str(4), lam.str(4)), flush=True)

def ldl(Ms):
    n = len(Ms)
    Dd = [None]*n; Lm = [[None]*n for _ in range(n)]
    for j in range(n):
        s_ = Ms[j][j]; Lj = Lm[j]
        for k in range(j):
            s_ = s_ - Lj[k]*Lj[k]*Dd[k]
        Dd[j] = s_
        if not (s_ > 0):
            return False, j, Dd
        for i in range(j+1, n):
            t = Ms[i][j]; Li = Lm[i]
            for k in range(j):
                t = t - Li[k]*Lj[k]*Dd[k]
            Li[j] = t/s_
    return True, n, Dd

ok_all = True
for key, name in (('G0', 'even'), ('G1', 'odd')):
    G, Gm = deser(D[key])
    n = G.nrows()
    Gf = np.array([[float(Gm[i, j]) for j in range(n)] for i in range(n)])
    w, V = np.linalg.eigh((Gf + Gf.T)/2)
    bulk = [i for i in range(n) if w[i] >= theta]; lad = [i for i in range(n) if w[i] < theta]
    K = len(lad)
    # refine the near-null subspace: X <- orthonormalise(Gm^{-1} X), at full precision (midpoints only)
    X = arb_mat([[arb(float(V[i, j])) for j in lad] for i in range(n)])
    for it in range(iters):
        Y = Gm.solve(X)
        # Gram-Schmidt on midpoints (mpmath-free, in arb then .mid())
        cols = [[Y[i, j].mid() for i in range(n)] for j in range(K)]
        Q = []
        for c in cols:
            v = list(c)
            for q in Q:
                dot = sum((v[i]*q[i] for i in range(n)), arb(0))
                v = [(v[i] - dot*q[i]).mid() for i in range(n)]
            nr = sum((vi*vi for vi in v), arb(0)).sqrt()
            Q.append([(vi/nr).mid() for vi in v])
        X = arb_mat([[Q[j][i] for j in range(K)] for i in range(n)])
    X = arb_mat([[X[i, K - 1 - j] for j in range(K)] for i in range(n)])   # largest first, smallest last
    Vb = [[arb(float(V[i, j])) for j in bulk] for i in range(n)]
    Vfull = arb_mat([[Vb[i][jj] for jj in range(len(bulk))] + [X[i, j] for j in range(K)] for i in range(n)])
    H = Vfull.transpose()*G*Vfull
    S = Vfull.transpose()*Vfull
    nb = len(bulk)
    A = arb_mat(n, n)
    for i in range(n):
        for j in range(n):
            A[i, j] = H[i, j] - lam*S[i, j]
    # bulk block by Weyl
    mind = min((A[i, i] for i in range(nb)), key=lambda q: float(q.lower()))
    off2 = arb(0)
    for i in range(nb):
        for j in range(nb):
            if i != j: off2 += A[i, j]*A[i, j]
    offF = off2.sqrt()
    bulk_ok = (mind - offF) > 0
    # Schur complement on the ladder block
    B = arb_mat([[A[i, j] for j in range(nb)] for i in range(nb)])
    Ct = arb_mat([[A[i, nb + j] for j in range(K)] for i in range(nb)])
    Z = B.solve(Ct)
    Sch = [[A[nb + a, nb + b] - sum((Ct[i, a]*Z[i, b] for i in range(nb)), arb(0)) for b in range(K)] for a in range(K)]
    ok, j, Dd = ldl(Sch)
    ok_all = ok_all and bool(bulk_ok) and ok
    piv = ', '.join(Dd[k].str(5) for k in range(max(0, K-4), K) if Dd[k] is not None)
    print('%s: n = %d, bulk %d (min diag %s, ||offdiag||_F %s, PD %s), ladder K = %d: Schur LDL %s, last pivots %s' %
          (name, n, nb, mind.str(4), offF.str(3), bulk_ok, K, ok, piv), flush=True)
print(('CERTIFIED: G - lam I > 0 in both sectors with lam(alpha - eps) >= eps^2  =>  Weil form >= 0 on the window x = %s' % D['x'])
      if ok_all else 'NOT certified', '(%.1f s)' % (time.time() - t_start), flush=True)
