src = open('wcert3_orig.py', encoding='utf-8').read()
def rep(old, new, cnt=1):
    global src
    assert src.count(old) == cnt, (old, src.count(old))
    src = src.replace(old, new)
HBLOCK = r'''
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
'''
rep("csum = sum((c for _, c in atoms), arb(0))\n", "csum = sum((c for _, c in atoms), arb(0))\n" + HBLOCK)
rep("""    v = (w1.digamma() + w2.digamma())/2 - acb(logpi)
    for ln, c in atoms:
        v = v - 2*c*(z*ln).cos()
    return v""", """    v = (w1.digamma() + w2.digamma())/2 - acb(logpi) - acb(Cr)
    return v""")
rep("alpha = arb((Minf_t(T0) - 2*csum).lower())", "alpha = arb((Minf_t(T0) - Cr).lower())")
rep("""        hval = Minf_t(t) - alpha
        for ln, c in atoms:
            hval = hval - 2*c*(t*ln).cos()
        wh = wgt*hval/pi""", """        if wgt is None:
            wh = hval_atom[i]
        else:
            hval = Minf_t(t) - Cr - alpha
            wh = wgt*hval/pi""")
rep("chunk = 400\n", "nodes = nodes + [(tt, None) for (tt, _) in hatoms]\nhval_atom = {}\nchunk = 400\n")
rep("""    part = nodes[c0:c0+chunk]
""", """    part = nodes[c0:c0+chunk]
    hval_atom = {i: ms/pi for i, (tt, ms) in enumerate([(None, None)]*0)}
    for i, (tt, wgt) in enumerate(part):
        if wgt is None:
            hval_atom[i] = dict(hatoms_map)[c0 + i]
""")
rep("nodes = nodes + [(tt, None) for (tt, _) in hatoms]\n", "nb_q = len(nodes)\nnodes = nodes + [(tt, None) for (tt, _) in hatoms]\nhatoms_map = [(nb_q + j, ms/pi) for j, (tt, ms) in enumerate(hatoms)]\n")
rep("hmax = max(abs(Minf_t(arb(0)) - 2*csum - alpha).upper(), abs(Minf_t(T0) + 2*csum - alpha).upper())",
    "hmax = max(abs(Minf_t(arb(0)) - Cr - alpha).upper(), abs(Minf_t(T0) - Cr - alpha).upper())")
rep("epsK = (2*a).sqrt()*(1/(2*pi))*2*T0*arb(hmax)*tau2.sqrt()",
    "epsK = (2*a).sqrt()*(1/(2*pi))*(2*T0*arb(hmax) + atom_mass_total)*tau2.sqrt()")
rep("open('G_x%s_N%d_T%s.pkl'", "open('GH_x%s_N%d_T%s.pkl'")
rep("open('tail_x%s_N%d_T%s.pkl'", "open('tailH_x%s_N%d_T%s.pkl'")
open('wcert_h.py', 'w', encoding='utf-8').write(src)
print('patched OK')
