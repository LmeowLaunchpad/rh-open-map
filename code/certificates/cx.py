# Separable invisible-harmonic (Caratheodory-Toeplitz) constant C_x vs Zhu's pointwise comb sup A_x.
import numpy as np, mpmath as mp
from math import log, sqrt, pi
def primes(n): return [p for p in range(2,n+1) if all(p%q for q in range(2,int(p**.5)+1))]
def Minf(t): return float(mp.re(mp.digamma(0.25+0.5j*t)) - mp.log(mp.pi))
def T_for(level):  # smallest T with Minf(T) >= level (Minf increasing for t>0)
    lo,hi=0.0,1e12
    for _ in range(200):
        mid=(lo+hi)/2
        if Minf(mid)>=level: hi=mid
        else: lo=mid
    return hi
print(" x     A_x(Zhu) C_x(sep) | T0 Zhu      T0 sep   (alpha=0.583)")
for x in [5,7,8,9,11,13,16,17,19,23,25,29,31]:
    X=x-1e-9
    A=0.0; C=0.0
    for p in primes(x):
        K=int(log(X)/log(p))
        if K<1: continue
        c=[log(p)*p**(-k/2) for k in range(1,K+1)]   # c_{p^k} = log p / p^{k/2}
        A+=2*sum(c)
        Tm=np.zeros((K+1,K+1))
        for i in range(K+1):
            for j in range(K+1):
                if i!=j: Tm[i,j]=-c[abs(i-j)-1]
        C+=-np.linalg.eigvalsh(Tm)[0]
    print(f"{x:3d}  {A:8.3f} {C:8.3f} | {T_for(A+0.583):10.0f} {T_for(C+0.583):9.0f}")
