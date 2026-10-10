# Where can the glued base reach?  Needs: T_C(x) = smallest T with M_inf(T) >= C_x + 1  (prime side safely positive),
# sampling region [2 pi x^4, R] inside the rigorously known zero range, R >= e*T_C.  (x not a prime power; C_x separable.)
import numpy as np, mpmath as mp
from math import log, pi, sqrt
def isp(p): return p>1 and all(p%q for q in range(2,int(p**.5)+1))
def Minf(t): return float(mp.re(mp.digamma(0.25+0.5j*t))-mp.log(mp.pi))
def T_for(level):
    lo,hi=1.0,1e40
    for _ in range(300):
        mid=(lo*hi)**.5
        if Minf(mid)>=level: hi=mid
        else: lo=mid
    return hi
def Cx(x):
    C=0.0
    for p in range(2,int(x)+1):
        if not isp(p): continue
        K=int(np.floor(log(x)/log(p)-1e-12))
        if K<1: continue
        c=[log(p)*p**(-k/2) for k in range(1,K+1)]
        T=np.zeros((K+1,K+1))
        for i in range(K+1):
            for j in range(K+1):
                if i!=j: T[i,j]=-c[abs(i-j)-1]
        C+=-np.linalg.eigvalsh(T)[0]
    return C
print("   x      C_x   C_x/sqrt(x)   T_C          2pi x^4      need R >= max(e*T_C, 10*2pi x^4)")
for x in [20.5,50.5,100.5,150.5,200.5,250.5,300.5,400.5]:
    C=Cx(x); TC=T_for(C+1)
    print(f"{x:6.1f} {C:8.3f} {C/sqrt(x):8.3f}   {TC:12.3e} {2*pi*x**4:12.3e}   {max(np.e*TC,20*pi*x**4):12.3e}")
