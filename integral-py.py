#!/usr/bin/env python3
"""Numerical integration: trapezoidal, Simpson's, Romberg, Gauss-Legendre."""
import sys,math

def trapezoid(f,a,b,n=1000):h=(b-a)/n;return h*(f(a)/2+sum(f(a+i*h)for i in range(1,n))+f(b)/2)
def simpson(f,a,b,n=1000):
    if n%2:n+=1
    h=(b-a)/n;s=f(a)+f(b)
    for i in range(1,n):s+=f(a+i*h)*(4 if i%2 else 2)
    return h*s/3
def romberg(f,a,b,n=8):
    R=[[0]*(n+1)for _ in range(n+1)]
    R[0][0]=trapezoid(f,a,b,1)
    for i in range(1,n+1):
        R[i][0]=trapezoid(f,a,b,2**i)
        for j in range(1,i+1):R[i][j]=(4**j*R[i][j-1]-R[i-1][j-1])/(4**j-1)
    return R[n][n]
def gauss_legendre(f,a,b,n=5):
    # 5-point rule
    nodes=[0,-0.5384693101,-0.9061798459,0.5384693101,0.9061798459]
    weights=[0.5688888889,0.4786286705,0.2369268851,0.4786286705,0.2369268851]
    mid=(b+a)/2;half=(b-a)/2
    return half*sum(w*f(mid+half*x)for x,w in zip(nodes,weights))

def main():
    if len(sys.argv)>1 and sys.argv[1]=="--test":
        f=lambda x:x**2
        assert abs(trapezoid(f,0,1)-1/3)<0.001
        assert abs(simpson(f,0,1)-1/3)<1e-10
        assert abs(romberg(f,0,1)-1/3)<1e-10
        assert abs(gauss_legendre(f,0,1)-1/3)<1e-6
        # sin(x) from 0 to pi = 2
        assert abs(simpson(math.sin,0,math.pi)-2)<1e-8
        assert abs(romberg(math.sin,0,math.pi)-2)<1e-10
        # e^x from 0 to 1 = e-1
        assert abs(simpson(math.exp,0,1)-(math.e-1))<1e-8
        print("All tests passed!")
    else:
        print(f"∫x² dx from 0 to 1 = {simpson(lambda x:x**2,0,1):.10f}")
if __name__=="__main__":main()
