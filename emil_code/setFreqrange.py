import numpy as n

def kappa1D(N,delta):
    
    kapparange = n.zeros((N,1))
    if N % 2 == 0:
	kapparange = n.arange(-N/2,N/2, dtype=float)/(delta*N)
    else:
        kapparange = n.arange(-(N-1)/2,(N-1)/2 + 1, dtype=float)/(delta*N)
    
    return kapparange