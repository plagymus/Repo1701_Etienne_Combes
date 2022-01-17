import time
import numpy as np
from math import pi
from numba import jit, prange, njit
import pyfftw

@jit
def discrete_frequency(k,n,h):
    if 2*k>n:
        k-=n
    return 2*pi*k/(n*h)#2*np.tan(pi*k/n)/h

@jit
def initialize_frequencies(N,K,H,filter_level):
    d=len(N)
    frequencies = np.zeros((d,max(K),filter_level))
    for i in range(d):
        for ki in range(K[i]):
            for p in range(filter_level):
                kip = ki+N[i]*p
                frequencies[i][ki][p] = discrete_frequency(kip,N[i]*filter_level,H[i])
    return frequencies


@jit
def initialize_filters(N,K,n):
    d=len(N)
    filters = np.zeros((d,max(K),n))
    for i in range(d):
        for ki in range(K[i]):
            for p in range(n):
                z=ki/float(N[i])+p
                if z==0:
                    filters[i][ki][p] = 1
                else:
                    filters[i][ki][p] = (np.sin(pi*z)/(n*np.sin(pi*z/n)))**2
    return filters

#La fonction initializeGreen prend en argument N, la taille de la grille/image (c'est un 'tuple' de taille d où d est la dimension. Ex : en dimension 2, pour une image 128x128, N=(128,128)) et renvoie les éléments dont on a besoin pour utiliser la fonction operate_field
def initializeGreen(N,filter_level=2):
    d = len(N)
    numComp=d
    H=1./N
    K=N.copy()
    K[-1] = N[-1]//2+1
    Npadded = N.copy()
    Npadded[-1] = K[-1]*2
    full_array = pyfftw.n_byte_align_empty(np.append(Npadded,numComp), pyfftw.simd_alignment,'float64')
    field = full_array[...,:N[-1],:]
    field_fourier = full_array.ravel().view('complex128').reshape(np.append(K,numComp))
    fft = pyfftw.FFTW(field, field_fourier, axes=range(d))
    ifft = pyfftw.FFTW(field_fourier, field, axes=range(d),direction='FFTW_BACKWARD')
    frequencies = initialize_frequencies(N,K,H,filter_level)
    filters = initialize_filters(N,K,filter_level)
    tupleK=tuple(K[1:])
    tupleFilters = np.zeros((filter_level**d,d),dtype=np.int64)#tuple(d*[filter_level])
    i=0
    for shift in np.ndindex(tuple(d*[filter_level])):
        tupleFilters[i] = np.array(shift)
        i+=1
    return field,field_fourier,fft,ifft,tupleK,frequencies,filters,tupleFilters

@jit
def get_matrix_fourier_filtered(k,N,frequencies,filter_level,filters,tupleFilters,k0):
    """
    k : multi-index in frequency grid
    """
    d=len(k)
    Gamma = np.zeros((d,d),dtype=np.complex128)
    q = np.empty(d)
    #for p in np.ndindex(tupleFilters):
    #    kp = k+N*np.array(p)
    for j in range(len(tupleFilters)):
        p = tupleFilters[j]
        kp = k+N*p
        f = 1
        for i in range(d):
            q[i] = frequencies[i][k[i]][p[i]]
            f*=filters[i][k[i]][p[i]]
        g=np.outer(q,q)
        qk0q = q.dot((k0).dot(q))
        if not(qk0q==0):# or nullify_frequency(kp/filter_level,N)):
            Gamma+=f*g/(qk0q)
    return Gamma

@jit
def operate_fourier(tau,k,N,frequencies,filter_level,filters,tupleFilters,k0):
    """
    tau : d-dimensionnal array of polarization
    k : multi-index in frequency grid
    """
    return -get_matrix_fourier_filtered(k,N,frequencies,filter_level,filters,tupleFilters,k0).dot(tau)


@jit(parallel=False)
def operate_fourier_field(x,y,tupleK,N,frequencies,filter_level,filters,tupleFilters,k0):
    """ 
    x : input, field of comp-dimensionnal array
    y : output, field of comp-dimensionnal array (operates in place if x = y)
    tupleK : a tuple of fourier grid dimensions
    """
    for kx in prange(N[0]):
        for kyz in np.ndindex(tupleK):
            k = (kx,)+kyz
            y[k] = operate_fourier(x[k],np.array(k),N,frequencies,filter_level,filters,tupleFilters,k0)
        #do not use tuple as arguments in parallel mode -> remove tupleFIlters
        #if not enough, remove tupleK and use scalar index to flattened x,y

#La fonction operate_field permet de calculer, pour un champ 'x' (qui représente un champ de polarisation), la déformation créée par la polarisation 'x'. C'est-à-dire, -Gamma(x). N est la taille de la grille/image. L'argument k0 est la conductivité du milieu de référence utilisé pour le calcul. Les autres arguments sont ceux qui sont issus de la fonction initializeGreen
def operate_field(x,yFourier,fft,ifft,tupleK,N,frequencies,filter_level,filters,tupleFilters,k0):
    xFourier=fft(x)
    operate_fourier_field(xFourier,yFourier,tupleK,N,frequencies,filter_level,filters,tupleFilters,k0)        
    return ifft(yFourier)

