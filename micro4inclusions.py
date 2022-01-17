from simplemicrostructuretri import variables0
import numpy as np
from figures import figures
from math import sqrt

#crée une microstrucure avec des inclusions de 2 contrastes différents. La microstructure est telle que khom=1 
#N0=taille de l'image,k1=condutivité inclusions 1, k2=conductivité inclusions 2
def variablesquad(N0,k1,k2):
    
    N = np.array([N0, N0])
    
    ka=(1/6)*(sqrt(25*k1**2-14*k1+25)-5*k1+5)
    kb=(1/6)*(sqrt(25*k2**2-14*k2+25)-5*k2+5)
    print("ka=",ka)
    print("kb=",kb)
    
    K1=variables0(int(N0/2),ka,k1,1)[1]
    K2=variables0(int(N0/2),kb,k2,1)[1]
    res1=np.concatenate([K1, K2], axis=1)
    res2=np.concatenate([K2, K1], axis=1)
    res=np.concatenate([res1, res2], axis=0)
    figures(res[:,:,0,0], "Microstructure")
    return N, res


