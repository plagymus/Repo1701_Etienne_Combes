from microstructuretri import booleanSpheres
import numpy as np
from figures import figures


#fonction qui cree une microstructure carrée avec deux sphères concentriques plongées dans un domaine carré
#N0=nombre de pixels du coté de l'image, k1=conductivité grande sphère (matrice), k2=conductivité petite sphère (inclusion), k3=conductivité du "complément" cad le milieu homogénéisé

def variables0(N0, k1, k2, k3):
    
    #k3=k1*((k1+k2)-(k1-k2)*0.25)/((k1+k2)+(k1-k2)*0.25)   #khom
    
    d = 2
    N = np.array([N0, N0])  # size of the image (N0 x N0)
    centers = [[N0 / 2, N0 / 2],[N0 / 2, N0 / 2]]  # list of center coordinates in image units
    radius = [N0 / 2, N0 / 4]  # list of radius for circles
    microstructure = booleanSpheres(centers, radius, N)

    # initialize conductivity field
    ki = np.zeros(tuple(N) + (d,d))
    
    for i in np.ndindex(tuple(N)):
        chi = microstructure.phasemap[i]  # characteristic function of pixel i
        if chi==0:      
            ki[i] = k3* np.eye(d)   #khom
        if chi==1:
            ki[i] = k1* np.eye(d)   #matrice
        if chi==2:
            ki[i] = k2* np.eye(d)   #inclusion
        
    figures(ki[:,:,0,0], "Microstructure")  #affiche la composante 11 de la microstructure
    
    print("K3=",k3)

    return  N, ki,microstructure.phasemap


