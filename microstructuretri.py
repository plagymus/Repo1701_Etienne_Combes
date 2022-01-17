import numpy as np

#crée une image carée avec deux spheres concentriques. le pixel vaut 0 dans la matrice, 1 dans la grande sphere et 2 dans la plus petite sphere

class booleanSpheres(object):
    def __init__(self,centers,radii,N):
        self.centers=centers
        self.radii=radii
        self.N = N
        self.d = len(N)
        self.phasemap = np.zeros(tuple(N),dtype=np.int8)
        self.initializePhaseMap()

    def initializePhaseMap(self):
        for i,center in enumerate(self.centers):
            r = self.radii[i]
            c = np.rint(center).astype(int)
            for s in np.ndindex(tuple(self.d*[2*int(r)+1])):
                ds = np.array(s)-int(r)
                n = c-ds
                if (n-center+0.5).dot(n-center+0.5)<= r*r:
                    self.phasemap[tuple(n % self.N)]=i+1
        
    def get_phase(self,n):
        return self.phasemap[n]  