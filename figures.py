# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 15:27:26 2021

@author: Etienne
"""

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

#la fonction figure permet de tracer une image (de la microstructure, du champ de température)
#vect=image (tableau 2D), label=légende

def figures(vect, label):
    imgplot = plt.imshow(vect)
    plt.title(label)
    plt.colorbar()
    plt.show()