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
    fig, (ax1) = plt.subplots(1)
    fig.subplots_adjust(wspace=0.5)
    im1 = ax1.imshow(vect , extent=[0, 1, 0, 1])
    ax1_divider = make_axes_locatable(ax1)
    cax1 = ax1_divider.append_axes("right", size="7%", pad="2%")
    cb1 = fig.colorbar(im1, cax=cax1)
    plt.title(label)
    plt.show()
    