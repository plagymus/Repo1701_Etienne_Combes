# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 11:44:21 2022

@author: Etienne
"""

#tracé des différents graphes

import matplotlib.pyplot as plt
from Ada1parametre import ada1
from math import sqrt
import numpy as np



p=128

k1=1e4
k2=1
ka=(1/6)*(sqrt(25*k1**2-14*k1+25)-5*k1+5)
kb=(1/6)*(sqrt(25*k2**2-14*k2+25)-5*k2+5)

    
it=500
micro=0

if micro==0:
   KA=max(k1,ka,1)
   KB=min(k1,ka,1) 

if micro==1:
    KA=max(k1,k2,ka,kb,1)
    KB=min(k1,k2,ka,kb,1)


kms=0.5*(KA+KB)
kem=sqrt(KA*KB)

K=["MS kms","EM kem","EM khom","EMada1 khom","MSada1 khom","MSada1bis khom","MSada1ter khom","MS GC khom"]
header1="MS kms, EM kem, EM khom, EMada1 khom, MSada1 khom, MSada1bis khom,MSada1ter khom, MS GC khom"

liste=[3,4,5,6,7] #algos que l'on souhaite faire tourner

K=[K[i] for i in liste]
header = ", ".join(K)

legend=str(k1)+" "+str(k2)+" "+str(p)+" Comparaison"

U0=[]
U1=[]
U2=[]
ER=[]
ALPHA=[]

for i in liste:
    
    if i==0:
        U=ada1(p,micro,k1,k2,kms,0,0,it)  #MS simple
    if i==1:
        U=ada1(p,micro,k1,k2,kem,1,0,it)  #EM simple kem
    if i==2:
        U=ada1(p,micro,k1,k2,1,1,0,it)  #EM simple khom
    if i==3:
        U=ada1(p,micro,k1,k2,1,1,1,it)    #EM ame 1 par
    if i==4:
        U=ada1(p,micro,k1,k2,1,0,1,it)    #MS ame 1 par 
    if i==5:
        U=ada1(p,micro,k1,k2,1,0,2,it)    #MS ame 1 par bis
    if i==6:
        U=ada1(p,micro,k1,k2,1,0,3,it)    #MS ame 1 par ter
    if i==7:
        U=ada1(p,micro,k1,k2,1,0,4,it)    #MS GC

    U0.append(U[0])
    U1.append(U[1])
    U2.append([abs(x-y) for x,y in zip(U[0],U[1])])
    ER.append(U[2])
    ALPHA.append(U[3])


itplot1=0
itplot2=5000


    
for i in range(len(K)):
    plt.plot([t for t in range(len(U0[i]))[itplot1:itplot2]],[abs(x-1) for x in U0[i][itplot1:itplot2]],label=K[i])
plt.title("(K-Kth/Kth) <CA>  Micro: "+str(micro)+" Discrétisation:"+str(p)+"  k1: "+str(k1)+"  k2: "+str(k2))
ax = plt.subplot(111)
ax.set_yscale('log')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()

for i in range(len(K)):
    plt.plot([t for t in range(len(U1[i]))[itplot1:itplot2]],[abs(x-1) for x in U1[i][itplot1:itplot2]],label=K[i])
plt.title("(K-Kth/Kth) <tACA>  Micro: "+str(micro)+" Discrétisation:"+str(p)+"  k1: "+str(k1)+"  k2: "+str(k2))
ax = plt.subplot(111)
ax.set_yscale('log')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()

for i in range(len(K)):
    plt.plot([t for t in range(len(U2[i]))[itplot1:itplot2]],U2[i][itplot1:itplot2],label=K[i])   
plt.title("Ecart Micro: "+str(micro)+"  Discrétisation:"+str(p)+"  k1: "+str(k1)+"  k2: "+str(k2))
ax = plt.subplot(111)
ax.set_yscale('log')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()


for i in range(len(K)):
    plt.plot([t for t in range(len(ER[i]))[itplot1:itplot2]],ER[i][itplot1:itplot2],label=K[i])   
plt.title("Erreur  Micro: "+str(micro)+"  Discrétisation:"+str(p)+"  k1: "+str(k1)+"  k2: "+str(k2))
ax = plt.subplot(111)
ax.set_yscale('log')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()

np.savetxt("calculs/"+legend+"CA.csv", np.array(np.transpose(U0)),header=header)
np.savetxt("calculs/"+legend+"tACA.csv", np.array(np.transpose(U1)),header=header)
np.savetxt("calculs/"+legend+"X.csv", np.array(np.transpose(ER)),header=header)
np.savetxt("calculs/"+legend+"ALPHA.csv", np.array(np.transpose(ALPHA)),header=header)




