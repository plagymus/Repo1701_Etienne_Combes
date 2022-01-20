# -*- coding: utf-8 -*-
"""
Created on Fri Jan 14 11:44:21 2022

@author: Etienne
"""

#tracé des différents graphes

import matplotlib.pyplot as plt
from Ada1parametre import ada1
from math import sqrt



p=256

c1=1E-5
c2=1E-5
k1=c1
k2=c2
ka=(1/6)*(sqrt(25*k1**2-14*k1+25)-5*k1+5)
kb=(1/6)*(sqrt(25*k2**2-14*k2+25)-5*k2+5)

    
it=400
micro=1

itplot1=0
itplot2=100

KA=max(c1,c2,ka,kb,1)
KB=min(c1,c2,ka,kb,1)


kms=0.5*(KA+KB)
kem=sqrt(KA*KB)

K=["MS kms","EM kem","EM khom","MSada1 khom","EMada1 khom"]



U0=[]
U1=[]
U2=[]
ER=[]

for i in [0,1,2,3,4]:
    
    if i==0:
        U=ada1(p,micro,c1,c2,kms,0,0,it)  #MS simple
    if i==1:
        U=ada1(p,micro,c1,c2,kem,1,0,it)  #EM simple
    if i==2:
        U=ada1(p,micro,c1,c2,1,1,0,it)  #EM simple
    if i==3:
        U=ada1(p,micro,c1,c2,1,0,1,it)    #MS ame 1 par
    if i==4:
        U=ada1(p,micro,c1,c2,1,1,1,it)    #EM ame 1 par

    U0.append(U[0])
    U1.append(U[1])
    U2.append([abs(x-y) for x,y in zip(U[0],U[1])])
    ER.append(U[2])


itplot1=0
itplot2=5000


    
for i in range(len(K)):
    plt.plot([t for t in range(len(U0[i]))[itplot1:itplot2]],[abs(x-1) for x in U0[i][itplot1:itplot2]],label=K[i])
plt.title("Evolutions de log(Cth-Chom/Cth) <CA> pour les différents algorithmes  Discrétisation:"+str(p)+"  C1: "+str(c1)+"  C2: "+str(c2)+"  k3: ")
ax = plt.subplot(111)
ax.set_yscale('log')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()

for i in range(len(K)):
    plt.plot([t for t in range(len(U1[i]))[itplot1:itplot2]],[abs(x-1) for x in U1[i][itplot1:itplot2]],label=K[i])
plt.title("Evolutions de log(Cth-Chom/Cth) <tACA>   Discrétisation:"+str(p)+"  C1: "+str(c1)+"  C2: "+str(c2)+"  k3: ")
ax = plt.subplot(111)
ax.set_yscale('log')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()

for i in range(len(K)):
    plt.plot([t for t in range(len(U2[i]))[itplot1:itplot2]],U2[i][itplot1:itplot2],label=K[i])   
plt.title("Evolutions de l'ecart    Discrétisation:"+str(p)+"  C1: "+str(c1)+"  C2: "+str(c2)+"  k3: ")
ax = plt.subplot(111)
ax.set_yscale('log')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()


for i in range(len(K)):
    plt.plot([t for t in range(len(ER[i]))[itplot1:itplot2]],ER[i][itplot1:itplot2],label=K[i])   
plt.title("Evolutions de l'erreur    Discrétisation:"+str(p)+"  C1: "+str(c1)+"  C2: "+str(c2)+"  k3: ")
ax = plt.subplot(111)
ax.set_yscale('log')
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.show()



