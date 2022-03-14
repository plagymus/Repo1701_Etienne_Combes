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

k1=1e-5
k2=1
ka=(1/6)*(sqrt(25*k1**2-14*k1+25)-5*k1+5)
kb=(1/6)*(sqrt(25*k2**2-14*k2+25)-5*k2+5)

    
it=2
micro=0
ini=0 #initialisation de Voigt ou Reuss

if micro==0:
   KA=max(k1,ka,1)
   KB=min(k1,ka,1) 

if micro==1:
    KA=max(k1,k2,ka,kb,1)
    KB=min(k1,k2,ka,kb,1)


kms=0.5*(KA+KB)
kem=sqrt(KA*KB)

K=["MS kms V","EM kem V","EM khom V","EMada1 khom V","MSada1 khom V","MSada1bis khom V","MSada1ter khom V","MS GC khom V","MS kms R","EM kem R","EM khom R","EMada1 khom R","MSada1 khom R","MSada1bis khom R","MSada1ter khom R","MS GC khom R","MS kms M","EM kem M","EM khom M","EMada1 khom M","MSada1 khom M","MSada1bis khom M","MSada1ter khom M","MS GC khom M"]

header1="MS kms, EM kem, EM khom, EMada1 khom, MSada1 khom, MSada1bis khom,MSada1ter khom, MS GC khom"
header1="MS kms V, EM kem V, EM khom V, EMada1 khom V, MSada1 khom V, MSada1bis khom V, MSada1ter khom V, MS GC khom V, MS kms R, EM kem R, EM khom R, EMada1 khom R, MSada1 khom R, MSada1bis khom R, MSada1ter khom R, MS GC khom R, MS kms M, EM kem M, EM khom M, EMada1 khom M, MSada1 khom M, MSada1bis khom M, MSada1ter khom M, MS GC khom M"

liste=[3,5] #algos que l'on souhaite faire tourner


K=[K[i] for i in liste]
header = ", ".join(K)
print(header)
legend=str(k1)+" "+str(k2)+" "+str(p)+" Comparaison"

U0=[]
U1=[]
U2=[]
ER=[]
ALPHA=[]

for i in liste:
    
    if i==0:
        U=ada1(p,micro,k1,k2,kms,0,0,it,0)  #MS simple V
    if i==1:
        U=ada1(p,micro,k1,k2,kem,1,0,it,0)  #EM simple kem V
    if i==2:
        U=ada1(p,micro,k1,k2,1,1,0,it,0)  #EM simple khom V
    if i==3:
        U=ada1(p,micro,k1,k2,1,1,1,it,0)    #EM ame 1 par V
    if i==4:
        U=ada1(p,micro,k1,k2,1,0,1,it,0)    #MS ame 1 par V 
    if i==5:
        U=ada1(p,micro,k1,k2,1,0,2,it,0)    #MS ame 1 par bis V
    if i==6:
        U=ada1(p,micro,k1,k2,1,0,3,it,0)    #MS ame 1 par ter V
    if i==7:
        U=ada1(p,micro,k1,k2,1,0,4,it,0)    #MS GC V
    if i==8:
        U=ada1(p,micro,k1,k2,kms,0,0,it,1)  #MS simple R
    if i==9:
        U=ada1(p,micro,k1,k2,kem,1,0,it,1)  #EM simple kem R
    if i==10:
        U=ada1(p,micro,k1,k2,1,1,0,it,1)  #EM simple khom R
    if i==11:
        U=ada1(p,micro,k1,k2,1,1,1,it,1)    #EM ame 1 par R
    if i==12:
        U=ada1(p,micro,k1,k2,1,0,1,it,1)    #MS ame 1 par  R
    if i==13:
        U=ada1(p,micro,k1,k2,1,0,2,it,1)    #MS ame 1 par bis R
    if i==14:
        U=ada1(p,micro,k1,k2,1,0,3,it,1)    #MS ame 1 par ter R
    if i==15:
        U=ada1(p,micro,k1,k2,1,0,4,it,1)    #MS GC R
    if i==16:
        U=ada1(p,micro,k1,k2,kms,0,0,it,2)  #MS simple M
    if i==17:
        U=ada1(p,micro,k1,k2,kem,1,0,it,2)  #EM simple kem M
    if i==18:
        U=ada1(p,micro,k1,k2,1,1,0,it,2)  #EM simple khom M
    if i==19:
        U=ada1(p,micro,k1,k2,1,1,1,it,2)    #EM ame 1 par M
    if i==20:
        U=ada1(p,micro,k1,k2,1,0,1,it,2)    #MS ame 1 par M 
    if i==21:
        U=ada1(p,micro,k1,k2,1,0,2,it,2)    #MS ame 1 par bis M
    if i==22:
        U=ada1(p,micro,k1,k2,1,0,3,it,2)    #MS ame 1 par ter M
    if i==23:
        U=ada1(p,micro,k1,k2,1,0,4,it,2)    #MS GC M

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




