#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
from greenOperatorConductionNumba import operate_field, initializeGreen
from simplemicrostructuretri import variables0
from figures import figures

#La fonction Ref est un algorithme qui permet de calculer la conductivit homogeneise (kref) avec EM 3 parametres. Ne marche que pour la microstructure simple avec 3 phases
#N0:nb pixels de l'image (carée)
#km=conductivité matrice
#k2,k2=conductivités inclusions
#Prec est la précision ou le nombre d'itération pour le critère d'arrêt de la récursion

def EMada3(N0,k1,k2,K0,Prec):

    
    CA,ACA,ER=[],[],[]
    
    khom=k1*((k1+k2)-(k1-k2)*0.25)/((k1+k2)+(k1-k2)*0.25)   #khom
    U=variables0(N0,k1,k2,khom)                             #initialisation de la microstructure
    
    N=U[0]
    ki=U[1]
    chi=U[2]
    d=len(N)
    filter_level=2
    k0=K0*np.eye(d)    
            
    A=np.zeros(tuple(N)+(d,d))
    
    for j in range(d):                   
        
        it=0
        
        Eps_field = np.zeros(tuple(N)+(d,))
        E_field = np.zeros(tuple(N)+(d,)) 
        X= np.zeros(tuple(N)+(d,))
        Y1= np.zeros(tuple(N)+(d,))
        Y2= np.zeros(tuple(N)+(d,))
        Y3= np.zeros(tuple(N)+(d,))
        chi1=np.zeros(tuple(N))   #fonctions indicatrices de la matrice ou des phases 1 et 2
        chi2=np.zeros(tuple(N))
        chi3=np.zeros(tuple(N))              

        
        for i in np.ndindex(tuple(N)):                     #Champ E (chargement)
                E_field[i+(j,)]=1.
        Eps_field=np.copy(E_field)                         #initialisation: Eps=chargement
        
        field,field_fourier,fft,ifft,tupleK,frequencies,filters,tupleFilters = initializeGreen(N,filter_level=filter_level)
        field1,field_fourier1,fft1,ifft1,tupleK1,frequencies1,filters1,tupleFilters1 = initializeGreen(N,filter_level=filter_level)
        field2,field_fourier2,fft2,ifft2,tupleK2,frequencies2,filters2,tupleFilters2 = initializeGreen(N,filter_level=filter_level)
        field3,field_fourier3,fft3,ifft3,tupleK3,frequencies3,filters3,tupleFilters3 = initializeGreen(N,filter_level=filter_level)

        for i in np.ndindex(tuple(N)):                         
                X[i]=np.dot(ki[i]-k0, Eps_field[i])
                if chi[i]==0:
                    chi3[i]=1
                if chi[i]==1:
                    chi1[i]=1
                if chi[i]==2:
                    chi2[i]=1
                
                
        X=operate_field(X,field_fourier,fft,ifft,tupleK,N,frequencies,filter_level,filters,tupleFilters,k0)
        Erreur0=np.linalg.norm(X)
        X_0X_0=0
        for i in np.ndindex(tuple(N)):                                                             
            X_0X_0+=np.dot(X[i],X[i])
        while it<Prec:                                    
            
            for i in np.ndindex(tuple(N)):
                Y1[i]=np.dot((ki[i]-k0)*chi1[i],X[i])
                Y2[i]=np.dot((ki[i]-k0)*chi2[i],X[i])
                Y3[i]=np.dot((ki[i]-k0)*chi3[i],X[i])
                
            Y_1=-operate_field(Y1,field_fourier1,fft1,ifft1,tupleK1,N,frequencies1,filter_level,filters,tupleFilters1,k0)   
            Y_2=-operate_field(Y2,field_fourier2,fft2,ifft2,tupleK2,N,frequencies2,filter_level,filters,tupleFilters2,k0) 
            Y_3=-operate_field(Y3,field_fourier2,fft2,ifft2,tupleK2,N,frequencies2,filter_level,filters,tupleFilters2,k0)  

            for i in np.ndindex(tuple(N)):
                Y1[i]=chi1[i]*X[i]+Y_1[i]
                Y2[i]=chi2[i]*X[i]+Y_2[i]
                Y3[i]=chi3[i]*X[i]+Y_3[i]
                
            Y3Y3=0
            Y1Y1=0
            Y2Y2=0
            
            Y3Y1=0
            Y3Y2=0
            Y1Y2=0

            XY3=0
            XY1=0
            XY2=0
            XX=0
            for i in np.ndindex(tuple(N)):                                                   
                Y3Y3+=np.dot(Y3[i],Y3[i])
                Y1Y1+=np.dot(Y1[i],Y1[i])
                Y2Y2+=np.dot(Y2[i],Y2[i])

                Y3Y1+=np.dot(Y3[i],Y1[i])
                Y3Y2+=np.dot(Y3[i],Y2[i])
                Y1Y2+=np.dot(Y1[i],Y2[i])

                XY3+=np.dot(X[i],Y3[i])                
                XY1+=np.dot(X[i],Y1[i])
                XY2+=np.dot(X[i],Y2[i])
                
                XX+=np.dot(X[i],X[i])
                
            M=np.zeros((3,3))
            M[0][0]=Y3Y3
            M[1][1]=Y1Y1
            M[2][2]=Y2Y2
            M[0][1]=M[1][0]=Y3Y1
            M[1][2]=M[2][1]=Y1Y2
            M[0][2]=M[2][0]=Y3Y2
            
            B=np.zeros((3))
            B[0],B[1],B[2]=XY3,XY1,XY2
            alpha3,alpha1,alpha2=np.linalg.solve(M,B)
            
    
        
            
            for i in np.ndindex(tuple(N)):
                Eps_field[i]=Eps_field[i]+np.dot(alpha3*chi3[i]+alpha1*chi1[i]+alpha2*chi2[i],X[i]) 
            
        
            Erreur=np.linalg.norm(X)/Erreur0
            X=X-alpha3*Y3-alpha1*Y1-alpha2*Y2
            
                
            print("Erreur ",Erreur)
            print("ALPHA1", alpha1)
            print("ALPHA2", alpha2)
            print("ALPHA3", alpha3)
            it+=1
            
            if j==0:
                CA11=np.zeros((d,))
                for i in np.ndindex(tuple(N)):
                    CA11+=np.dot(ki[i],Eps_field[i])/(N[0]**d)
                CA.append(CA11[0])
                ACA11=np.zeros((d,))
                for i in np.ndindex(tuple(N)):
                    ACA11+=np.dot(np.transpose(Eps_field[i]),np.dot(ki[i],Eps_field[i]))/(N[0]**d)
                ACA.append(ACA11[0])
                ER.append(Erreur)
                
        for i in np.ndindex(tuple(N)+(d,)):
            A[i+(j,)]=Eps_field[i]        #on remplit les colonnes de la matrice de localisation

    krefCA=np.zeros((d,d))
    for i in np.ndindex(tuple(N)):
        krefCA+=np.dot(ki[i],A[i])/(N[0]**d)              
                
    print("khomCA ",krefCA)
    
    krefACA=np.zeros((d,d))
    for i in np.ndindex(tuple(N)):
        krefACA+=np.dot(np.transpose(A[i]),np.dot(ki[i],A[i]))/(N[0]**d)

    
    print("khomACA ",krefACA)
    
    E_moy=np.zeros((d,d))
    for i in np.ndindex(tuple(N)):
        E_moy+=A[i]/(N[0]**d)
    
    print("Emoy ",E_moy)
    
    deltaCA=(CA[len(CA)-1]-CA[len(CA)-2])/CA[len(CA)-2]
    print("deltaCA ",deltaCA)
    deltaACA=(ACA[len(ACA)-1]-ACA[len(ACA)-2])/ACA[len(ACA)-2]
    print("deltaACA ",deltaACA)
    
    
    
    figures(Eps_field[:,:,1], "e22")

    
    return CA,ACA,ER,krefCA[0][0],krefACA[0][0]