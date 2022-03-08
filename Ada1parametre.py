import numpy as np
from greenOperatorConductionNumba import initializeGreen,operate_field
from simplemicrostructuretri import variables0
from micro4inclusions import variablesquad
from figures import figures
from math import sqrt


#La fonction ada1 permet de calculer la conductivité homognisée (khom) avec Moulinec et Suquet adaptatif 1 parametre ou Eyre Milton adaptatif 1 paramètre ou les algoritmes MS et EM basiques
#N0:nb pixels de l'image (carrée)
#k1=conductivité inclusion 1
#k2=conductivité inclusion 2  (uniquement utile pour la microstructure 1 avec 4 inclusions)
#K0 pour le choix du k0
#Algo le choix de l'algo: 0 pour conditionnement MS  ou 1 pour EM
# AME=0 pour les algoritmes basiques,AME=1 pour les algoritmes améliorés, AME=2 ppur les algoritmes améliorés bis, 3 pour MS ada ter, AME=4 pour le GC (marche avec MS)
#Prec est le nombre d'itérations
#Ini=0 pour Voigt et 1 pour initialisation de Reuss

def ada1(N0,Micro,k1,k2,K0,Algo,AME,Prec,Ini):
    
    CA,ACA,ER,ALPHA=[],[],[],[]

    ka=(1/6)*(sqrt(25*k1**2-14*k1+25)-5*k1+5)    #utile pour la microstructure simple, pour avoir khom=1
    
    if Micro==0:
        U=variables0(N0,ka,k1,1)                             #initialisation de la microstructure basique
    if Micro==1:
        U=variablesquad(N0,k1,k2)                               #initialisation de la microstructure avec 4 inclusions
    
    N=U[0]
    ki=U[1]
    d=len(N)
    filter_level=2
    k0=K0*np.eye(d)
    
    kmoyinv=0*np.eye(d)
    for i in np.ndindex(tuple(N)):
        kmoyinv+=np.linalg.inv(ki[i])/(N[0]**d)     #<K^-1>, utile pour initialisation de Reuss
    print("kmoyinv= ",kmoyinv)
    
    
    alpha=np.zeros(tuple(N)+(d,d))  #variable stockant le champ identité (conditionnement MS) ou 2*(C+C0)^(-1).C0 (conditionnement EM)
    for i in np.ndindex(tuple(N)):
        if Algo==0:
            alpha[i]=np.eye(d)
        if Algo==1:
            alpha[i]=2*np.dot(np.linalg.inv(ki[i]+k0), k0)

    A=np.zeros(tuple(N)+(d,d))   #matrice de localisation, tableau de dimension N[0]*2 de stockage des gradients de température (vecteurs), pour les deux cas de chargements [1,0] et [0,1]
            
    for j in range(d):                             #on itère sur les 2 directions
        
        it=0
        
        Eps_field = np.zeros(tuple(N)+(d,))          #initialisation du champ de déformation
        Eps_field_V = np.zeros(tuple(N)+(d,))        #initialisation du champ de déformation Voigt
        Eps_field_R = np.zeros(tuple(N)+(d,))        #initialisation du champ de déformation Reuss
        E_field = np.zeros(tuple(N)+(d,))            #initialisation du champ qui contient le chargement
        X_V= np.zeros(tuple(N)+(d,))
        X_R= np.zeros(tuple(N)+(d,))
        X= np.zeros(tuple(N)+(d,))
        Y= np.zeros(tuple(N)+(d,))
        
        for i in np.ndindex(tuple(N)):                     #Champ E (chargement)
                E_field[i+(j,)]=1.
                
        Eps_field_V=np.copy(E_field)                         #initialisation Voigt: Eps=chargement
        for i in np.ndindex(tuple(N)):
            X_V[i]=np.dot(ki[i]-k0, E_field[i])             #calcul de tau (début de la récursion)
        for i in np.ndindex(tuple(N)):                      #initialisation Reuss
            Eps_field_R[i]=np.dot(np.dot(np.linalg.inv(ki[i]),np.linalg.inv(kmoyinv)),E_field[i])
        for i in np.ndindex(tuple(N)):
            X_R[i]=np.dot(ki[i]-k0, Eps_field_R[i]) 
            
        field,field_fourier,fft,ifft,tupleK,frequencies,filters,tupleFilters = initializeGreen(N,filter_level=filter_level)
        X_V=operate_field(X_V,field_fourier,fft,ifft,tupleK,N,frequencies,filter_level,filters,tupleFilters,k0)   #calcul de - la convolée de tau par gamma0
        field,field_fourier,fft,ifft,tupleK,frequencies,filters,tupleFilters = initializeGreen(N,filter_level=filter_level)
        X_R=operate_field(X_R,field_fourier,fft,ifft,tupleK,N,frequencies,filter_level,filters,tupleFilters,k0)   #calcul de - la convolée de tau par gamma0
        X_R=E_field-Eps_field_R+X_R
        
        if Ini==0:
            X=np.copy(X_V)
            Eps_field=np.copy(Eps_field_V)
        if Ini==1:
            X=np.copy(X_R)
            Eps_field=np.copy(Eps_field_R)
        if Ini==2:
            l_V=0
            l_R=0
            X_RR,X_RV,X_X=0,0,0
            for i in np.ndindex(tuple(N)):
                X_RR+=np.dot(X_R[i],X_R[i])
                X_RV+=np.dot(X_R[i],X_V[i])
                X_X+=np.dot(X_V[i]-X_R[i],X_V[i]-X_R[i])
            l_V=(X_RR-X_RV)/X_X
            l_R=1-l_V
            print("l_V, l_R :", l_V,l_R)
            Eps_field=l_V*Eps_field_V+l_R*Eps_field_R
            X=l_V*X_V+l_R*X_R
            
        
        if AME==4:
            P=np.copy(X)
            

        Erreur0=np.linalg.norm(X)                #Erreur en norme L2 de X (itération 0, pour normaliser)
        
        while it<Prec:                                     
            print("it: ",it)            
            
            if AME==4:
                for i in np.ndindex(tuple(N)):
                    Y[i]=np.dot(ki[i]-k0, P[i])
            else:
                aX= np.zeros(tuple(N)+(d,))
                for i in np.ndindex(tuple(N)):                                                            
                    aX[i]=np.dot(alpha[i], X[i])
                for i in np.ndindex(tuple(N)):
                    Y[i]=np.dot(ki[i]-k0, aX[i])
            
            field,field_fourier,fft,ifft,tupleK,frequencies,filters,tupleFilters = initializeGreen(N,filter_level=filter_level)
            Y=operate_field(Y,field_fourier,fft,ifft,tupleK,N,frequencies,filter_level,filters,tupleFilters,k0)   
            
            if AME==4: 
                Y=-Y+P
            else:
                Y=-Y+aX
            
            
            if AME==1:                #amelioration: calcul du alpha
                E1=0
                E2=0
                
                for i in np.ndindex(tuple(N)):                                                             
                    E1+=(np.dot(X[i], Y[i]))
                    E2+=(np.dot(Y[i], Y[i]))
                al=E1/E2
            
            if AME==2:                #amelioration: calcul du alpha bis
                E1=0
                E2=0
                
                for i in np.ndindex(tuple(N)):                                                             
                    E1+=(np.dot(X[i], np.dot(ki[i],X[i])))
                    E2+=(np.dot(Y[i], Y[i]))
                al=E1/(E2*K0)
            
            if AME==3:                #amelioration: calcul du alpha bis
                E1=0
                E2=0    
                for i in np.ndindex(tuple(N)):                                                             
                    E1+=(np.dot(X[i],X[i]))
                    E2+=(np.dot(X[i], np.dot(ki[i],X[i])))
                al=K0*E1/E2
            
            if AME==4:                #amelioration: calcul du alpha GC
                E1=0
                E2=0
                
                for i in np.ndindex(tuple(N)):                                                             
                    E1+=(np.dot(X[i], X[i]))
                    E2+=(np.dot(P[i], Y[i]))
                al=E1/E2
                
            if AME==0:
                al=1
                
            
            print("ALPHA", al)
            ALPHA.append(al)
            if AME==4:
                Eps_field=Eps_field+al*P
            else:
                Eps_field=Eps_field+al*aX
            
            
            Erreur=np.linalg.norm(X)/Erreur0
            if AME==4:   
                prevX=np.copy(X)          #sauvegarde du Xn pour le GC
            X=X-al*Y
            
            if AME==4:                #amelioration: calcul du beta GC
                E1=0
                E2=0
                
                for i in np.ndindex(tuple(N)):                                                             
                    E1+=(np.dot(X[i], X[i]))
                    E2+=(np.dot(prevX[i], prevX[i]))
                beta=E1/E2
                P=X+beta*P
            
            it+=1
            
            if j==0:                          #on calcule les grandeurs d interet (CA, tACA) uniquement pour un cas de chargement (dans une direction, j=0, E=[1,0]) et on met les valeurs dans les listes dédiées
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
                
    print("krefCA ",krefCA)
    
    krefACA=np.zeros((d,d))
    for i in np.ndindex(tuple(N)):
        krefACA+=np.dot(np.transpose(A[i]),np.dot(ki[i],A[i]))/(N[0]**d)


    print("krefACA ",krefACA)
    
    E_moy=np.zeros((d,d))
    for i in np.ndindex(tuple(N)):
        E_moy+=A[i]/(N[0]**d)
    
    print("Emoy ",E_moy)
    
    deltaCA=(CA[len(CA)-1]-CA[len(CA)-2])/CA[len(CA)-2]
    print("deltaCA ",deltaCA)
    deltaACA=(ACA[len(ACA)-1]-ACA[len(ACA)-2])/ACA[len(ACA)-2]
    print("deltaACA ",deltaACA)
    legend="Quel algo?"
    if Algo==0 and AME==0:
        legend="MS"
    if Algo==1 and AME==0:
        legend="EM"
    if Algo==0 and AME==1:
        legend="MSada1"
    if Algo==1 and AME==1:
        legend="EMada1"
    if Algo==0 and AME==2:
        legend="MSada1bis"
    if Algo==0 and AME==3:
        legend="MSada1ter"
    if Algo==0 and AME==4:
        legend="MS GC"
    
    if Ini==0:
        initial="Ini Voigt"
    if Ini==1:
        initial="Ini Reuss"
    if Ini==2:
        initial="Ini Moyenne"
        
    figures(Eps_field[:,:,1], "e22 "+legend+"  k0="+str(K0)+" "+initial) 
    
    return CA,ACA,ER,ALPHA,krefCA[0][0],krefACA[0][0]