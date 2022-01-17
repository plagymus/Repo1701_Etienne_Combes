import numpy as np
from greenOperatorConductionNumba import initializeGreen,operate_field
from simplemicrostructuretri import variables0
from micro4inclusions import variablesquad
from figures import figures



#La fonction ada1 permet de calculer la conductivité homognisée (khom) avec Moulinec et Suquet adaptatif 1 parametre ou Eyre Milton adaptatif 1 paramètre
#N0:nb pixels de l'image (carrée)
#k1=conductivité matrice
#k2=conductivité inclusion
#Prec est le nombre d'itération ou la precision sur l'erreur
#K0 pour le choix du k0
#Algo le choix de l'algo: 0 pour conditionnement MS  ou 1 pour EM

def ada1(N0,Micro,k1,k2,K0,Algo,AME,Prec):
    
    CA,ACA,ER,ALPHA=[],[],[],[]

    khom=k1*((k1+k2)-(k1-k2)*0.25)/((k1+k2)+(k1-k2)*0.25)   #khom
    
    if Micro==0:
        U=variables0(N0,k1,k2,khom)                             #initialisation de la microstructure
    if Micro==1:
        U=variablesquad(N0,k1,k2)                               #initialisation de la microstructure
    
    N=U[0]
    ki=U[1]
    d=len(N)
    filter_level=2
    k0=K0*np.eye(d)
    
    alpha=np.zeros(tuple(N)+(d,d))  #variable stockant le champ identité (conditionnement MS) ou 2*(C+C0)^(-1).C0 (conditionnement EM)
    for i in np.ndindex(tuple(N)):
        if Algo==0:
            alpha[i]=np.eye(d)
        if Algo==1:
            alpha[i]=2*np.dot(np.linalg.inv(ki[i]+k0), k0)

    
    A=np.zeros(tuple(N)+(d,d))   #matrice de localisation, tableau de dimension N[0]*2 de stockage des gradients de température (vecteurs), pour les deux cas de chargements [1,0] et [0,1]
            
    for j in range(d):                             #on itère sur les 2 directions
        
        it=0

        relativeChange = 1
        
        Eps_field = np.zeros(tuple(N)+(d,))        #initialisation du champ de déformation
        E_field = np.zeros(tuple(N)+(d,))          #initialisation du champ qui contient le chargement
        X= np.zeros(tuple(N)+(d,))
        Y= np.zeros(tuple(N)+(d,))
        
        for i in np.ndindex(tuple(N)):                     #Champ E (chargement)
                E_field[i+(j,)]=1.
        Eps_field=np.copy(E_field)                         #initialisation: Eps=chargement
        
        for i in np.ndindex(tuple(N)):
                X[i]=np.dot(ki[i]-k0, Eps_field[i])             #calcul de tau (début de la récursion)
            
        field,field_fourier,fft,ifft,tupleK,frequencies,filters,tupleFilters = initializeGreen(N,filter_level=filter_level)
        X=operate_field(X,field_fourier,fft,ifft,tupleK,N,frequencies,filter_level,filters,tupleFilters,k0)   #calcul de - la convolée de tau par gamma0
        

        Erreur0=np.linalg.norm(X)                #Erreur en norme L2 de X (itération 0, pour normaliser)
        
        while it<Prec:                                     
            print("it: ",it,"      RelativeChange: ",relativeChange)            
            
            aX= np.zeros(tuple(N)+(d,))
            for i in np.ndindex(tuple(N)):                                                            
                aX[i]=np.dot(alpha[i], X[i])
            
            for i in np.ndindex(tuple(N)):
                Y[i]=np.dot(ki[i]-k0, aX[i])
            
            field,field_fourier,fft,ifft,tupleK,frequencies,filters,tupleFilters = initializeGreen(N,filter_level=filter_level)
            Y=operate_field(Y,field_fourier,fft,ifft,tupleK,N,frequencies,filter_level,filters,tupleFilters,k0)   

            Y=-Y+aX
            
            
            if AME==1:                #amelioration calcul du alpha ou algo basiques
                E1=0
                E2=0
                
                for i in np.ndindex(tuple(N)):                                                             
                    E1+=(np.dot(X[i], Y[i]))
                    E2+=(np.dot(Y[i], Y[i]))
                al=E1/E2
            else:
                al=1
                
            
            print("ALPHA", al)
            ALPHA.append(al)
            
            Eps_field=Eps_field+al*aX
            
            Erreur=np.linalg.norm(X)/Erreur0
            
            X=X-al*Y
                
            
            it+=1
            
            if j==0:                          #on calcule les grandeurs d interet (CA, tACA) uniquement pour un cas de chargement (dans une direction) et on met les valeurs dans les listes dédiées
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
    
    if Algo==0:
        legend="MSada1"
    else:
        legend="EMada1"
    
    figures(Eps_field[:,:,1], "e22 "+legend+"  k0="+str(K0)) 
    
    return CA,ACA,ER,krefCA[0][0],krefACA[0][0]