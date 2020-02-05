# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 14:25:16 2019

@author: daniel
"""

import numpy as np
import scipy.linalg as sc
import matplotlib.pyplot as plt
import pandas as pd
from pandas import ExcelWriter
from numba import jit
from mpl_toolkits.mplot3d import Axes3D 

###############################################################################
################ MONTAGEM DAS MATRIZES DE RIGIDEZE MASSA ######################
###############################################################################

def Matriz(Nome,Nn):
    ## Nome : Nome do arquivo(lembrar de colocar ".xlsx").
    ## Nn : Número de nós.
    Arquivo = pd.read_excel(Nome)
    cx = list(Arquivo['Cx'])[0:Nn]
    cy = list(Arquivo['Cy'])[0:Nn]
    

    Id1 = list(Arquivo['barra (nó 1)'])
    Id2 = list(Arquivo['barra (nó 2)'])
    Nb = len(Id1)
                                ## Importando arquivos Excel

    ID =np.zeros((2,Nb))
    ID[0,:] = Id1
    ID[1,:] = Id2

    A = list(Arquivo['Area(m2)'])
    I =list(Arquivo['Inércia(m4)'])
    RHO = list(Arquivo['Densidade'])

    E = 28*10**9##[N/m2]


 
    ## Matriz identidade - identificar a conexão entre os nós
    lb = ID[0,:]
    l1 = lb*3-2
    l2 =lb*3-1
    l3 =lb*3
    lb2 =ID[1,:]      ##Montagem da matriz Id em relação aos graus de liberdade
    l4 = 3*lb2-2
    l5= 3*lb2 -1
    l6= 3*lb2
    IDG = np.zeros((6,70))
    IDG[0,:]=l1
    IDG[1,:]=l2
    IDG[2,:]=l3
    IDG[3,:]=l4
    IDG[4,:]=l5
    IDG[5,:]=l6
        
    
## Comprimento de cada barra 


    Lx = np.zeros(Nb)
    Ly = np.zeros(Nb)
    cosx = np.zeros(Nb)
    cosy = np.zeros(Nb)

    L = np.zeros(Nb)

    for n in range (Nb):
    
        k1 = int(ID[0,n] -1)  ## BUSCANDO OS NÓS A PARTIR DA MATRIZ ID
        k2 = int(ID[1,n] -1)
        Lx[n] = cx[k2] - cx[k1]
        Ly[n] = cy[k2] - cy[k1]
        L[n] = np.sqrt(Lx[n]**2 + Ly[n]**2)
        cosx[n] = Lx[n]/L[n]
        cosy[n] = Ly[n]/L[n]
    
        
   ## Matriz de rigidez global.

    K = np.zeros((132,132))
    M = np.zeros((132,132))

    for i in range (Nb):
    
    ## Matriz de rigidez local da barra
        Ke =np.array([[E*A[i]/L[i], 0, 0, -E*A[i]/L[i],0 ,0 ],
                  [0, 12*E*I[i]/(L[i]**3), 6*E*I[i]/(L[i]**2), 0,
                   -12*E*I[i]/(L[i]**3),6*E*I[i]/(L[i]**2)],
                  [0,6*E*I[i]/(L[i]**2), 4*E*I[i]/L[i], 0, 
                   -6*E*I[i]/(L[i]**2), 2*E*I[i]/L[i] ],
                  [-E*A[i]/L[i], 0, 0, E*A[i]/L[i],0 ,0 ],
                  [0, -12*E*I[i]/(L[i]**3), -6*E*I[i]/(L[i]**2),
                   0,12*E*I[i]/(L[i]**3),-6*E*I[i]/(L[i]**2)],
                  [0,6*E*I[i]/(L[i]**2), 2*E*I[i]/L[i], 0,
                   -6*E*I[i]/(L[i]**2), 4*E*I[i]/L[i] ]])

    ## Matriz de massa local da barra
        Me = ((RHO[i]*A[i]*L[i])/420)*np.array([[140, 0, 0, 70, 0, 0],
         [0, 156, 22*L[i], 0, 54, -13*L[i]],
         [0, 22*L[i], 4*(L[i]**2), 0, 13*L[i], -3*(L[i]**2)],
         [70, 0, 0, 140, 0, 0],
         [0, 54, 13*L[i], 0, 156, -22*L[i]],
         [0, -13*L[i], -3*(L[i]**2), 0, -22*L[i], 4*(L[i]**2)]])
   
    
        R = np.array([[cosx[i], cosy[i], 0, 0 ,0 ,0],
                  [-cosy[i], cosx[i],0, 0, 0, 0],
                  [0,0,1,0,0,0],                     ## Matriz de rotação
                  [0,0,0,cosx[i], cosy[i], 0],
                  [0, 0, 0,-cosy[i], cosx[i],0],
                  [0,0,0,0,0,1]])
                     
        KT = np.dot(np.dot(R.T, Ke),R)            ## Rotação das matrizes 
        MT = np.dot(np.dot(R.T, Me),R)
    
        k_temp1 = np.zeros((132,132))
    
                                      ## Matrizes temporárias para alocação
    
        m_temp1 = np.zeros((132,132))
   
   
    
        m= int(IDG[0,i]-1)
        n= int(IDG[2,i])
        o= int(IDG[3,i]-1)
        p= int(IDG[5,i])
    
        k_temp1[m:n,m:n] = KT[0:3,0:3]
        k_temp1[o:p,m:n] = KT[3:6,0:3]
        k_temp1[m:n,o:p] = KT[0:3,3:6]
        k_temp1[o:p,o:p] = KT[3:6,3:6]
    
        K += k_temp1 
    
        m_temp1[m:n,m:n] = MT[0:3,0:3]
        m_temp1[o:p,m:n] = MT[3:6,0:3]
        m_temp1[m:n,o:p] = MT[0:3,3:6]
        m_temp1[o:p,o:p] = MT[3:6,3:6]
    
        M += m_temp1 
        
    return K,M
###############################################################################
############### Eliminação dos graus de liberdade restritos####################
###############################################################################    
def Restr(K,M,Nr): 
    
    ## K : Matriz de rigidez
    ## M : Matriz de massa
    ## Nr : Lista com os graus de liberdade restritos
    Kr_1 = np.delete(K,Nr,0)
    Kr  = np.delete(Kr_1,Nr,1)
    Mr_1 = np.delete(M,Nr,0)
    Mr  = np.delete(Mr_1,Nr,1)

    df = pd.DataFrame(Kr)
    writer = ExcelWriter('Matriz de rigidez.xlsx')
    df.to_excel(writer,'Sheet1', index=False) 
    writer.save()
                               ##Exportação das matrizes de massa e rigidez
    df1 = pd.DataFrame(Mr)
    writer = ExcelWriter('Matriz de massa.xlsx')
    df1.to_excel(writer,'Sheet1', index=False)
    writer.save()
    
    return Kr,Mr
###############################################################################
########### Cálculo das frequências naturais e modos de vibração###############
###############################################################################
def Eig(Kr,Mr,N,d,Na):
    ## Kr: Matriz de rigidez restringida
    ## Mr : Matriz de rigidez restringida
    ## N : Número de modos desejados
    ## d : distancia entre pilares
    ## Na : número de andares do pórtico
    w21,Phi1 = sc.eig(Kr,Mr)

    iw = w21.argsort()
    w21 = w21[iw]          ## Garantindo a ordem dos autovalores e autovetores
    Phi1 = Phi1[:,iw]

    wr = np.real(w21)
    wk = np.sqrt(wr)
    fk = np.real(wk/(2*np.pi))
    wk = 2*np.pi*fk

    for k in range(N):
        print(k+1, "ª frequencia natural = {0:3.2f}Hz".format(fk[k]),"\n")
    
    plt.figure(1, figsize=(8,8))
    x = np.arange(Na+1)
    Phi = Phi1[::12,0:N]


    for k in range(N):
        pk = np.zeros(Na+1)
        pk[1:] = Phi[:,k]
        pk /=np.max(np.abs(pk))
        plt.subplot(1,N,k+1)
    
        for n in range(Na):
        
            o = np.linspace(pk[n+1],pk[n+1]+3*d,10)
            y1 = np.ones(Na)*n+1          ##Criação dos andares horizontais
        
            plt.plot(o, y1, 'b')
    
        plt.plot(pk[1:],x[1:],'bo')
        plt.plot(pk[1:]+d, x[1:], 'bo')
        plt.plot(pk[1:]+2*d, x[1:], 'bo')   ## Plotagem das linhs e dos nós
        plt.plot(pk[1:]+3*d, x[1:], 'bo')
        plt.plot(pk,x,'b')
        plt.plot(pk+d, x,'b')
        plt.plot(pk+2*d, x,'b')
        plt.plot(pk+3*d, x,'b')
    
    
    
        plt.xlim(-2, 12);  plt.ylabel(str(k+1));
        plt.ylim( 0.0, 11);
        plt.title('f= {0:3.2f}Hz'.format(fk[k]));
        plt.grid(True)
    
    return fk,wk
###############################################################################
######################### Matriz de amortecimento##############################
###############################################################################
def Rayleigh(Kr,Mr,wk,z1,z2): 
    ## Kr: Matriz de rigidez restringida
    ## Mr : Matriz de rigidez restringida
    ## wk : Frequencias naturais(rad/s)
    ## z1 : Fator de amortecimento do primeiro modo 
    ## z2 : Fator de amortecimento do segundo modo
    zeta =np.zeros(2)
    zeta[0] = z1
    zeta[1] = z2
    a0 = -2*(wk[0]*wk[1])/(wk[0]**2 - wk[1]**2)
    a1 = a0*( zeta[0]*wk[1] - zeta[1]*wk[0])
    a2 = a0*(-zeta[0]/wk[1] + zeta[1]/wk[0])      


    Cr  = a1*Mr + a2*Kr              



    df2 = pd.DataFrame(Cr)
    writer = ExcelWriter('Matriz de Amortecimento.xlsx')
    df2.to_excel(writer,'Sheet1', index=False)
    writer.save()
    
    return Cr
###############################################################################
#########################Espectro de Kanai-Tajimi##############################
###############################################################################

def Kanai_Tajimi(Ap,tipo,duraçao,dt):
    ## Ap : Peek Ground Acceleration
    ## tipo : tipo de solo 
    ## duraçao : tempo de duração do sinal
    ## dt : Passo temporal
    g = 9.806
    Ap *= g
    pg = 3
    tf = int(duraçao/dt)
    if tipo == 'rocha':
        wg = 8* np.pi
        zg = 0.6
    elif tipo == 'solo_rígido':
        wg = 5* np.pi
        zg = 0.6
    elif tipo == 'solo_mole':
         wg = 2.4* np.pi
         zg = 0.85
        
    f = np.linspace(0,25,tf)
    df = f[1]-f[0]
    w = 2*np.pi*f
    S0 = (Ap**2)/((pg**2)*(np.pi*wg*((1/(2*zg))+2*zg)))
    Sg = S0*((1+4*(zg**2)*(w/wg)**2)/(((1-(w/wg)**2)**2)+4*(zg**2)*(w/wg)**2))


    plt.figure(2, figsize=(8,4)) 
    plt.plot(f,Sg,'b')
    plt.xlabel('frequência(Hz)'); plt.ylabel('Densidade espectral(m²/s³)');
    plt.xlim(0,20); plt.ylim(0,max(Sg)*2);plt.title(' Espectro de aceleração')
    plt.grid(True)

    import random 
    P = np.zeros(tf)
    for i in range(tf):
        P1 = random.uniform(0,2*np.pi)
        P[i] = P1
   
    ag= np.zeros(tf)
    t = np.linspace(0,duraçao,tf)
    S = np.zeros(tf)



    for i in range(tf):

        S =np.sqrt(2*Sg*df)*np.cos(w*t[i]+ P)
        ag[i] = sum(S)

    ag*= Ap/np.max(abs(ag))  ## Normalização das acelerações
    
## Função de envoltória

    env =np.ones(tf)
    env1 = np.arange(int(0.05*tf))/(0.05*tf)
    env2 = (1.148698355**t[0:int(0.8*tf)])/64/4
    env[0:int(0.05*tf)] = env1
    env[int(0.2*tf):tf] = env2[::-1]

    plt.figure(3,figsize=(8,4))

    plt.plot(t,ag,'b')
    plt.xlabel('Tempo (s)'); plt.ylabel('Aceleração (m/s²)');
    plt.xlim(0,duraçao); plt.ylim(-Ap,Ap);plt.title(' Aceleração do solo')
    plt.grid(True)

    plt.figure(4,figsize=(8,4))
    plt.plot(t,ag,'c')
    plt.plot(t,env,'r--',t,-env,'r--')
    plt.xlabel('Tempo (s)'); plt.ylabel('Aceleração (m/s²)');
    plt.xlim(0,duraçao); plt.ylim(-Ap,Ap);plt.title(' Função de envoltória')
    plt.grid(True)

    age = ag*env
    age*= Ap/np.max(abs(age))
    plt.figure(5,figsize=(8,4))
    plt.plot(t,age,'b')
    plt.xlabel('Tempo (s)'); plt.ylabel('Aceleração (m/s²)');
    plt.xlim(0,duraçao); plt.ylim(-Ap,Ap);
    plt.title('Aceleração do solo parametrizada')
    plt.grid(True)
    
    return t,age

###############################################################################
####################### Criação do vetor de forças#############################
###############################################################################
def Sismo(Mr,age,t):
    ## Mr : Matriz de rigidez restrita 
    ## t : lista do tempo discretizado 
    ## age : sinal de aceleração discetizado 
    
    n = int (len(Mr[0,:]))
    B = np.zeros((n,1))
    B[::3,0] = np.ones(int(n/3))
    ag = np.zeros((1,len(t)))
    ag[0,:] = age
    F1 = np.dot(B.T,Mr)
    F = np.dot(F1.T,ag) 
    
    plt.figure(6,figsize=(8,4))   
    plt.plot(t,F[n-3,:],'b')
    plt.xlabel('Tempo (s)'); plt.ylabel('Força N')
    plt.xlim(0,max(t)); plt.ylim(-max(F[n-3])*1.2,max(F[n-3])*1.2);
    plt.title('Força no 10º pavimento')
    plt.grid(True)
    
    return F
###############################################################################   
######################### MÉTODO DE NEWMARK ###################################
###############################################################################

def Newmark(Kr,Mr,Cr,F,t):
     ## Kr: Matriz de rigidez restringida
     ## Mr : Matriz de rigidez restringida
     ## Cr : Matriz de amortecimento
     ## F : Vetor de força discretizado no tempo
     ## t : lista do tempo discretizado
    tf = int(len(t))
    n = len(F[:,0])
    A = np.zeros((n,tf))
    v = np.zeros((n,tf))
    d = np.zeros((n,tf))
    dt =t[1]-t[0]

    delta = 0.5
    alfa = 0.25
    a0 = 1/(alfa*(dt**2))
    a1 = 1/(alfa*dt)
    a2 = (1/(2*alfa))-1
    a3 = delta/(dt*alfa)
    a4 = delta/alfa - 1
    a5 = (dt/2)*(delta/alfa - 2)

    A[:,0] = np.dot(np.linalg.inv(Mr),(F[:,0]-np.dot(Cr,
     v[:,0])-np.dot(Kr,d[:,0])))
    d4 = a0*Mr + a3*Cr + Kr
    D = np.linalg.inv(d4)

    for i in range(tf-1):
        d1 = np.dot(Mr,(a0*d[:,i]+ a1*v[:,i] + a2*A[:,i]))
        d2 = np.dot(Cr,(a3*d[:,i]+ a4*v[:,i] + a5*A[:,i]))
        d3 = F[:,i+1]+ d1 + d2
        d[:,i+1] = np.dot(D,d3)
        v[:,i+1] = a3*(d[:,i+1] - d[:,i]) - a4*v[:,i] - a5*A[:,i]
        A[:,i+1] = a0*(d[:,i+1] - d[:,i]) - a1*v[:,i] - a2*A[:,i]
    
    return d,v,A

###############################################################################
#################### Captação de sinal real.csv ###############################
###############################################################################
def Signal(Nome_Arquivo):
    # Nome_Arquivo : nome do arquivo que será importado ( lembrar do ".csv")
   teste = pd.read_csv(Nome_Arquivo,delimiter = ',')
   t = list(teste['Time (s)'])

   X =[]
   X.append(list(teste['X-Axis (m/s2)']))
   X.append(list(teste['Y- Axis (m/s2)']))
   X.append(list(teste['Z-Axis (m/s2)']))

   plt.figure(200,figsize=(8,28))
   for n in range(3):
        plt.figure(figsize=(8,12))
        plt.subplot(3,1,n+1)
        plt.plot(t,X[n],'green')
        plt.title('aceleração')
        plt.xlim(0,50); plt.ylim(-0.4,0.4)
        plt.xlabel('Tempo (s)'); plt.ylabel('Aceleração (g)');
        plt.grid(True)
      
       
   return t,X

##############################################################################


###############################################################################
################ MONTAGEM DAS MATRIZES DE RIGIDEZE MASSA 3d####################
###############################################################################

def Matriz3D(Nome,Nn):
    ## Nome : Nome do arquivo(lembrar de colocar ".xlsx").
    ## Nn : Número de nós.
    Arquivo = pd.read_excel(Nome)

    cx = list(Arquivo['Cx'])[0:Nn]
    cy = list(Arquivo['Cy'])[0:Nn]
    cz = list(Arquivo['Cz'])[0:Nn]
    
    Id1 = list(Arquivo['barra (nó 1)'])
    Id2 = list(Arquivo['barra (nó 2)'])
    Nb = len(Id1)
  ## Importando arquivos Excel

    ID =np.zeros((2,Nb))
    ID[0,:] = Id1
    ID[1,:] = Id2

    A = list(Arquivo['Area(m2)'])
    Iz =list(Arquivo['Inércia z (m4)'])
    Iy =list(Arquivo['Inércia y (m4)'])
    RHO = list(Arquivo['Densidade'])

    E = 28*10**9##[N/m2]
    G = 1630*10**6
    J = 2*Iz
## Matriz identidade - identificar a conexão entre os nós
    lb = ID[0,:]
    l1 = lb*6-5
    l2 =lb*6-4
    l3 =lb*6-3
    l4 = lb*6-2
    l5 =lb*6-1
    l6 =lb*6
    lb2 =ID[1,:]          ##Montagem da matriz Id em relação aos graus de liberdade
    l7 = lb2*6-5
    l8 =lb2*6-4
    l9 =lb2*6-3
    l10 = lb2*6-2
    l11 =lb2*6-1
    l12 =lb2*6

    IDG = np.zeros((12,Nb))
    IDG[0,:]=l1
    IDG[1,:]=l2
    IDG[2,:]=l3
    IDG[3,:]=l4
    IDG[4,:]=l5
    IDG[5,:]=l6
    IDG[6,:]=l7
    IDG[7,:]=l8
    IDG[8,:]=l9
    IDG[9,:]=l10
    IDG[10,:]=l11
    IDG[11,:]=l12


    K = np.zeros((Nn*6,Nn*6))
    M = np.zeros((Nn*6,Nn*6))

    for i in range (Nb):
    
        k1 = int(ID[0,i] -1)  ## BUSCANDO OS NÓS A PARTIR DA MATRIZ ID
        k2 = int(ID[1,i] -1)
        Lx = cx[k2] - cx[k1]
        Ly = cy[k2] - cy[k1]
        Lz = cz[k2] - cz[k1]
        L = np.sqrt(Lx**2 + Ly**2 +Lz**2)
        l = Lx/L
        m = Ly/L
        n = Lz/L
        D = np.sqrt(l**2 + m**2)
        
        if D == 0:
            R = np.array([[0,0,1,0,0,0,0,0,0,0,0,0],
                  [0,1,0,0,0,0,0,0,0,0,0,0],
                  [1,0,0,0,0,0,0,0,0,0,0,0],
                  [0,0,0,0,0,1,0,0,0,0,0,0],
                  [0,0,0,0,1,0,0,0,0,0,0,0],
                  [0,0,0,-1,0,0,0,0,0,0,0,0],
                  [0,0,0,0,0,0,0,0,1,0,0,0],
                  [0,0,0,0,0,0,0,1,0,0,0,0],
                  [0,0,0,0,0,0,-1,0,0,0,0,0],
                  [0,0,0,0,0,0,0,0,0,0,0,1],
                  [0,0,0,0,0,0,0,0,0,0,1,0],
                  [0,0,0,0,0,0,0,0,0,-1,0,0]])
        else:
        
            R = np.array([[l,m,n,0,0,0,0,0,0,0,0,0],
                  [(-m/D),l/D,0,0,0,0,0,0,0,0,0,0],
                  [(-l*n)/D,-m*n/D,D,0,0,0,0,0,0,0,0,0],
                  [0,0,0,l,m,n,0,0,0,0,0,0],
                  [0,0,0,-m/D,l/D,0,0,0,0,0,0,0],
                  [0,0,0,-l*n/D,-m*n/D,D,0,0,0,0,0,0],
                  [0,0,0,0,0,0,l,m,n,0,0,0],
                  [0,0,0,0,0,0,-m/D,l/D,0,0,0,0],
                  [0,0,0,0,0,0,-l*n/D,-m*n/D,D,0,0,0],
                  [0,0,0,0,0,0,0,0,0,l,m,n],
                  [0,0,0,0,0,0,0,0,0,-m/D,l/D,0],
                  [0,0,0,0,0,0,0,0,0,-l*n/D,-m*n/D,D]])
    
    
       
   
    ## Matriz de rigidez local da barra
        Ke =np.array(([[E*A[i]/L, 0, 0,0,0,0, -E*A[i]/L,0 ,0,0,0,0],
                  [0, 12*E*Iz[i]/(L**3),0,0,0, 6*E*Iz[i]/(L**2), 0,-12*E*Iz[i]/(L**3),0,0,0,6*E*Iz[i]/(L**2)],
                  [0,0,12*E*Iy[i]/(L**3),0,-6*E*Iy[i]/(L**2),0,0,0,-12*E*Iy[i]/(L**3),0,-6*E*Iy[i]/(L**2),0],
                  [0,0,0,G*J[i]/L,0,0,0,0,0,-G*J[i]/L,0,0],
                  [0,0,-6*E*Iy[i]/(L**2),0,4*E*Iy[i]/L,0,0,0,6*E*Iy[i]/(L**2),0,2*E*Iy[i]/L,0],
                  [0,6*E*Iz[i]/(L**2),0,0,0, 4*E*Iz[i]/L,0,-6*E*Iz[i]/(L**2),0,0,0, 2*E*Iz[i]/L],
                  [-E*A[i]/L, 0, 0,0,0,0, E*A[i]/L,0 ,0,0,0,0],
                  [0, -12*E*Iz[i]/(L**3),0,0,0, -6*E*Iz[i]/(L**2), 0,12*E*Iz[i]/(L**3),0,0,0,-6*E*Iz[i]/(L**2)],
                  [0,0,-12*E*Iy[i]/(L**3),0,6*E*Iy[i]/(L**2),0,0,0,12*E*Iy[i]/(L**3),0,6*E*Iy[i]/(L**2),0],
                  [0,0,0,-G*J[i]/L,0,0,0,0,0,G*J[i]/L,0,0],
                  [0,0,-6*E*Iy[i]/(L**2),0,2*E*Iy[i]/L,0,0,0,6*E*Iy[i]/(L**2),0,4*E*Iy[i]/L,0],
                  [0,6*E*Iz[i]/(L**2),0,0,0, 2*E*Iz[i]/L,0,-6*E*Iz[i]/(L**2),0,0,0, 4*E*Iz[i]/L]]))
    
    ## Matriz de massa local da barra
        rx = (0.3*L**3)/(12*A[i])
    
        L = L/2
        Me = ((RHO[i]*A[i]*L)/105)*np.array([[70, 0, 0,0,0,0,35, 0, 0,0,0,0],
         [0, 78,0,0,0, 22*L, 0, 27,0,0,0, -13*L],
         [0,0,78,0, -22*L,0,0,0, 27, 0, 13*L, 0],
         [0,0,0,70*rx, 0, 0,0,0,0,-35*rx, 0, 0],
         [0,0,-22*L,0,8*L**2,0,0,0,-13*L,0,-6*L**2,0],
         [0,22*L,0,0,0,8*L**2,0,13*L,0,0,0,-6*L**2],
         [35, 0, 0,0,0,0,70, 0, 0,0,0,0],
         [0, 27,0,0,0, 13*L, 0, 78,0,0,0, -22*L],
         [0,0,27,0, -13*L,0,0,0, 78, 0, 22*L, 0],
         [0,0,0,-35*rx, 0, 0,0,0,0,70*rx, 0, 0],
         [0,0,13*L,0,-6*L**2,0,0,0,22*L,0,8*L**2,0],
         [0,-13*L,0,0,0,-6*L**2,0,-22*L,0,0,0,8*L**2]])
   
   
    
                     
        KT = np.dot(np.dot(R.T, Ke),R)            ## Rotação das matrizes 
        MT = np.dot(np.dot(R.T, Me),R)
    
        k_temp1 = np.zeros((Nn*6,Nn*6))
    
                                      ## Matrizes temporárias para alocação
    
        m_temp1 = np.zeros((Nn*6,Nn*6))
   
   
    
        j= int(IDG[0,i]-1)
        f= int(IDG[5,i])
        o= int(IDG[6,i]-1)
        p= int(IDG[11,i])
    
        k_temp1[j:f,j:f] = KT[0:6,0:6]
        k_temp1[o:p,j:f] = KT[6:12,0:6]
        k_temp1[j:f,o:p] = KT[0:6,6:12]
        k_temp1[o:p,o:p] = KT[6:12,6:12]
    
        K += k_temp1 
    
        m_temp1[j:f,j:f] = MT[0:6,0:6]
        m_temp1[o:p,j:f] = MT[6:12,0:6]
        m_temp1[j:f,o:p] = MT[0:6,6:12]
        m_temp1[o:p,o:p] = MT[6:12,6:12]
    
        M += m_temp1 
        
    return K,M
###############################################################################
########### Cálculo das frequências naturais e modos de vibração###############
###############################################################################
def Eig3D(Kr,Mr,N,d,Na):
    ## Kr: Matriz de rigidez restringida
    ## Mr : Matriz de rigidez restringida
    ## N : Número de modos desejados
    ## d : distancia entre pilares
    ## Na : número de andares do pórtico
    w21,Phi1 = sc.eig(Kr,Mr)

    iw = w21.argsort()
    w21 = w21[iw]              ## Garantindo a ordem dos autovalores e autovetores
    Phi1 = Phi1[:,iw]

    wr = np.real(w21)
    wk = np.sqrt(wr)
    fk = np.real(wk/(2*np.pi))
    wk = 2*np.pi*fk

    for k in range(N):
        print(k+1, "ª frequencia natural = {0:3.2f}Hz".format(fk[k]),"\n")


  
    fig = plt.figure(figsize=(10,15))

    ax = fig.add_subplot(111,projection = '3d')

    x = np.arange(Na+1)

    Phi = Phi1[::24,N-1]
    Phi2 = Phi[0:Na]
    Phi3 = Phi1[2:,N-1]
    Phi4 = Phi3[::24]
    Phi5 = Phi4[0:Na]
    
    pk = np.zeros(Na+1)
    pk[1:] = Phi2
    pk /=np.max(np.abs(pk))
    
    pk1 = np.zeros(Na+1)
    pk1[1:] = Phi5
    pk1 /=np.max(np.abs(pk1))
    
    
    for n in range(Na):
        
        o = np.linspace(pk[n+1],pk[n+1]+9,10)
        y1 = np.ones(Na)*n+1            ##Criação dos 10 andares horizontais
        z = np.zeros(Na)
        ax.plot(o, z+pk1[n+1],y1,'b')
        ax.plot(o, z+pk1[n+1]+d,y1,'b')
        ax.plot(o, z+pk1[n+1]+2*d,y1,'b')
    
    for n in range(Na):
        
        y1 = np.ones(Na)*n+1            ##Criação dos 10 andares horizontais
        o = np.linspace(pk1[n+1],pk1[n+1]+6,10)
        k = np.zeros(Na)
        ax.plot(k+pk[n+1], o,y1,'b')
        ax.plot(k+3+pk[n+1], o,y1,'b')
        ax.plot(k+6+pk[n+1], o,y1,'b')
        ax.plot(k+9+pk[n+1], o,y1,'b')
    
    z = pk1[1:]
    ax.plot(pk[1:], z, x[1:],'bo')
    ax.plot(pk[1:]+d, z, x[1:], 'bo')
    ax.plot(pk[1:]+2*d, z, x[1:], 'bo')   ## Plotagem das linhs e dos nós
    ax.plot(pk[1:]+3*d, z , x[1:],'bo')

    ax.plot(pk[1:],z+d,x[1:],'bo')
    ax.plot(pk[1:]+d,z+d,x[1:], 'bo')
    ax.plot(pk[1:]+2*d,z+d,x[1:], 'bo')   ## Plotagem das linhs e dos nós
    ax.plot(pk[1:]+3*d,z+d,x[1:],'bo')

    ax.plot(pk[1:],z+6,x[1:],'bo')
    ax.plot(pk[1:]+3,z+6,x[1:], 'bo')
    ax.plot(pk[1:]+2*3,z+6,x[1:], 'bo')   ## Plotagem das linhs e dos nós
    ax.plot(pk[1:]+3*3,z+6,x[1:],'bo')


    z = pk1
    ax.plot(pk,z,x,'b')
    ax.plot(pk+3, z,x,'b')
    ax.plot(pk+2*3,z,x,'b')
    ax.plot(pk+3*3,z,x,'b')


    ax.plot(pk,z+3,x,'b')
    ax.plot(pk+3, z+3,x,'b')
    ax.plot(pk+2*3,z+3,x,'b')
    ax.plot(pk+3*3,z+3,x,'b')


    ax.plot(pk,z+6,x,'b')
    ax.plot(pk+3, z+6,x,'b')
    ax.plot(pk+2*3,z+6,x,'b')
    ax.plot(pk+3*3,z+6,x,'b')
    
    
    

    plt.title('Frequência natural =1.33 Hz')
    ax.set_xlabel('Eixo X (m)')
    ax.set_ylabel('Eixo Y (m)')
    ax.set_zlabel('Andar')
    plt.grid(True)
    plt.show()
    return fk,wk
###############################################################################
###############################################################################
####################### Criação do vetor de forças#############################
###############################################################################
def Sismo3D(Mr,age,t):
    ## Mr : Matriz de rigidez restrita 
    ## t : lista do tempo discretizado 
    ## age : sinal de aceleração discetizado 
    tf = int(len(t))
    n = int (len(Mr[0,:]))
    B = np.zeros(n)
    B[::6] = np.ones(int(n/6))
    F = np.zeros((n,tf))
    F1 = np.dot(Mr,B)

    for i in range (tf):
        F[:,i] = F1* age[i] 
    
    plt.figure(6,figsize=(8,4))   
    plt.plot(t,F[n-6 ,:],'b')
    plt.xlabel('Tempo (s)'); plt.ylabel('Força N')
    plt.xlim(0,max(t)); plt.ylim(-max(F[n-6])*1.5,max(F[n-6])*1.5);
    plt.title('Força no 10º pavimento')
    plt.grid(True)
    
    return F
##############################################################################
