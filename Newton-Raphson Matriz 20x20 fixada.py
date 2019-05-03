# -*- coding: utf-8 -*-
"""
Created on Thu May  2 03:20:56 2019

@author: Bruno Felipe
"""
#É assumido que o sistema a ser resolvido é puramente linear dentro do conjunto IR^n->IR^m. Portanto já é colocado da forma J(x) (Jacobiano) a matriz convergente gerada.
#É gerado polinômios do tipo Fk(xi)=Somatório[1:n](aix^i)=bk=>Pk(xi)=Fk(xi)-bk=Somatório[1:n](aix^i)-bk=0
#Será adotado bk=(k+1)*10 (0<=k<n) para todos os polinômios, para facilitar a geração automática.
def jacobiana(Linf,Lsup):
    import random
#Baseado no critério de sessenfeld det[A]{ord n}>0
    n=20
    matQ=[]
    beta=[]
    bet=[]
    somE=[]
    for i in range(n):
        matQ.append([0]*n)
        beta.append([0]*n)
    for i in range(n):
        somE.append(0)
    for i in range(n):
        for j in range(n):
            if i!=j:
                matQ[i][j]=random.randint(Linf,Lsup)
    for i in range(n):
        a=random.uniform(0.1,0.9)
        bet.append(a)
    for i in range(n):
        for j in range(n):
            beta[j][i]=bet[i]    
    for i in range(n):
        for j in range(n):
            if i<j:
                beta[i][j]=1
            elif i==j:
                beta[i][j]=0
    for i in range(n):
        for j in range(n):
            somE[i]=somE[i]+abs(matQ[i][j]*beta[i][j])
    for i in range(n):
        for j in range(n):
            if i==j:
                matQ[i][j]=int(round(somE[i]/bet[i],0))
    return(matQ)
#------------------------------------------------------------------------------------------------------------
def chute(n):
    import random
    A=[]
    for i in range(n):
        A.append(random.randint(-3,3))
    return(A)
#------------------------------------------------------------------------------------------------------------
def pol(jac,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20):
    a=[]
    for i in range(len(jac)):
        if i==0:
            a.append(jac[i][0]*x1+jac[i][1]*x2+jac[i][2]*x3+jac[i][3]*x4+jac[i][4]*x5+jac[i][5]*x6+jac[i][6]*x7+jac[i][7]*x8+jac[i][8]*x9+jac[i][9]*x10+jac[i][10]*x11+jac[i][11]*x12+jac[i][12]*x13+jac[i][13]*x14+jac[i][14]*x15+jac[i][15]*x16+jac[i][16]*x17+jac[i][17]*x18+jac[i][18]*x19+jac[i][19]*x20-10)
        if i==1:
            a.append(jac[i][0]*x1+jac[i][1]*x2+jac[i][2]*x3+jac[i][3]*x4+jac[i][4]*x5+jac[i][5]*x6+jac[i][6]*x7+jac[i][7]*x8+jac[i][8]*x9+jac[i][9]*x10+jac[i][10]*x11+jac[i][11]*x12+jac[i][12]*x13+jac[i][13]*x14+jac[i][14]*x15+jac[i][15]*x16+jac[i][16]*x17+jac[i][17]*x18+jac[i][18]*x19+jac[i][19]*x20-20)
        if i==2:
            a.append(jac[i][0]*x1+jac[i][1]*x2+jac[i][2]*x3+jac[i][3]*x4+jac[i][4]*x5+jac[i][5]*x6+jac[i][6]*x7+jac[i][7]*x8+jac[i][8]*x9+jac[i][9]*x10+jac[i][10]*x11+jac[i][11]*x12+jac[i][12]*x13+jac[i][13]*x14+jac[i][14]*x15+jac[i][15]*x16+jac[i][16]*x17+jac[i][17]*x18+jac[i][18]*x19+jac[i][19]*x20-30)
        if i==3:
            a.append(jac[i][0]*x1+jac[i][1]*x2+jac[i][2]*x3+jac[i][3]*x4+jac[i][4]*x5+jac[i][5]*x6+jac[i][6]*x7+jac[i][7]*x8+jac[i][8]*x9+jac[i][9]*x10+jac[i][10]*x11+jac[i][11]*x12+jac[i][12]*x13+jac[i][13]*x14+jac[i][14]*x15+jac[i][15]*x16+jac[i][16]*x17+jac[i][17]*x18+jac[i][18]*x19+jac[i][19]*x20-40)
        if i==4:
            a.append(jac[i][0]*x1+jac[i][1]*x2+jac[i][2]*x3+jac[i][3]*x4+jac[i][4]*x5+jac[i][5]*x6+jac[i][6]*x7+jac[i][7]*x8+jac[i][8]*x9+jac[i][9]*x10+jac[i][10]*x11+jac[i][11]*x12+jac[i][12]*x13+jac[i][13]*x14+jac[i][14]*x15+jac[i][15]*x16+jac[i][16]*x17+jac[i][17]*x18+jac[i][18]*x19+jac[i][19]*x20-50)
        if i==5:
            a.append(jac[i][0]*x1+jac[i][1]*x2+jac[i][2]*x3+jac[i][3]*x4+jac[i][4]*x5+jac[i][5]*x6+jac[i][6]*x7+jac[i][7]*x8+jac[i][8]*x9+jac[i][9]*x10+jac[i][10]*x11+jac[i][11]*x12+jac[i][12]*x13+jac[i][13]*x14+jac[i][14]*x15+jac[i][15]*x16+jac[i][16]*x17+jac[i][17]*x18+jac[i][18]*x19+jac[i][19]*x20-60)
        if i==6:
            a.append(jac[i][0]*x1+jac[i][1]*x2+jac[i][2]*x3+jac[i][3]*x4+jac[i][4]*x5+jac[i][5]*x6+jac[i][6]*x7+jac[i][7]*x8+jac[i][8]*x9+jac[i][9]*x10+jac[i][10]*x11+jac[i][11]*x12+jac[i][12]*x13+jac[i][13]*x14+jac[i][14]*x15+jac[i][15]*x16+jac[i][16]*x17+jac[i][17]*x18+jac[i][18]*x19+jac[i][19]*x20-70)
        if i==7:
            a.append(jac[i][0]*x1+jac[i][1]*x2+jac[i][2]*x3+jac[i][3]*x4+jac[i][4]*x5+jac[i][5]*x6+jac[i][6]*x7+jac[i][7]*x8+jac[i][8]*x9+jac[i][9]*x10+jac[i][10]*x11+jac[i][11]*x12+jac[i][12]*x13+jac[i][13]*x14+jac[i][14]*x15+jac[i][15]*x16+jac[i][16]*x17+jac[i][17]*x18+jac[i][18]*x19+jac[i][19]*x20-80)
        if i==8:
            a.append(jac[i][0]*x1+jac[i][1]*x2+jac[i][2]*x3+jac[i][3]*x4+jac[i][4]*x5+jac[i][5]*x6+jac[i][6]*x7+jac[i][7]*x8+jac[i][8]*x9+jac[i][9]*x10+jac[i][10]*x11+jac[i][11]*x12+jac[i][12]*x13+jac[i][13]*x14+jac[i][14]*x15+jac[i][15]*x16+jac[i][16]*x17+jac[i][17]*x18+jac[i][18]*x19+jac[i][19]*x20-90)
        if i==9:
            a.append(jac[i][0]*x1+jac[i][1]*x2+jac[i][2]*x3+jac[i][3]*x4+jac[i][4]*x5+jac[i][5]*x6+jac[i][6]*x7+jac[i][7]*x8+jac[i][8]*x9+jac[i][9]*x10+jac[i][10]*x11+jac[i][11]*x12+jac[i][12]*x13+jac[i][13]*x14+jac[i][14]*x15+jac[i][15]*x16+jac[i][16]*x17+jac[i][17]*x18+jac[i][18]*x19+jac[i][19]*x20-100)
        if i==10:
            a.append(jac[i][0]*x1+jac[i][1]*x2+jac[i][2]*x3+jac[i][3]*x4+jac[i][4]*x5+jac[i][5]*x6+jac[i][6]*x7+jac[i][7]*x8+jac[i][8]*x9+jac[i][9]*x10+jac[i][10]*x11+jac[i][11]*x12+jac[i][12]*x13+jac[i][13]*x14+jac[i][14]*x15+jac[i][15]*x16+jac[i][16]*x17+jac[i][17]*x18+jac[i][18]*x19+jac[i][19]*x20-110)
        if i==11:
            a.append(jac[i][0]*x1+jac[i][1]*x2+jac[i][2]*x3+jac[i][3]*x4+jac[i][4]*x5+jac[i][5]*x6+jac[i][6]*x7+jac[i][7]*x8+jac[i][8]*x9+jac[i][9]*x10+jac[i][10]*x11+jac[i][11]*x12+jac[i][12]*x13+jac[i][13]*x14+jac[i][14]*x15+jac[i][15]*x16+jac[i][16]*x17+jac[i][17]*x18+jac[i][18]*x19+jac[i][19]*x20-120)
        if i==12:
            a.append(jac[i][0]*x1+jac[i][1]*x2+jac[i][2]*x3+jac[i][3]*x4+jac[i][4]*x5+jac[i][5]*x6+jac[i][6]*x7+jac[i][7]*x8+jac[i][8]*x9+jac[i][9]*x10+jac[i][10]*x11+jac[i][11]*x12+jac[i][12]*x13+jac[i][13]*x14+jac[i][14]*x15+jac[i][15]*x16+jac[i][16]*x17+jac[i][17]*x18+jac[i][18]*x19+jac[i][19]*x20-130)
        if i==13:
            a.append(jac[i][0]*x1+jac[i][1]*x2+jac[i][2]*x3+jac[i][3]*x4+jac[i][4]*x5+jac[i][5]*x6+jac[i][6]*x7+jac[i][7]*x8+jac[i][8]*x9+jac[i][9]*x10+jac[i][10]*x11+jac[i][11]*x12+jac[i][12]*x13+jac[i][13]*x14+jac[i][14]*x15+jac[i][15]*x16+jac[i][16]*x17+jac[i][17]*x18+jac[i][18]*x19+jac[i][19]*x20-140)
        if i==14:
            a.append(jac[i][0]*x1+jac[i][1]*x2+jac[i][2]*x3+jac[i][3]*x4+jac[i][4]*x5+jac[i][5]*x6+jac[i][6]*x7+jac[i][7]*x8+jac[i][8]*x9+jac[i][9]*x10+jac[i][10]*x11+jac[i][11]*x12+jac[i][12]*x13+jac[i][13]*x14+jac[i][14]*x15+jac[i][15]*x16+jac[i][16]*x17+jac[i][17]*x18+jac[i][18]*x19+jac[i][19]*x20-150)
        if i==15:
            a.append(jac[i][0]*x1+jac[i][1]*x2+jac[i][2]*x3+jac[i][3]*x4+jac[i][4]*x5+jac[i][5]*x6+jac[i][6]*x7+jac[i][7]*x8+jac[i][8]*x9+jac[i][9]*x10+jac[i][10]*x11+jac[i][11]*x12+jac[i][12]*x13+jac[i][13]*x14+jac[i][14]*x15+jac[i][15]*x16+jac[i][16]*x17+jac[i][17]*x18+jac[i][18]*x19+jac[i][19]*x20-160)
        if i==16:
            a.append(jac[i][0]*x1+jac[i][1]*x2+jac[i][2]*x3+jac[i][3]*x4+jac[i][4]*x5+jac[i][5]*x6+jac[i][6]*x7+jac[i][7]*x8+jac[i][8]*x9+jac[i][9]*x10+jac[i][10]*x11+jac[i][11]*x12+jac[i][12]*x13+jac[i][13]*x14+jac[i][14]*x15+jac[i][15]*x16+jac[i][16]*x17+jac[i][17]*x18+jac[i][18]*x19+jac[i][19]*x20-170)
        if i==17:
            a.append(jac[i][0]*x1+jac[i][1]*x2+jac[i][2]*x3+jac[i][3]*x4+jac[i][4]*x5+jac[i][5]*x6+jac[i][6]*x7+jac[i][7]*x8+jac[i][8]*x9+jac[i][9]*x10+jac[i][10]*x11+jac[i][11]*x12+jac[i][12]*x13+jac[i][13]*x14+jac[i][14]*x15+jac[i][15]*x16+jac[i][16]*x17+jac[i][17]*x18+jac[i][18]*x19+jac[i][19]*x20-180)
        if i==18:
            a.append(jac[i][0]*x1+jac[i][1]*x2+jac[i][2]*x3+jac[i][3]*x4+jac[i][4]*x5+jac[i][5]*x6+jac[i][6]*x7+jac[i][7]*x8+jac[i][8]*x9+jac[i][9]*x10+jac[i][10]*x11+jac[i][11]*x12+jac[i][12]*x13+jac[i][13]*x14+jac[i][14]*x15+jac[i][15]*x16+jac[i][16]*x17+jac[i][17]*x18+jac[i][18]*x19+jac[i][19]*x20-190)
        if i==19:
            a.append(jac[i][0]*x1+jac[i][1]*x2+jac[i][2]*x3+jac[i][3]*x4+jac[i][4]*x5+jac[i][5]*x6+jac[i][6]*x7+jac[i][7]*x8+jac[i][8]*x9+jac[i][9]*x10+jac[i][10]*x11+jac[i][11]*x12+jac[i][12]*x13+jac[i][13]*x14+jac[i][14]*x15+jac[i][15]*x16+jac[i][16]*x17+jac[i][17]*x18+jac[i][18]*x19+jac[i][19]*x20-200)
    return(a)
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
import numpy as np
import timeit
jac=jacobiana(-5,5)
tol=float(input('Tolerância: '))
n=int(input('Número de interações máxima: '))
x0=chute(20)
cont=0
x=x0
cond=False
inicio=timeit.default_timer()
while cont<n and cond==False:
    cont=cont+1
    auxjac=-np.linalg.inv(jac)
    poli=pol(jac,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15],x[16],x[17],x[18],x[19])
    increm=np.dot(auxjac,poli)
    x=x+increm
    if np.linalg.norm(increm)<tol:
        print('Os valores de x que satisfazem os coeficientes são: ')
        print(x)
        cond=True
print('Os novos valores de P(x) em x solução são: ')
k=pol(jac,x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15],x[16],x[17],x[18],x[19])
for i in range(len(k)):
    k[i]=round(k[i],10)
print(k)
fim = timeit.default_timer()
print('duracao algoritmo de Newton: %f' % (fim - inicio))
