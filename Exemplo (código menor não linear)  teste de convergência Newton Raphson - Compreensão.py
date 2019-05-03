# -*- coding: utf-8 -*-
"""
Created on Thu May  2 05:10:40 2019

@author: Bruno Felipe
"""
#Exemplo do site de cálculo numérico UFRGS para uma dimensão pequena RESPOSTA [1,34681;0,46]
def F(x):
    import numpy as np
    y = np.zeros(2)  
 
    y[0] = x[0]**2 - np.cos(x[0]*x[1]) - 1  
    y[1] = np.sin(x[1]) - 2*np.cos(x[0])  
 
    return y  
 
def JF(x):
    import numpy as np
    y = np.zeros((2,2))  
 
    y[0,0] = 2*x[0] + x[1]*np.sin(x[0]*x[1])  
    y[0,1] = x[0]*np.sin(x[0]*x[1])  
 
    y[1,0] =  2*np.sin(x[0])  
    y[1,1] = np.cos(x[1])
    return y

import numpy as np
tol=float(input('Tolerância: '))
n=int(input('Número de interações máxima: '))
x0=[1.5,0.5]
cont=0
x=x0
cond=False
while cont<n and cond==False:
    cont=cont+1
    auxjac=-np.linalg.inv(JF(x))
    poli=F(x)
    increm=np.dot(auxjac,poli)
    x=x+increm
    if np.linalg.norm(increm)<tol:
        print(x)
        cond=True