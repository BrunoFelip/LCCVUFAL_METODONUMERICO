# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 08:52:51 2019

@author: Bruno Felipe
"""

#integração numérica Regra do Trapézio Repetido
import numpy as np
import timeit
from random import randint
pol=np.array([0]*7,dtype=float)
for i in range(7):
    pol[i]=randint(-15,15)
p=np.poly1d(pol)
diff=np.poly1d.deriv(p,m=2)
h=-1
while h<=0:
    xi=int(input('xi: '))
    xf=int(input('xf: '))
    n=int(input('Espaçamento: '))
    h=(xf-xi)/n
    print('Valor de h: {}'.format(h))
inicio=timeit.default_timer()
if n==0:
    print('O espaçamento n é 0.')
else:
    if n<0:
        print('Espaçamento n é menor que 0.')
    else:
        x=xi
        I=0
        der=[]
        s=(xf-xi)/0.1
        for i in range(int(s+1)):
            aux=xi+i*s
            der.append(diff(x))
        for i in range(n+1):
            x=xi+i*h
            if i==0:
                I=I+p(x)
            elif i==n:
                I=I+p(x)
            else:
                I=I+2*p(x)
        I=I*h/2
        Err=((xf-xi)**3/(12*(n**2)))*max(der)
        IErrsup=I+Err
        IErrinf=I-Err
fim = timeit.default_timer()
print('Função: ')
print(np.poly1d(p))
print('A integral da função: {}. O erro estimado: {}. Intervalo de solução: {}<{}<{}'.format(I,Err,IErrinf,I,IErrsup))
print('duracao algoritmo trapézio: %f' % (fim - inicio))
print('')
print('Módulo Scipy Integração')
inicio1=timeit.default_timer()
from scipy.integrate import quad
Integr = quad(p, xi, xf)
print(Integr)
fim1 = timeit.default_timer()
print('duracao algoritmo scipy: %f' % (fim1 - inicio1))