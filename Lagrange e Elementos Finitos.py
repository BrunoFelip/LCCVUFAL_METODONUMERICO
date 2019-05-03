# -*- coding: utf-8 -*-
"""
Created on Wed May  1 17:22:04 2019

@author: Bruno Felipe
"""

#Tomado por base na apostila do prof. William UFAL (Dezembro 2009) Pág 59 - Algoritmo Base
def lagr(x,fx,tipe='x'):
    n=len(x)
    import sympy as sy
    p=0
    var=sy.symbols(tipe)
    for i in range(n):
        s=1
        for j in range(n):
            if j!=i:
                s=s*(var-x[j])/(x[i]-x[j])
        p=p+fx[i]*s
    p=sy.factor(p)
    return(p)
#--------------------------------------------------------
#Atenção ao conferir adequadamente o polinômio simbólico gerado com
#os dados das equações de forma. Exemplo: (x-1)(y-1)=-(1-x)(-(1-y))=(1-x)(1-y)
import timeit
import sympy as sy
inicio=timeit.default_timer()
print('Ponto 1 N1')

#Parábola Trecho (1,5,2)
N1x=lagr([-1,0,1],[1,0,0],'r')
print(N1x)
#Parábola Trecho (1,8,4)
N1y=lagr([-1,0,1],[1,0,0],'s')
print(N1y)
N1xy=N1x*N1y
print(N1xy)
print('Expandido: ', end='')
EN1xy=sy.expand(N1xy)
print(EN1xy)
print('')

print('Ponto 2 N2')

#Parábola Trecho (1,5,2)
N2x=lagr([-1,0,1],[0,0,1],'r')
print(N2x)
#Parábola Trecho (2,6,3)
N2y=lagr([-1,0,1],[1,0,0],'s')
print(N2y)
N2xy=N2x*N2y
print(N2xy)
print('Expandido: ', end='')
EN2xy=sy.expand(N2xy)
print(EN2xy)
print('')

print('Ponto 3 N3')

#Parábola Trecho [2,6,3]
N3x=lagr([-1,0,1],[0,0,1],'r')
print(N3x)
#Parábola Trecho [4,7,3]
N3y=lagr([-1,0,1],[0,0,1],'s')
print(N3y)
N3xy=N3x*N3y
print(N3xy)
print('Expandido: ', end='')
EN3xy=sy.expand(N3xy)
print(EN3xy)
print('')

print('Ponto 4 N4')

#Parábola Trecho [4,7,5]
N4x=lagr([-1,0,1],[1,0,0],'r')
print(N4x)
#Parábola Trecho [1,8,4]
N4y=lagr([-1,0,1],[0,0,1],'s')
print(N4y)
N4xy=N4x*N4y
print(N4xy)
print('Expandido: ', end='')
EN4xy=sy.expand(N4xy)
print(EN4xy)
print('')

print('Ponto 5 N5')
    
#Parábola de segundo Grau Trecho (1,5,2)=>x(-1,0,1)->fx(0,1,0)
N5x=lagr([-1,0,1],[0,1,0],'r')
print(N5x)
#Parábola(5,9,7)=>x(-1,1,1)=>fx(1,0,0) 
N5y=lagr([-1,0,1],[1,0,0],'s')
print(N5y)
N5xy=N5x*N5y
print(N5xy)
print('Expandido: ', end='')
EN5xy=sy.expand(N5xy)
print(EN5xy)
print('')

print('Ponto 6 N6')

#Parábola Trecho [8,9,6]
N6x=lagr([-1,0,1],[0,0,1],'r')
print(N6x)
#Parábola Trecho [2,6,3]
N6y=lagr([-1,0,1],[0,1,0],'s')
print(N6y)
N6xy=N6x*N6y
print(N6xy)
print('Expandido: ', end='')
EN6xy=sy.expand(N6xy)
print(EN6xy)
print('')

print('Ponto 7 N7')

#Parábola Trecho [4,7,3]
N7x=lagr([-1,0,1],[0,1,0],'r')
print(N7x)
#Parábola Trecho [5,9,7]
N7y=lagr([-1,0,1],[0,0,1],'s')
print(N7y)
N7xy=N7x*N7y
print(N7xy)
print('Expandido: ', end='')
EN7xy=sy.expand(N7xy)
print(EN7xy)
print('')

print('Ponto 8 N8')

#Parábola Trecho [6,9,8]
N8x=lagr([-1,0,1],[-2,0,0],'r')
print(N8x)
#Parábola Trecho [1,8,4]
N8y=lagr([-1,0,1],[1,0,0],'s')
print(N8y)
N8xy=N8x*N8y
print(N8xy)
print('Expandido: ', end='')
EN8xy=sy.expand(N8xy)
print(EN8xy)
print('')

print('Ponto 9 N9')

#Parábola Trecho [5,9,7]
N9x=lagr([-1,0,1],[0,1,0],'r')
print(N9x)
#Parábola Trecho [8,9,6]
N9y=lagr([-1,0,1],[0,1,0],'s')
print(N9y)
N9xy=N9x*N9y
print(N9xy)
print('Expandido: ', end='')
EN9xy=sy.expand(N9xy)
print(EN9xy)
print('')

fim = timeit.default_timer()
print('duracao algoritmo de lagrange programada: %f' % (fim - inicio))
#-------------------------------------------------------------------------------------------------
from scipy.interpolate import lagrange

inicio1=timeit.default_timer()
print('Com scipy e sympy')
print('')
print('Ponto 1 N1')

r=sy.Symbol('r')
s=sy.Symbol('s')

cN1x=lagrange([-1,0,1],[1,0,0])
aN1x=cN1x.c
vN1x=cN1x[2]*r**2+cN1x[1]*r**1+cN1x[0]

cN1y=lagrange([-1,0,1],[1,0,0])
aN1y=cN1y.c
vN1y=cN1y[2]*s**2+cN1y[1]*s**1+cN1y[0]
vN1xy=vN1x*vN1y
print('Expandido: ', end='')
VN1xy=sy.expand(vN1xy)
print(VN1xy)
print('')

print('Ponto 2 N2')


cN2x=lagrange([-1,0,1],[0,0,1])
aN2x=cN2x.c
vN2x=cN2x[2]*r**2+cN2x[1]*r**1+cN2x[0]

cN2y=lagrange([-1,0,1],[1,0,0])
aN2y=cN2y.c
vN2y=cN2y[2]*s**2+cN2y[1]*s**1+cN2y[0]
vN2xy=vN2x*vN2y
print('Expandido: ', end='')
VN2xy=sy.expand(vN2xy)
print(VN2xy)
print('')

print('Ponto 3 N3')


cN3x=lagrange([-1,0,1],[0,0,1])
aN3x=cN3x.c
vN3x=cN3x[2]*r**2+cN3x[1]*r**1+cN3x[0]

cN3y=lagrange([-1,0,1],[0,0,1])
aN3y=cN3y.c
vN3y=cN3y[2]*s**2+cN3y[1]*s**1+cN3y[0]
vN3xy=vN3x*vN3y
print('Expandido: ', end='')
VN3xy=sy.expand(vN3xy)
print(VN3xy)
print('')

print('Ponto 4 N4')


cN4x=lagrange([-1,0,1],[1,0,0])
aN4x=cN4x.c
vN4x=cN4x[2]*r**2+cN4x[1]*r**1+cN4x[0]

cN4y=lagrange([-1,0,1],[0,0,1])
aN4y=cN4y.c
vN4y=cN4y[2]*s**2+cN4y[1]*s**1+cN4y[0]
vN4xy=vN4x*vN4y
print('Expandido: ', end='')
VN4xy=sy.expand(vN4xy)
print(VN4xy)
print('')

print('Ponto 5 N5')


cN5x=lagrange([-1,0,1],[0,1,0])
aN5x=cN5x.c
vN5x=cN5x[2]*r**2+cN5x[1]*r**1+cN5x[0]

cN5y=lagrange([-1,0,1],[1,0,0])
aN5y=cN5y.c
vN5y=cN5y[2]*s**2+cN5y[1]*s**1+cN5y[0]
vN5xy=vN5x*vN5y
print('Expandido: ', end='')
VN5xy=sy.expand(vN5xy)
print(VN5xy)
print('')

print('Ponto 6 N6')


cN6x=lagrange([-1,0,1],[0,0,1])
aN6x=cN5x.c
vN6x=cN6x[2]*r**2+cN6x[1]*r**1+cN6x[0]

cN6y=lagrange([-1,0,1],[0,1,0])
aN6y=cN6y.c
vN6y=cN6y[2]*s**2+cN6y[1]*s**1+cN6y[0]
vN6xy=vN6x*vN6y
print('Expandido: ', end='')
VN6xy=sy.expand(vN6xy)
print(VN6xy)
print('')

print('Ponto 7 N7')


cN7x=lagrange([-1,0,1],[0,1,0])
aN7x=cN7x.c
vN7x=cN7x[2]*r**2+cN7x[1]*r**1+cN7x[0]

cN7y=lagrange([-1,0,1],[0,0,1])
aN7y=cN7y.c
vN7y=cN7y[2]*s**2+cN7y[1]*s**1+cN7y[0]
vN7xy=vN7x*vN7y
print('Expandido: ', end='')
VN7xy=sy.expand(vN7xy)
print(VN7xy)
print('')

print('Ponto 8 N8')


cN8x=lagrange([-1,0,1],[-2,0,0])
aN8x=cN8x.c
vN8x=cN8x[2]*r**2+cN8x[1]*r**1+cN8x[0]

cN8y=lagrange([-1,0,1],[1,0,0])
aN8y=cN8y.c
vN8y=cN8y[2]*s**2+cN8y[1]*s**1+cN8y[0]
vN8xy=vN8x*vN8y
print('Expandido: ', end='')
VN8xy=sy.expand(vN8xy)
print(VN8xy)
print('')

print('Ponto 9 N9')


cN9x=lagrange([-1,0,1],[0,1,0])
aN9x=cN9x.c
vN9x=cN9x[2]*r**2+cN9x[1]*r**1+cN9x[0]

cN9y=lagrange([-1,0,1],[0,1,0])
aN9y=cN9y.c
vN9y=cN9y[2]*s**2+cN9y[1]*s**1+cN9y[0]
vN9xy=vN9x*vN9y
print('Expandido: ', end='')
VN9xy=sy.expand(vN9xy)
print(VN9xy)
print('')

fim1 = timeit.default_timer()
print('duracao algoritmo scipy.interpolate: %f' % (fim1 - inicio1))