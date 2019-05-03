# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 02:36:36 2019

@author: Bruno Felipe
"""
#Ajuste Polinomial
#n referencia o grau do polinômio desejado
def pol(n,x,y):
    import numpy as np
    C=np.zeros((n,n),dtype=np.float64)
    b=np.array([0]*n,float)
    for i in range(n):
        for j in range(n):
            C[i][j]=0
            for k in range(len(x)):
                C[i][j]=C[i][j]+x[k]**(i+j)
            C[j][i]=C[i][j]
        b[i]=0
        for k in range(len(x)):
            b[i]=b[i]+y[k]*x[k]**(i)
    D=np.linalg.inv(C)
    r=np.dot(D,b)
    s=np.array([0]*len(r),dtype=float)
    for i in range(len(r)):
        s[i]=round(r[len(r)-1-i],10)
    p=np.poly1d(s)
    return(p)
#-------------------------------------------------------------------------
def qualid(p,x,fx):
    medm=0
    for i in range(len(x)):
        medm=medm+fx[i]/len(x)
    num=0
    dem=0
    for i in range(len(x)):
            num=num+(p(x[i])-medm)**2
            dem=dem+(fx[i]-medm)**2
    if dem!=0:
        rquad=(num/dem)
        return(rquad)
    else:
        return('Divisão por zero.')
#-------------------------------------------------------------------------
import random
import numpy as np
x=[]
y=[]
n=random.randint(1,60)
for i in range(n):
    x.append(random.randint(1,500))
    y.append(random.randint(1,500))
import timeit
inicio=timeit.default_timer()
print('x= ',end='')
print(x)
print('y= ',end='')
print(y)
print('Algoritmo elaborado: ')
for i in range(1,6,1):
    p=pol(i,x,y)
    r=qualid(p,x,y)
    print(np.poly1d(p))
    print('r^2={}'.format(r))
    if i==1:
        p0=p
    elif i==2:
        p1=p
    elif i==3:
        p2=p
    elif i==4:
        p3=p
    elif i==5:
        p4=p
fim = timeit.default_timer()
print('duracão do algoritmo elaborado: %f' % (fim - inicio))
print('#-----------------------------------------------------')
print('#-----------------------------------------------------')
#----------------------------------------------------------------------
inicio1=timeit.default_timer()
a0=np.polyfit(x,y,0)
print(np.poly1d(a0))
a1=np.polyfit(x,y,1)
print(np.poly1d(a1))
a2=np.polyfit(x,y,2)
print(np.poly1d(a2))
a3=np.polyfit(x,y,3)
print(np.poly1d(a3))
a4=np.polyfit(x,y,4)
print(np.poly1d(a4))
fim1 = timeit.default_timer()
print('duracão da interpolação com numpy: %f' % (fim1 - inicio1))
print('#-----------------------------------------------------')
print('#-----------------------------------------------------')
#---------------------------------------------------------------------------------
from scipy.optimize import curve_fit
inicio2=timeit.default_timer()
def func4(x,a,b,c,d,e):
    return a*x**4+b*x**3+c*x**2+d*x+e
popt4,pcov4=curve_fit(func4, x, y)
b4=np.poly1d(popt4)
def func3(x,a,b,c,d):
    return a*x**3+b*x**2+c*x+d
popt3,pcov3=curve_fit(func3, x, y)
b3=np.poly1d(popt3)
def func2(x,a,b,c):
    return a*x**2+b*x+c
popt2,pcov2=curve_fit(func2, x, y)
b2=np.poly1d(popt2)
def func1(x,a,b):
    return a*x+b
popt1,pcov1=curve_fit(func1, x, y)
b1=np.poly1d(popt1)
def func0(x,a):
    return a*x**0
popt0,pcov0=curve_fit(func0, x, y)
b0=np.poly1d(popt0)
print(np.poly1d(b0))
print(np.poly1d(b1))
print(np.poly1d(b2))
print(np.poly1d(b3))
print(np.poly1d(b4))
fim2 = timeit.default_timer()
print('duracão da interpolação com scipy: %f' % (fim2 - inicio2))



