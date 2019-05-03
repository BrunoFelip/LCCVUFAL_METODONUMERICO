# -*- coding: utf-8 -*-
"""
Created on Sat Apr 27 2019

@author: Bruno Felipe
"""    
#--------------------------------------------------------------------------------------------------------------------
def matconverg(n,Linf,Lsup):
    import random
    #valor de n>1 e quanto maior n (n>1000), maior o tempo de processamento (pode travar ou não o programa) e a precisão das raízes entre numpy e o método de Gauss podem divergir.
    #Linf e Lsup são os limitadores dos elementos aleatórios da matriz de coeficientes
    #Tomado por base a teoria do livro de cálculo numérico ASPECTOS TEÓRICOS E COMPUTACIONAIS 2° EDIÇÃO
    #será adotado um elemento matricial do subconjunto de matrizes do conj. IR^n->IR^n que satisfaçam o critério de Sassenfeld para eq. {x}=[C]{x}+{g}:
    #O algoritmo para gerar essa matriz deverá satisfazer as seguintes condições:
    #1)Será gerado uma lista delementos de ordem n:
    matQ=[]
    beta=[]
    bet=[]
    somE=[]
    for i in range(n):
        matQ.append([0]*n)
        beta.append([0]*n)
    for i in range(n):
        somE.append(0)
    #2)Todos os elementos, exceto a diag. principal, recebem qualquer elemento da função randint.
    for i in range(n):
        for j in range(n):
            if i!=j:
                matQ[i][j]=random.randint(Linf,Lsup)
    #3)É montado a matriz beta por meio do vetor bet, onde cada elemento será múltiplicado pelos coef. da matriz matQ, obtendo cada elemento das somas da Linha matQ.
    #Exemplo
    #n=4
    # L0: somE0=a00*0+a01*1+a02*1+a03*1
    # L1: somE1=a10*bet0+a11*0+a12*1+a13*1
    # L2: somE2=a20*bet0+a21*bet1+a22*0+a23*1
    # L3: somE3=a30*bet0+a31*bet1+a32*bet2+a33*0
    #bet=[bet0,bet1,bet2,bet3->adotado]
    #matQ=[[a00,a01,a02,a03],[a10,a11,a12,a13],[a20,a21,a22,a23],[a30,a31,a32,a33]]
    #beta=[[0,1,1,1],[bet0,0,1,1],[bet0,bet1,0,1],[bet0,bet1,bet2,0]]
    for i in range(n):
        #adotado valores aleatórios para beta entre 0.1 e 0.9, levando em consideração o arredondamento dos elementos da diagonal principal.
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
    #Cada elemento da diagonal principal é determinado pela eq. beti=somEi/aii.
    for i in range(n):
        for j in range(n):
            somE[i]=somE[i]+abs(matQ[i][j]*beta[i][j])
    for i in range(n):
        for j in range(n):
            if i==j:
                matQ[i][j]=int(round(somE[i]/bet[i],0))
    #Garantido que qualquer chute inicial {X} convergirá para uma única solução, já que det(matQ)>0.
    return(matQ)
#----------------------------------------------------------------------------------------------------------------------------
#Será considerado que a matriz gerada é uma matriz qualquer, onde não se sabe se o critério de senssefeld é ou não satisfeito.
def sassenfeld(matQ):
    n=len(matQ)
    #vetor de elementos neutros aditivos
    bet=[]
    #vetor de elementos neutros multiplicativos e vetor resposta
    beta=[]
    for i in range(n):
        bet.append(0)
        beta.append(1)
    for i in range(n):
        div=matQ[i][i]
        for j in range(n):
            if i!=j:
                bet[i]=bet[i]+abs(matQ[i][j]*beta[j]/div)
        beta[i]=bet[i]
    return(beta)
#---------------------------------------------------------------------------------------------------------------------------
def resolv(matQ,vetB,lim,interac):
#seja [A]{X}=[B], onde [A]=matQ,[B]=vetB e {x}=vetX, o sistema a ser resolvido, onde a matriz satisfaz o critério de sessenfeld, logo será adotado Xij=0, qualquer i,j natural.
    X=[]
    x=[]
    r=[]
    s=[]
    t=[]
    for i in range(len(matQ)):
        #Elementos neutros da adição
        x.append(0)
        r.append(0)
        s.append(0)
        t.append(0)
        X.append(0)
    cond=False
    cont=0
    contlim=0
    while cond==False and contlim<interac:
            if cont>0:
                for i in range(len(r)):
                    r[i]=X[i]
            for i in range(len(matQ)):
                x[i]=vetB[i]/matQ[i][i]
                for j in range(len(matQ)):
                    if i!=j:
                        x[i]=x[i]-matQ[i][j]*X[j]/matQ[i][i]
                X[i]=x[i]
            for i in range(len(r)):
                t[i]=abs(X[i])
                s[i]=abs(X[i]-r[i])
            a=max(s,key=float)
            b=max(t,key=float)
            incr=a/b
            if incr<lim:
                cond=True
            cont=1
            contlim=contlim+1
    print('Número de tentativas: {}'.format(contlim))
    return(X)
#----------------------------------------------------------------------------------------------------------------------------
import timeit
import random
import numpy as np
n=input('Digite a ordem da matriz n de coeficientes: ')
if n.isnumeric():
    n=int(n)
    if n>1:
        matQ=matconverg(n,-50,50)
        b=[]
        for i in range(n):
            b.append(random.randint(-500,500))
        a=sassenfeld(matQ)
        a=max(a)
        if a<1:
            print('Matriz de Coeficientes Convergentes.')
            inicio=timeit.default_timer()
            c=resolv(matQ,b,0.0000001,100000)
            fim = timeit.default_timer()
            print('duracao gauss-seidl: %f' % (fim - inicio))
            inicio1=timeit.default_timer()
            d=np.linalg.solve(matQ,b)
            fim1 = timeit.default_timer()
            print('duracao np.linalg.solve: %f' % (fim1 - inicio1))
            #Verificar o quando se o Gauss-Seidl apresenta uma resposta coerente:
            e=np.round(np.dot(matQ,c)-b, decimals=7)
        else:
            print('Matriz não convergente.')
    else:
        print('Digite um valor de n válido: ')
else:
    print('Opção inválida.')
print('Resultado pelo algoritmo desenvolvido: ')
print(c)
print('verificando se os resultados obtidos pelo meu algoritmo são precisos(quando todos forem nulos): ')
print(e)
print('Resultado pelo algoritmo linalg.solve: ')
print(d)