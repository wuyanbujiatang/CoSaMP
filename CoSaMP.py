#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 30 15:49:55 2018
@author: wsj

"""
import numpy as np
import random
import matplotlib.pyplot as plt
import copy
def CoSaMP(y,phi,psi,K):
    A = phi*psi.H
    M,N = np.shape(A)
    v = y
    s = np.zeros((N,1))
    pos_reserved = []
    for i in range(K):
        product = A.H*v
        pos = np.argsort(-abs(product.T))[0,0:2*K]
        pos_total = pos_reserved
        for j in range(2*K):
            if pos[0,j] not in pos_reserved:
                pos_total.append(pos[0,j])
        if len(pos_total)<=M:
            At = copy.copy(A[:,pos_total])
        else:
            if i == 1: 
                ls_reserved = 0
            break
        ls = (At.H*At).I*At.H*y
#        print(ls)
        pos_theta = np.argsort(-(abs(ls.T)))
        pos_k = []
        for kk in range(K):
            pos_k.append(pos_theta[0,kk])
        pos_reserved = []
        for kk in range(len(pos_k)):
            pos_reserved.append(pos_total[pos_k[kk]])
        print(pos_reserved)
        ls_reserved = ls[pos_k]
        v = y - At[:,pos_k]*ls_reserved
        if np.linalg.norm(v)<1e-12:
            break
    s[pos_reserved] = ls_reserved
    return psi.H*s
def my_plot(x,x_r):
    fig,axes = plt.subplots(nrows=1,ncols=1,figsize=(8,6))  
    line1, = axes.plot(x_r,'r-s')
    line2, = axes.plot(x,'b--o')
    axes.legend((line1,line2),('recovered signal','src signal'),loc = 'upper right')
    axes.set_xlabel(u'time')
    axes.set_ylabel(u'signal')
    axes.set_title('CoSaMP algorithm')
if __name__ == '__main__':
    N = 256
    M = 64
    K = 8
    x = np.zeros((N,1))
    index_k = random.sample(range(N),K)
    x[index_k] = 5*np.random.randn(K,1)
    psi = np.eye(N,N)
    phi = np.random.randn(M,N)
    phi = np.mat(phi)
    psi = np.mat(psi)
    x = np.mat(x)
    y = phi*psi.H*x
    my_plot(CoSaMP(y,phi,psi,K),x)
#    print(mp(y,phi,psi,K)-x)
    

