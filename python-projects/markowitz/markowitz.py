#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 23:14:49 2019

@author: geoffreylau
"""
## initialization
import numpy as np
import cvxopt as cv
import matplotlib.pyplot as plt

## the following code is equivalent to the matlab code provide in the file 
## hw_2_data.m

n = 8
corr = np.zeros((n,n))
for i in range (n):
    for j in range (n):
        corr[i][j] = pow(-1,abs(i-j))/(abs(i-j)+1)
sigma = np.zeros((n,1))
mu = np.zeros((n,1))
sigma[0] = 2
mu[0] = 3
for i in range(0,n-1):
    sigma[i+1] = sigma[i] + 2 * np.random.random_sample()
    mu[i+1] = mu[i] + 1
D = np.diagflat(sigma)
C2 = D @ corr @ D
C = 0.5 * (C2 + np.transpose(C2))

## end of translation of given matlab code

##stopping cvxopt showing optimization progress
cv.solvers.options['show_progress'] = False 

## Exercise 1 begins
## setting up the 25 values for r, from 3 to 9 with steps of 0.25
r = np.linspace(3.,9.,25)
## formatting matrix into cvxopt matrix, preparing for qp function
u = cv.matrix(mu)
p = cv.matrix(C)
q = cv.matrix([0.0 for i in range(n)])
G = - cv.matrix(np.eye(n))
h = cv.matrix(0.0, (n ,1)) 
A = cv.matrix([[1.0 for i in range(n)],[u]]).trans()
## sig and mean stores the standard deviation and mean of the 25 values
sig = np.zeros(len(r))
mean = np.zeros(len(r))
## running qp for 25 times, storing into sig and mean array each time
for i in range(len(r)):
    b = cv.matrix([1.0,r[i]])
    temp = cv.solvers.qp(p,q,G,h,A,b)
    temp2 = np.matrix(temp['x'])
    sig[i] = np.sqrt(temp2.transpose() @ C @ temp2)
    mean[i] = mu.transpose() @ temp2

plt.plot(sig,mean,'o-')
plt.ylabel('mean')  
plt.xlabel('std')
plt.title('Exercise 1')
plt.show()
## Excerise 1 ends

## Exercise 2 begins

## reformatting matrix with the changes in constraints

e = cv.matrix(1.0,(1,n))
G2 = cv.matrix([G,e],(n+1,n))
h2 = cv.matrix([h,1.0],(n+1,1))
A2 = cv.matrix([u]).trans()
## sig2 and mean2 stores the standard deviation and mean of the 25 values
sig2 = np.zeros(len(r))
mean2 = np.zeros(len(r))
for i in range(len(r)):
    b2 = cv.matrix([r[i]])
    temp = cv.solvers.qp(p,q,G2,h2,A2,b2)
    temp2 = np.matrix(temp['x'])
    sig2[i] = np.sqrt(temp2.transpose() @ C @ temp2)
    mean2[i] = mu.transpose() @ temp2
## plotting the original exercise with this exercise
plt.plot(sig,mean,'o-')
plt.plot(sig2,mean2,'yo-')
plt.ylabel('mean')  
plt.xlabel('std')
plt.title('Exercise 2')
plt.show()
## Exercise 2 ends

## Exercise 3 begins

## reformatting due to constraint changes
G3 = cv.matrix([G,-u.trans()],(n+1,n))
A3 = e
b3 = cv.matrix(1.0)
sig3 = np.zeros(len(r))
mean3 = np.zeros(len(r))
## sig3 and mean3 stores the standard deviation and mean of the 25 values
for i in range(len(r)):
    h3 = cv.matrix([h,-r[i]],(n+1,1))
    temp = cv.solvers.qp(p,q,G3,h3,A3,b3)
    temp2 = np.matrix(temp['x'])
    sig3[i] = np.sqrt(temp2.transpose() @ C @ temp2)
    mean3[i] = mu.transpose() @ temp2
## plotting the original exercise with this exercise
plt.plot(sig,mean,'o-')
plt.plot(sig3,mean3,'go-')
plt.ylabel('mean')  
plt.xlabel('std')
plt.title('Exercise 3')
plt.show()
## Exercise 3 ends

## Exercise 4 begins
## reformatting due to constraint changes
G4 = cv.matrix(0.0,(n,n))
sig4 = np.zeros(len(r))
mean4 = np.zeros(len(r))
## sig4 and mean4 stores the standard deviation and mean of the 25 values
for i in range(len(r)):
    b = cv.matrix([1.0,r[i]])
    temp = cv.solvers.qp(p,q,G4,h,A,b)
    temp2 = np.matrix(temp['x'])
    sig4[i] = np.sqrt(temp2.transpose() @ C @ temp2)
    mean4[i] = mu.transpose() @ temp2
## plotting the original exercise with this exercise
plt.plot(sig,mean,'o-')
plt.plot(sig4,mean4,'ro-')
plt.ylabel('mean')  
plt.xlabel('std')
plt.title('Exercise 4')
plt.show()
## Exercise 4 ends