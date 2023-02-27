#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 23:48:50 2019

@author: geoffreylau
"""

# initialization

import numpy as np
from scipy import optimize

# part 1

# input of all 3 variables, declaring how many slack variables are introduced

A = np.array([[1/20,1/16,1/11,1/5,1,0],[1/10,1/15,1/24,1/10,0,1]],dtype = 'f')
B = np.array([8,8],dtype = 'f')
C = np.array([-12,-9,-8,-10,0,0],dtype = 'f')

# amount of slack variables, as we can assume the first basis as slack variables in this assignment

slack = 2


# finding size of A

numrowsA = len(A)
numcolsA = len(A[0])




# part 2

# setting up beta and mu, choosing the last columns as basis

mu = np.arange(0,numcolsA - slack)
beta = np.arange(numcolsA - slack, numcolsA)

# A beta

A_beta = np.zeros((numrowsA,slack),dtype = 'f')
for y in range(slack):
    for x in range(numrowsA):
        A_beta[x][y] = A[x][beta[y]]
        x = x + 1
    y = y + 1
    
# A mu
A_mu = np.zeros((numrowsA,numcolsA - slack),dtype = 'f')
for y in range(numcolsA-slack):
    for x in range(numrowsA):
        A_mu[x][y] = A[x][mu[y]]
        x = x + 1
    y = y + 1
    
# C beta

C_beta = np.zeros(slack,dtype = 'f')
for x in range(slack):
    C_beta[x]=C[beta[x]]
    x = x + 1
    
# C mu
    
C_mu = np.zeros(numcolsA - slack,dtype = 'f')
for x in range(numcolsA - slack):
    C_mu[x]=C[mu[x]]
    x = x + 1

# begin iteration
    
optimal = False
feasible = True

while optimal == False:

    
    
# part 3    
    
    
# calculate b bar
    
    A_beta_inv = np.linalg.inv(A_beta)
    b_bar = np.matmul(A_beta_inv,B)
    count = 0
    for x in range(len(b_bar)):
        if(b_bar[x]<0):
            count = count + 1
    x_current = np.array([0]*numcolsA,dtype = 'f')
    for x in range(len(beta)):
        x_current[beta[x]] = b_bar[x]
    print(x_current, "is being tested.")
    
# checking if it is a BFS
    
    if(count == 0):
        print ("We have basic feasible solution!")
        print ("Now we check for optimality.")
        feasible = True
    else:
        print ("Not feasible!")
        
# If not BFS we stop immediately by triggering boolean
        
        optimal = True
        feasible = False

    if(optimal == False and feasible == True):
        
# part 4                

# Finding simplex multipliers
        
        y = np.matmul(np.linalg.inv(np.transpose(A_beta)),C_beta)
        
# Computing reduced cost
        
        r = np.subtract(C_mu,np.matmul(np.transpose(A_mu),y))
        


        count2 = 0
        for x in range(len(r)):
            if(r[x]>=0):
                count2 = count2 + 1
        if (count2 == len(r)):
            
# If optimal we can leave loop by triggering boolean
            
            print("Optimal!")
            optimal = True
        else:
            print("Not optimal!")
            
# Finding the smallest value in the reduced cost, locating its position in the matrix
        
        entering_index = np.argmin(r)
        
    if(optimal == False and feasible == True):

# part 5

        
# Finding the smallest value in the basis, locating its position in the matrix
# Using ratio to calculate the maximum t values, and choosing the smallest ratio to leave the basis
        
        a_temp = np.array([0]*numrowsA,dtype = 'f')
        x_beta_ratio = np.array([0]*numrowsA,dtype = 'f')
        for x in range(numrowsA):
            a_temp[x] = A_mu[x][entering_index]
        a_bar = np.matmul(np.linalg.inv(A_beta),a_temp)
        for x in range(numrowsA):
            x_beta_ratio[x] = b_bar[x] / a_bar[x]
            if (x_beta_ratio[x]<= 0):
                
# since equality, we dont want any ratio that is negative or 0, so we just randomly append... 
# ... a huge value so it wont be the smallest ratio that we are choosing
                
                x_beta_ratio[x] = 10000000
        leaving_index = np.argmin(x_beta_ratio)
        # Replacing mu and beta with the two locations we calculated
        C_beta[leaving_index], C_mu[entering_index] = C_mu[entering_index], C_beta[leaving_index]
        beta[leaving_index], mu[entering_index] = mu[entering_index], beta[leaving_index]
        for x in range(numrowsA):
            A_beta[x][leaving_index], A_mu[x][entering_index] = A_mu[x][entering_index],A_beta[x][leaving_index]

# end of iteration



# printing our optimal solution 
            
x_opt = np.array([0]*numcolsA,dtype = 'f')
for x in range(len(beta)):
    x_opt[beta[x]] = b_bar[x]
print(x_opt, "is the optimal solution.")

# checking our program with the linear programming of scipystats

print("Now we use the scipy linprog function to check if our coding is correct.")
check = optimize.linprog(C,A_eq=A, b_eq=B, method = 'simplex')
print(check)