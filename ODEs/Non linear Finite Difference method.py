# -*- coding: utf-8 -*-
"""
Created on Sun Nov 15 21:39:07 2020

@author: Ahmad Aiman B. Mohd Nazir
"""
#import the libraries 
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve

#define the exact solution , u(x)
def uexact(x):
    return np.log((np.e - 1)*x + 1)
#define the f = u"
def f(x,u,uprime):
    return -uprime**2

# set the boundary conditions where u(x0) = alpha and u(xn) = beta 
x0 = 0
xn = 1
alpha = 0
beta = 1
#set the number of iterations
n = 19
h = (xn - x0)/(n+1)     #calculate the step size

#create the array of list for xi values
x = []
for i in range(0,n+2):
    x.append(x0 + i*h)
    i+=1
#calculate the exact solution and create the list of array for u(xi)
yexact = []
for i in x:
    yexact.append(uexact(i))
#calculate the initial approximation
m = (beta - alpha)/(xn - x0)    #calculate the gradient
#create array of p and calculate the values of u(xi) using Newton-Raphson method y = mx + c
p = []
for i in range(0,n+2):
    p.append(alpha + i*h*m)
    i+=1
    
p1 = [p]    #create another array with same values in array p

#create the matrix A with size nxn and input the values using for loop
A = np.zeros((n,n))
for i in range(0,n):
    for j in range(0,n):
        if i == j:
            A[i][j] = -2
        elif i<j and j-i==1:
            A[i][j] = 1
        elif i>j and i-j==1:
            A[i][j] = 1
        else:
            A[i][j] = 0
    j+=1
i+=1
# set the initial value k for while loop, maximum number of iterations 
#and tolerance
k = 0
limit = 10
TOL = 1e-05

while k < limit:
#create the list B that same as matrix size nx1 
    B = [f(x[i],p[i],(p[i+1]-p[i-1])/(2*h))*h**2 for i in range(1,n+1)]
    B[ 0] -= alpha # input the alpha value in array B at the first position 
    B[-1] -= beta  #input the beta value in array B at the last position
    
    # Solve the matrix equation where Au = B for u 
    u = list(solve(A, B))
     #input the alpha and beta values in first and last position of array u
    u = [alpha] + u + [beta]    #list values of y_approximation
# append the values of u into array p1 and set p as u 
    p1.append(u)
    p = u
    k += 1      #update the value of k
    #set the condition whether the approximation difference less than TOL
    if all(p1[k][i]-p1[k-1][i] < TOL for i in range(len(p))):
        break

#print the results
print('')
print('n = {} and h = {}'.format(n,h))  #print the number of iterations and step size
print('Number of iterartions : ',k) #print the number of iterations in while loop
print('{:12}{:12}{:12}'.format('x','u(x)','error')) #print xi,u(xi) = u and error = u - y_approximation
for i in range(0,n+2):
    print('{:.4f}{:12.6f}{:12.6f}'.format(x[i],u[i],abs(u[i]-yexact[i])))

#plot the graph for exact and approximation solution
xtrue = np.arange(0,1.01,0.01)
plt.plot(xtrue,uexact(xtrue),'-g') #exact solution
plt.plot(x,u,'ok')  #approximation solution
plt.title('Non Linear Finite Method')
plt.xlabel('x')
plt.ylabel('u(x)')






