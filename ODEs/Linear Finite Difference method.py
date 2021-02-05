# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 13:19:09 2020

@author: Ahmad Aiman Mohd Nazir
"""
#import the relevent libraries
from scipy.linalg import solve
import  numpy as np
import matplotlib.pyplot as plt
x = [0]         #set the initial value of x and define all parameters
f = []          #set array of f for f(x)
pi = 3.142      #define pi value
n = 8           #the number of step size that we want       
h = 1/n         #step size or increment
p = 0           #p(x)
q = pow(np.pi,2)    #q(x)
u = 1               #u(0)
u_1 = -1            #u(1)
#define matrix A with size (n-1)x(n-1)
A = np.zeros(((n-1,n-1)))
#append x values into array x by using loop
for i in range (0,n):
    x.append(x[i]+h)
    i+=1
print('x = {}'.format(x))    #print list of x
print('')

#append the values of f(x) for each value of x[j]
for j in range(0,n+1):
    f.append(2*pow(np.pi,2)*np.cos(np.pi*x[j])) 
    j += j
print('f = ')
print(f)
print('')
#calculate the value of a,b and c
a = 2+ q*h**2
b = -(1-1/2*p*h)
c = -(1+1/2*p*h)
print('a = {}'.format(a))
print('b = {}'.format(b))
print('c = {}'.format(c))
print('')

#input the values obtain from a,b and c into matrix A
for i in range(0,n-1):
    for j in range(0,n-1):
        if i == j:
            A[i][j] = a
        elif i<j and j-i==1:
            A[i][j] = b
        elif i>j and i-j==1:
            A[i][j] = c
        else:
            A[i][j] = 0
    j+=1
i+=1
print('matrix A =')         #print matrix A
print(A)
print('')
#input the values into matrix Z
Z = np.zeros((n-1,1))
for i in range(0,n-1):
    for j in range(0,1):
        if i == j:
            Z[i][j] = (1+(1/2)*p*h)*u + h**2*f[i+1]
        elif i > j and i<=n-3:
            Z[i][j] = h**2*f[i+1]
        elif i>j and i-j == n-2:
            Z[i][j] = (1-(1/2)*p*h)*u_1+(h**2)*f[i+1]
    j+=1
i+=1
print('matrix Z =')         #print marix Z
print(Z)
print('')
X = solve(A,Z)               #solve matrix AX = Z
print('X =')
print(X)
print('')
xnew = np.arange(0,1.01,h)      #define new array for x, ranging from 0 to 1 with increment = h
print('xnew = {}'.format(xnew)) #print values xnew
#y = X.transpose()
y1 = np.insert(X,0,1)            #insert value y = 1 into y1[0]
ynew = np.append(y1,-1)          #append value y = -1 in y1 array
print('ynew = {}'.format(ynew))
#plot the graph ynew vs xnew
plt.plot(xnew,ynew,'ok')
plt.xlabel('x')
plt.ylabel('U(x)')
plt.title('U(x) versus x for finite difference method')
plt.show()

