# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 14:38:50 2020

@author: Ahmad Aiman B. Mohd Nazir
"""
#import the libraries
import numpy as np
import matplotlib.pyplot as plt
#define the exact solution u(x)
def uexact(x):
    return np.log((np.e - 1)*x + 1)
#define the function f = u", f1 = df/du and f2 = df/du' 
def f(x,u,uprime):
    return -uprime**2

def f1(x,u,uprime):
    return 0

def f2(x,u,uprime):
    return -2*uprime
#set the boundary condition u(x0) = alpha and u(xn) = beta
x0 = 0
alpha = 0
xn =1 
beta = 1

#calculate the step size/increment
n = 20      #number of iterations
h = (xn - x0)/n

#calculate the initial guess value, s0
s = (beta - alpha)/(xn-x0)
#calculate and create an array of value for xi
xtrue = []
for i in range(0,n+1):
    xtrue.append(x0 +i*h)
    i+=1
#calculate and create an array value for uexact(xi)    
yexact = []
for val in xtrue:
    yexact.append(uexact(val))

#initilize value for while loop and set the tolerance 
k = 1
TOL = 1e-5

#RK4 algorithm to solve non linear eqn
while k <= n:
    y1 = [alpha]        #define y1 as list of y_approximation
    y2 = [s]            #define y2 as list of initial guess
    y3 = 0
    y4 = 1
    
    for i in range(0,n):
        x = x0 +h*(i)
        
        k11 = h*y2[i]
        k12 = h*f(x,y1[i],y2[i])
        
        k21 = h*(y2[i] + k12/2)
        k22 = h*f(x+h/2,y1[i]+k11/2,y2[i]+k12/2)
        
        k31 = h*(y2[i]+k22/2)
        k32 = h*f(x+h/2,y1[i]+k21/2,y2[i]+k22/2)
        
        k41 = h*(y2[i]+k32)
        k42 = h*f(x+h,y1[i]+k31,y2[i]+k32)
        
        y1.append(y1[i] + (k11 +2*k21 + 2*k31 + k41)/6)
        y2.append(y2[i] + (k12 + 2*k22 + 2*k32 + k42)/6)
        
        j11 = h*y4
        j12 = h*(y3*f1(x,y1[i],y2[i]) + y4*f2(x,y1[i],y2[i]))
        
        j21 = h*(y4+j12/2)
        j22 = h*((y3+j11/2)*f1(x+h/2,y1[i],y2[i]) + (y4+j12/2)*f2(x+h/2,y1[i],y2[i]))
        
        j31 = h*(y4+j22/2)
        j32 = h*((y3+j21/2)*f1(x+h/2,y1[i],y2[i]) + (y4+j22/2)*f2(x+h/2,y1[i],y2[i]))
        
        j41 = h*(y4+j32)
        j42 = h*((y3+j31)*f1(x+h/2,y1[i],y2[i]) + (y4+j32)*f2(x+h/2,y1[i],y2[i]))

        y3 +=(j11 + 2*j21+ 2*j31 + j41)/6
        y4 +=(j12 + 2*j22+ 2*j32 + j42)/6
        
    if abs(y1[n]-beta) <= TOL:      #check the y_approximation wheather close to beta
        break
    else:                           #If not continue with the next iteration to improve the s value
        s -= (y1[n]-beta)/y3
        k+=1
#print the results
print('n = {} and h = {}'.format(n,h)) #print number of iteration with step size
print('Number of iterations : ',k)  #print number of iteration in while loop
#print the result where x = xtrue array, u(x) = y1 array and error = y1 - yexact 
print('{:12}{:12}{:12}'.format('x','u(x)','error'))
for i in range(n+1):
    print('{:.4f}{:12.6f}{:12.6f}'.format(xtrue[i],y1[i],abs(y1[i]-yexact[i])))

#plot the exact solution with its approximation
xactual = np.arange(x0,xn+0.01,0.01)
plt.plot(xactual,uexact(xactual),'-g') #exact solution
plt.plot(xtrue,y1,'ok') #approximation solution
plt.title('Non linear Shooting Method')
plt.xlabel('x')
plt.ylabel('u(x)')



        
