# -*- coding: utf-8 -*-
"""
Created on Sun Nov  8 11:52:13 2020

@author: Ahmad Aiman Mohd Nazir
"""
# Solving ODE linear BVP using shooting method
# import the relevence libraries
import numpy as np
import matplotlib.pyplot as plt

# Define the function f1,f2,f3 and f4
def f1(x,y1,y2):
    return y2
def f2(x,y1,y2):
    return y1*np.pi**2 - (2*np.pi**2)*np.cos(x*np.pi)
def f3(x,y3,y4):
    return y4
def f4(x,y3,y4):
    return y3*np.pi**2

# set the boundary condition
x0 = 0
xn = 1
alpha = 1
beta = -1

# choose the number of grid points
n = 8
h = (xn-x0)/n

# create a list number of x from x0,x1,...,xn
x = []
for i in range(0,n+1):
    x.append(x0+i*h)
    i+=1

#print list of values x
print('x =',x)
print('')   #for spacing
# define list1(y1i), list2(y2i), list3(y3i) and list4(y4i)
list1 = []
list2 = []
list3 = []
list4 = []
y1i = 1         #set the initial values for each iteration
y2i = 0         # y1i = y1(0) , y2i = y2(0), y3i = y3(0), y4i = y4(0)
y3i = 0 
y4i = 1

print('{:10}{:10}{:10}{:10}{:10}'.format('x0','y1i','y2i','y3i','y4i'))

#RK4 loop for solve two coupled IVPs
for i in range(0,n+1):
        K11 = h*f1(x0,y1i,y2i)
        K12 = h*f2(x0,y1i,y2i)
        K13 = h*f3(x0,y3i,y4i)
        K14 = h*f4(x0,y3i,y4i)
    
        K21 = h*f1(x0+0.5*h, y1i+0.5*K11, y2i+0.5*K12)
        K22 = h*f2(x0+0.5*h, y1i+0.5*K11, y2i+0.5*K12)
        K23 = h*f3(x0+0.5*h, y3i+0.5*K13, y4i+0.5*K14)
        K24 = h*f4(x0+0.5*h, y3i+0.5*K13, y4i+0.5*K14)
        
        K31 = h*f1(x0+0.5*h, y1i+0.5*K21, y2i+0.5*K22)
        K32 = h*f2(x0+0.5*h, y1i+0.5*K21, y2i+0.5*K22)
        K33 = h*f3(x0+0.5*h, y3i+0.5*K23, y4i+0.5*K24)
        K34 = h*f4(x0+0.5*h, y3i+0.5*K23, y4i+0.5*K24)
        
        K41 = h*f1(x0+h, y1i+K31, y2i+K32)
        K42 = h*f2(x0+h, y1i+K31, y2i+K32)
        K43 = h*f3(x0+h, y3i+K33, y4i+K34)
        K44 = h*f4(x0+h, y3i+K33, y4i+K34)
        
        list1.append(y1i)           #add the values into the list
        list2.append(y2i)
        list3.append(y3i)
        list4.append(y4i)
        
        print('{:.4f}{:10.6f}{:10.6f}{:10.6f}{:10.6f}'.format(x0,y1i,y2i,y3i,y4i))
        y1i += (K11+2*K21+2*K31+K41)/6          #update the values for y1i,y2i,y3i and y4i
        y2i += (K12+2*K22+2*K32+K42)/6
        y3i += (K13+2*K23+2*K33+K43)/6
        y4i += (K14+2*K24+2*K34+K44)/6
        
        x0 += h                 #increment of x0 or update new value of x0


#obtain the values for y1(b), y3(b) and print it
y1n = list1[-1]              #change the value in [] according to n 
y3n = list3[-1]
print('')
print('y1(b) = {} and y3(b) = {}'.format(y1n,y3n))
print('')
#calculate the coeffiecent of y3(x)
c = (beta-y1n)/y3n

#obtain the approximate solution of BVP, u(x) = y1(x) + c*y3(x)
u = []
for i in range(0,n+1):
    u.append(list1[i]+c*list3[i])
    i+=1

#print the result of u(x) and x
print('{:13}{}'.format('x','u(x)'))
for i in range(0,n+1):
    print('{:.6f}{:12.6f}'.format(x[i],u[i]))
    

#plot the graph 
plt.title('Linear Shooting Method for BVP')
plt.xlabel('x')
plt.ylabel('u(x)')
plt.plot(x,u,'ok')
plt.show()



    
    


