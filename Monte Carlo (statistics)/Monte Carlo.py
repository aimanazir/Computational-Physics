# -*- coding: utf-8 -*-
"""
Created on Tue Jan  5 18:33:03 2021

@author: Ahmad Aiman Mohd Nazir
"""
#import the relevent libraries 
from sympy import integrate, sqrt
from sympy.abc import x
from scipy.integrate import simps
import random
import numpy as np

#define the function
def f(x):
  return 1/np.sqrt(1+x**4)

#monte carlo calculation
n = 10000
sum1 = 0.0
sum2 = 0.0
exactvalue = 0.92703733865069
for i in range(1,n):
  xm = random.uniform(0,1)
  sum1 = sum1 + f(xm)
  sum2 = sum2 + f(xm)**2

s = sum1/n
ds = np.sqrt(abs(sum2/n - (sum1/n)**2)/n)
error = abs(s - exactvalue)/exactvalue

print('Monte Carlo Calculation')
print('Exact value = {}'.format(exactvalue))
print('Monte Carlo value = {}'.format(s))
print('Error = {}'.format(error))
print('')

#Symbolic calculation
a = integrate((1/sqrt(1+x**4)),(x,0,1))
print('Symbolic Method')
print(float(a))
print('')

#Simpson's method
xs = np.arange(0,1,0.01)
y = f(xs)
area = simps(y,x=xs)
print('Simpson Method')
print(area)
