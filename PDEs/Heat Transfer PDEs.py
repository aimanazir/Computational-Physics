# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 21:36:31 2020

@author: Ahmad Aiman Mohd Nazir
"""

#import numpy as np
import matplotlib.pyplot as mpl


#fixed parameters
T0 = 25.0
D = 1e-4
L = 1.0

#input parameters
dx = L*float(input('grid spacing in units of wire lenght (L) ->'))
dt = dx**2/D*float(input('time step in units of (dx/c) ->'))
tmax = L**2/D*float(input('evolution time in units of (L/c) ->'))

#construct initial data
N = int(L/dx)
x = [0.0]*(N+1)
u0 =[0.0]*(N+1)
v0 = [0.0]*(N+1)
u1 =[0.0]*(N+1)

for j in range(N+1):
   x[j] = j*dx


#prepare animated plot
mpl.ion()
(line,)=mpl.plot(x,u0,'-k')
mpl.ylim(0,T0)
mpl.xlabel('x(m)')
mpl.ylabel('Temperature (Celcius)')

#perform evolution
t = 0.0
while t < tmax:
     #update plot
     line.set_ydata(u0)
     mpl.title('t=%5f'%t)
     mpl.draw()
     mpl.pause(0.1)
     #derivatives at interior points
     for j in range(1,N):
         u1[j] = u0[j]+dt*D*(u0[j+1]-2.0*u0[j]+u0[j-1])/dx**2

     #boundary conditions
     u1[0]=T0
     u1[N]=0.0

     #swap old and new lists
     (u0,u1)=(u1,u0)
     t += dt

#freeze final plot
mpl.ioff()
mpl.show()