# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 16:28:11 2020

@author: Ahmad Aiman Mohd Nazir 
"""

import math
import matplotlib.pyplot as plt

hbar = 1.0 # Reduced Planck's constant
m    = 1.0 # Mass
k    = 1.0 # Spring constant

# Grid and time intervals
dx   = 0.01 #spatial time intervals
dt   = 0.05 #time step size
tmax = 5.0 #max time range  #Here I change for grid x and t
xmin = -3.0 #staring point for x    
xmax = 3.0 #final point for x
N    = int((xmax-xmin)/dx)  #number of mesh points

# Initial data
x = [0.0]*(N+1) #initial x for N+1 data points
u = [0.0]*(N+1) #initial u for N+1 data points
v = [0.0]*(N+1) #initial v for N+1 data points
p = [0.0]*(N+1) #initial p for N+1 data points

omega0 = (k/m)**0.5 #calculate the angular frequency 
sigma  = (hbar/(2.0*m*omega0))**0.5 #1/coeffiecient 

for j in range(N+1):
    x[j] = xmin+j*dx #x mesh points
    u[j] = x[j]*math.exp(-x[j]**2/(4.0*sigma**2))
    u[j] = math.sqrt(2*m*omega0/(hbar))*u[j]/(2.0*math.pi*sigma**2)**0.25  #calculate for eigenfunctions at n =1 
    #Here I change for question no.3(i) where n =1
    p[j] = u[j].real**2+u[j].imag**2 #prob density 

# Potential
E0 = 1.5*hbar*omega0 #ground state, n=1 energy H0 
#Change the value for eigenvalue calculation
V  = [0.0]*(N+1) #initial potential for each mesh point

for j in range(N+1):
    V[j] = 0.5*k*x[j]**2 #calculate potential at each mesh points 

#Crank-Nicolson method 
# Setup coefficients of the tridiagonal matrix
alpha = gamma = -1j*hbar*dt/(4.0*m*dx**2) #assign value of alpha and gamma coeffiecient
beta = [0.0]*(N+1)
for j in range(N):
    beta[j] = 1.0 - 2.0*alpha + 1j*(V[j]/(2.0*hbar))*dt #calculate the beta coeffiecient

# Prepare animated plot
plt.ion()
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_xlim(xmin, xmax)
ax1.set_ylim(0.0, 1.1*max(p))
ax1.set_xlabel('x')
ax1.set_ylabel('Probability density')
ax2 = ax1.twinx()
ax2.set_xlim(xmin, xmax)
ax2.set_ylim(0.0, 1.1*max(V)/E0)
ax2.set_ylabel('V / E0')

# Plot potential function and wave function
ax2.plot(x, [Vj/E0 for Vj in V], 'b')   #blue line on the plot for potential and x, for Q3(ii)
(line, ) = ax1.plot(x, p, 'k-') #black line on plot for prob density and x, for Q3(iii)

# Preform the evolution
t = 0.0
while t-tmax < 0.5*dt:
    # Update plot
    for j in range(N+1):
        p[j] = u[j].real**2+u[j].imag**2
    line.set_ydata(p)
    plt.title(f't = {t:5f}') 
    plt.draw()
    plt.pause(0.1)

    # Set the values of the RHS
    for j in range(1, N):
        v[j] = -alpha*u[j-1]+(2.0-beta[j])*u[j]-gamma*u[j+1]
    v[1] -= alpha*u[0]
    v[N-1] -= gamma*u[N]

    # Forward sweep
    u[1] = v[1]/beta[1]
    v[1] = gamma/beta[1]
    for j in range(2, N):
        den = beta[j]-alpha*v[j-1]
        u[j] = (v[j]-alpha*u[j-1])/den
        v[j] = gamma/den

    # Backward sweep
    for j in reversed(range(1, N)):
        u[j] -= u[j+1]*v[j]
    t += dt

# Freeze final plot
plt.ioff()
plt.draw()
plt.show()





