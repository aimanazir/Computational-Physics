# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 22:46:11 2021

@author: Ahmad Aiman Mohd Nazir
"""

import math
import cmath
import matplotlib.pyplot as mpl

hbar = 1.0 # Reduced Planck's constant
m    = 1.0 # Mass
k0   = 1.0 # Initial wavenumber

# Grid and time intervals
dx   = 0.25/k0
dt   = 0.25*m/(hbar*k0**2)
tmax = 200.0*m/(hbar*k0**2)
xmin = -200.0/k0
xmax = 200.0/k0
N = int((xmax-xmin)/dx)

# Initial data
x  = [0.0]*(N+1)
u  = [0.0]*(N+1)
v  = [0.0]*(N+1)
p  = [0.0]*(N+1)
x0 = -100.0/k0
sigma = 10.0/k0

for j in range(N+1):
    x[j] = xmin+j*dx
    u[j] = cmath.exp(-(x[j]-x0)**2/(4.0*sigma**2)+1j*k0*x[j])
    u[j] = u[j]/(2.0*math.pi*sigma**2)**0.25
    p[j] = u[j].real**2+u[j].imag**2

# Potential
E0 = (hbar*k0)**2/(2.0*m)
a = 5.0/k0
V = [0.0]*(N+1)
for j in range(N+1):
    if abs(x[j]) < a: V[j] = E0

# Setup coefficients of the tridiagonal matrix
alpha = gamma = -1j*hbar*dt/(4.0*m*dx**2)
beta = [0.0]*(N+1)
for j in range(N):
    beta[j] = 1.0 - 2.0*alpha + 1j*(V[j]/(2.0*hbar))*dt

# Prepare animated plot
mpl.ion()
fig = mpl.figure()
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
ax2.plot(x, [Vj/E0 for Vj in V], 'b')
(line, ) = ax1.plot(x, p, 'k-')

# Preform the evolution
t = 0.0
while (t-tmax) < 0.5*dt:
    # Update plot
    for j in range(N+1):
        p[j] = u[j].real**2+u[j].imag**2
    line.set_ydata(p)
    mpl.title(f't = {t:5f}')
    mpl.draw()
    mpl.pause(0.01)

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
mpl.ioff()
mpl.draw()
mpl.show()