This code was written by Ahmad Aiman Mohd Nazir (17060060) for computational physics course mini project.

In this code we are going to solve for nonhomogenous linear ODEs where 
it genaral form can be written as follow <latex> $my'' + by' + ky = f(x)$</latex>
for this ODE, the physics problem that we are going to solve is related to 
forced vibration on the spring.  

Linear ODEs is <latex>$y'' + 8y' +16y = 8sin(4x)$</latex>.  
The exact solution of this ODEs is <latex>$y(x)$ = $\frac{1}{4}$$e^{-4x}$ + $x$$e^{-4x}$ -$\frac{1}{4}$$cos(4x)$</latex>

The method that we are going to use here is finite difference method. 
The general expression of ODEs can be represented as <latex>$x'' + P(x)x' + q(x)x = f(x)$</latex>