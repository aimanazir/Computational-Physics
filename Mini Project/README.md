This code was written by Ahmad Aiman Mohd Nazir (17060060) for computational physics course mini project.

In this code we are going to solve for nonhomogenous linear ODEs where 
it genaral form can be written as follow <img src="https://render.githubusercontent.com/render/math?math=my'' %2B by' %2B ky = f(x)"> and
for this ODE, the physics problem that we are going to solve is related to 
forced vibration on the spring.  

Linear ODEs is <img src="https://render.githubusercontent.com/render/math?math=y'' %2B 8y' %2B 16y = 8sin(4x)">.   
The exact solution of this ODEs is <img src="https://render.githubusercontent.com/render/math?math=$y(x) = \frac{1}{4}e^{-4x} %2B xe^{-4x} -\frac{1}{4}cos(4x)$">

The method that we are going to use here is finite difference method. 
The general expression of ODEs can be represented as <img src="https://render.githubusercontent.com/render/math?math=x'' %2B P(x)x' %2B q(x)x = f(x)">
