# meshLess_1D v1

# Author(s): Thomas Aston
# Date updated: 27/01/2022

# Main file for calling functions for solving 1D problems using MFS

# To do list:
# 
# 
#

import numpy as np

import sympy as sym
from sympy import *
from sympy.plotting import plot

from matplotlib import pyplot as plt
init_printing( use_latex='mathjax' )  # use pretty math output
# Set up the problem...

l = 1  # Domain length
N = 6   # Number of nodes 
r = 0.2 # Sphere radius

X = np.linspace(0, l, N) # Generate nodes

print('List of nodal coordinates: ',X)

u_s = 1 # Essential boundary at x = 0

f_s = 2 # Natural boundary at x = 1

# Analytical solution...
x = sym.Symbol('x')

u_analytical = 0.5*(x-x**3/3)+2*x+1

# print(u_analytical)
p1 = plot(u_analytical,xlim=(0,1), ylim=(0,3.5), show=False)

# Generate the shape functions...
phi = []
x_i = 0

W_i = []
sum_W_j = []
phi = []
h = []
for i in range(N):
    s_i = ((x - x_i)**2)**(0.5)/r

    W_i.append(sym.Piecewise((1 - 6*s_i**2 + 8*s_i**3 - 3*s_i**4, s_i<=1) , (0, s_i>1)))
    
    x_j = 0
    W_j = []
    for j in range(N):        
        s_j = ((x - x_j)**2)**(0.5)/r
        W_j.append(sym.Piecewise((1 - 6*s_j**2 + 8*s_j**3 - 3*s_j**4, s_j<=1) , (0, s_j>1)))
        x_j+=1/(N-1)
    
    sum_W_j.append(np.sum(W_j))
    
    phi.append(W_i[i]/sum_W_j[i])

    h.append(phi[i])
    h.append(phi[i]*(x-x_i)/r)

    x_i+=l/(N-1)

h = np.reshape(h, (N, 2))

# p2 = plot(h[1,0],h[1,1],xlim=(0,1),show=True)

# Assemble the stiffness matrix...

# k = np.zeros([2,2])
K = np.zeros([N*2, N*2])
# print(K)

print(diff(h[0,0],x))
K[0,0] = integrate(diff(h[0,0],x)*diff(h[0,0],x),(x,0,1))

print(K)

# plt.show()