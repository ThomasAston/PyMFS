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
from matplotlib import pyplot as plt

# Set up the problem...

l = 1  # Domain length
N = 6   # Number of nodes 
r = 0.2 # Sphere radius

x = np.linspace(0, l, N) # Generate nodes

print('List of nodal coordinates: ',x)

u_s = 1 # Essential boundary at x = 0

f_s = 2 # Natural boundary at x = 1

# Analytical solution...
resolution = 100
dx = l/resolution
X = np.linspace(0,1,resolution)

u_analytical = []
for i in range(0,len(X)):
    u_analytical.append(0.5*(X[i]-(X[i]**3)/3)+2*X[i]+1)


# Generate shape functions at each node...

# Start by initialising empty phi array, which includes rows for each node
phi = np.zeros([len(x),resolution])
W_i = np.zeros([len(x),resolution])
W_j = np.zeros([len(x),resolution])
h = np.zeros([len(x),2,resolution])

for i in range(N):
    print('Current node: ', i+1)
    for a in range(len(X)):         
        # Calculate the normalised distance from node i
        s = np.abs(X[a] - x[i])/r

        # Calculate W_i
        if s <= 1:
            W_i[i,a] = 1-6*s**2+8*s**3-3*s**4
        else:
            W_i[i,a] = 0

        # Calculate W_j for each of the other nodes
        for j in range(0,N):
            s_j = np.abs(X[a] - x[j])/r
            
            if s_j<=1:
                W_j[j,a] = 1-6*s_j**2+8*s_j**3-3*s_j**4
            else:
                W_j[j,a] = 0

        phi[i,a] = W_i[i,a]/sum(W_j[:,a])
        
        # Compute the 2 shape functions
        h[i,0,a] = phi[i,a]
        h[i,1,a] = phi[i,a] * (X[a] - x[i])/r



# Plot the shape functions
ax = plt.axes()
plt.plot(X,h[0,0,:])
plt.plot(X,h[0,1,:])
plt.plot(X,h[1,0,:])
plt.plot(X,h[1,1,:])
plt.plot(X,h[2,0,:])
plt.plot(X,h[2,1,:])
plt.plot(X,h[3,0,:])
plt.plot(X,h[3,1,:])
plt.plot(X,h[4,0,:])
plt.plot(X,h[4,1,:])
plt.plot(X,h[5,0,:])
plt.plot(X,h[5,1,:])

plt.xlabel("x")
plt.ylabel("Shape functions")
plt.show()

# Stiffness matrix
K = np.zeros([N,2,N])


# Post-processing...

# ax = plt.axes()
# plt.plot(X,u_analytical)
# plt.xlabel("x")
# plt.ylabel("u(x)")

# plt.show()