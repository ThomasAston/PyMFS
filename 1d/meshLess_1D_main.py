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
resolution = 1000
dx = l/resolution
X = np.linspace(0,1,resolution)

u_analytical = []
for i in range(0,len(X)):
    u_analytical.append(0.5*(X[i]-(X[i]**3)/3)+2*X[i]+1)

# Generate shape functions at each node...
# Start by initialising empty phi array, which includes rows for each node
phi = np.zeros((len(x),resolution))
W_i = np.zeros((len(x),resolution))
h = np.zeros((len(x),2,resolution))

# print(W_i[5][100])
for i in range(N):
    print('Current node: ', i+1)
    for a in range(len(X)):         
        # Calculate the normalised distance from node i
        s = np.abs(X[a] - x[i])/r

        # Calculate W_i
        if s <= 1:
            W_i[i][a] = 1-6*s**2+8*s**3-3*s**4
        else:
            W_i[i][a] = 0

        # Calculate W_j for each of the other nodes
        sum_W_j =0
        for j in range(N):
            s_j = np.abs(X[a] - x[j])/r
            
            if s_j<=1:
                W_j = 1-6*s_j**2+8*s_j**3-3*s_j**4
            else:
                W_j = 0

            sum_W_j += W_j

        phi[i][a] = W_i[i][a]/sum_W_j
        
        # Compute the 2 shape functions
        h[i][0][a] = phi[i][a]
        h[i][1][a] = phi[i][a] * (X[a] - x[i])/r

# Differentiate the shape functions
dh_dx = np.zeros((len(x),2,resolution))
for i in range(N):
    for j in range(2):
        for a in range(resolution-1):
            dh_dx[i][j][a] = (h[i][j][a+1] - h[i][j][a])/(X[a+1]-X[a])


# Initialise the stiffness matrix...
K = np.zeros((N,N,2,2))

# Loop over each node and each degree of freedom...
for i in range(N):
    for j in range(N):
        for m in range(2):
            for n in range(2):    
                
                # Loop from x1 to x2 to numerically integrate...
                for a in range(0, resolution):
                    if a == 0 or a==resolution:#(resolution/(N-1)):
                        f_x = dh_dx[i][m][a] * dh_dx[j][n][a]
                    else:
                        f_x = 2*dh_dx[i][m][a] * dh_dx[j][n][a]

                    K[i][j][m][n] += dx/2*f_x

# Force vector...
f_Im = np.zeros((N, 2))
fHat_Im = np.zeros((N,2))
KU_ImJn = np.zeros((N,N,2,2))

for i in range(N):
    for m in range(2):
        
        f_Im[i][m] = np.trapz(X[:]*h[i][m][:], x=None, dx=dx)

        # Evaluate fHat_Im 
        if i==5: # node is on neumann boundary
            fHat_Im[i][m] = f_s * h[i][m][resolution-1]
        elif i==0: # node is on dirichlet boundary
            sum_KU_ImJn = 0
            for j in range(N):
                for n in range(2):
                    KU_ImJn[i][j][m][n] = -((h[i][m][1]*h[j][n][1] - h[i][m][0]*h[j][n][0]))/(X[1]-X[0])
                    # sum_KU_ImJn += KU_ImJn
                    # print(sum_KU_ImJn)
            
            fU_Im = u_s*dh_dx[i][m][0]

            fHat_Im[i][m] = sum_KU_ImJn - fU_Im
        else: 
            fHat_Im[i][m] = 0


f = f_Im + fHat_Im
print(f)
# DoF vector...

# q = 

# for i in range(N):

# Plot the shape functions
ax = plt.axes()
plt.plot(X,h[0][0][:])
plt.plot(X,h[0][1][:])
plt.plot(X,h[1][0][:])
plt.plot(X,h[1][1][:])
plt.plot(X,h[2][0][:])
plt.plot(X,h[2][1][:])
plt.plot(X,h[3][0][:])
plt.plot(X,h[3][1][:])
# plt.plot(X,h[5][0][:])
# plt.plot(X,h[5][0][:])
# plt.plot(X,h[5][0][:])
plt.xlabel("x")
plt.ylabel("Shape functions")
plt.show()

# Stiffness matrix


# Post-processing...

# ax = plt.axes()
# plt.plot(X,u_analytical)
# plt.xlabel("x")
# plt.ylabel("u(x)")

# plt.show()