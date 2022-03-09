# meshLess_1D v1

# Author(s): Thomas Aston
# Date updated: 08/02/2022

# Main file for solving 1D axial bar problem using MFS

# To do list:
# 
# 
#

import numpy as np
import sympy as sym
from matplotlib import pyplot as plt

##############################################
# PROBLEM DEFINITION
##############################################

# We want to solve Poisson's equation for a 1D bar subjected to distributed axial
# loading. Analytical displacement field is chosen such that the following domain
# properties and boundary values are valid.

### Domain ###
l = 1  # Domain length

### Discretisation ###
N = 6   # Number of nodes 
r = 0.2 # Sphere radius
x = np.linspace(0, l, N) # Generate nodes
print('List of nodal coordinates: ',x)
order = 2 # Number of terms in the local basis (= number of required shape functions)

### BCs ###
u_s = 1 # Essential boundary at x = 0
f_s = 2 # Natural boundary at x = 1

### Analytical displacement field ###
resolution = 1001
dx = l/resolution
X = np.linspace(0,1,resolution)
u_analytical = []
for i in range(0,len(X)):
    u_analytical.append(0.5*(X[i]-(X[i]**3)/3)+2*X[i]+1)


##############################################
# GENERATE SHAPE FUNCTIONS
##############################################

### Initialise empty arrays ###
phi = np.zeros((len(x),resolution)) # Shepard PU functions
W_i = np.zeros((len(x),resolution)) # Weight functions
h = np.zeros((len(x),order,resolution)) # Shape functions

sum_W_j = np.zeros(len(X))
### Generate by looping over nodes ###
for i in range(N): 
    # print('Current node:', x[i])
    for a in range(len(X)):         
        # Calculate the normalised distance from node i 
        s = np.abs(X[a] - x[i])/r
    
        # Calculate W_i
        if s <= 1:
            W_i[i][a] = 1-6*s**2+8*s**3-3*s**4
        else:
            W_i[i][a] = 0

        # Calculate W_j for each of the other nodes
        sum_W_j_a =0
        for j in range(N):
            s_j = np.abs(X[a] - x[j])/r
            
            if s_j<=1:
                W_j = 1-6*s_j**2+8*s_j**3-3*s_j**4
            else:
                W_j = 0

            sum_W_j_a += W_j

        sum_W_j[a] = sum_W_j_a
    
        phi[i][a] = W_i[i][a]/sum_W_j[a]

        # Compute the 2 shape functions
        if order == 2:
            h[i][0][a] = phi[i][a]
            h[i][1][a] = phi[i][a] * (X[a] - x[i])/r
    

# Differentiate the shape functions
dh_dx = np.zeros((len(x),order,resolution))
for i in range(N):
    for j in range(order):
        for a in range(resolution):
            if a == resolution-1:
                dh_dx[i][j][a] = (h[i][j][a] - h[i][j][a-1])/(X[a]-X[a-1])
            else:
                dh_dx[i][j][a] = (h[i][j][a+1] - h[i][j][a])/(X[a+1]-X[a])


##############################################
# ASSEMBLE STIFFNESS MATRIX
##############################################

### Initialise empty arrays ###
K_nested = np.zeros((N,N,order,order)) # Stiffness matrix
K_open = np.zeros((N*order,N*order)) # Unravelled stiffness matrix

### Loop over nodes and DoFs ###
for i in range(N):
    for j in range(N):
        for m in range(order):
            for n in range(order):    
                
                if x[i]-r<0:
                    x1 = 0
                    x2 = x[i]+r
                elif x[i]+r>l:
                    x1 = x[i]-r
                    x2 = l
                else:
                    x1 = x[i]-r
                    x2 = x[i]+r
                
                # x1 = round(x1,2)
                # x2 = round(x2,2)
   
                x1_loc = int(x1/l * (resolution-1))
                x2_loc = int(x2/l * (resolution-1))

                K_nested[i][j][m][n] = np.trapz(dh_dx[i,m,x1_loc:x2_loc] * dh_dx[j,n,x1_loc:x2_loc], x=None, dx=dx)
                K_open[order*i:order+order*i,order*j:order+order*j] = K_nested[i,j,:,:]


##############################################
# IMPLEMENT BOUNDARY CONDITIONS
##############################################

f_Im = np.zeros((N, order))
fHat_Im = np.zeros((N,order))
KU_ImJn = np.zeros((N,N,order,order))
KU_ImJn_open = np.zeros((N*order,N*order))

for i in range(N):
    for m in range(order):
        
        if x[i]-r<0:
            x1 = 0
            x2 = x[i]+r
        elif x[i]+r>l:
            x1 = x[i]-r
            x2 = l
        else:
            x1 = x[i]-r
            x2 = x[i]+r
        
        x1 = round(x1,2)
        x2 = round(x2,2)

        x1_loc = int(x1/l * (resolution-1))
        x2_loc = int(x2/l * (resolution-1))

        f_Im[i][m] = np.trapz(X[x1_loc:x2_loc]*h[i,m,x1_loc:x2_loc], x=None, dx=dx)

        # Evaluate fHat_Im 
        if x[i] + r > l: # sphere intersects neumann boundary
            fHat_Im[i][m] = f_s * h[i][m][resolution-1]
        elif x[i] - r < 0: # sphere intersects dirichlet boundary
            fU_Im = -u_s*dh_dx[i][m][0]
            fHat_Im[i][m] =  -fU_Im
            
            for j in range(N):
                for n in range(order):
                    KU_ImJn[i][j][m][n] = -((h[i][m][1]*h[j][n][1] - h[i][m][0]*h[j][n][0]))/(X[1]-X[0])
                    KU_ImJn_open[order*i:order+order*i,order*j:order+order*j] = KU_ImJn[i,j,:,:]
        else: 
            fHat_Im[i][m] = 0

##############################################
# SOLVE SYSTEM OF EQUATIONS
##############################################
K_open = K_open-KU_ImJn_open

f = f_Im + fHat_Im
print(f)
f = np.ravel(f)

q = np.linalg.solve(K_open,f)

##############################################
# REASSEMBLE DISPLACEMENT FIELD
##############################################
q = np.reshape(q, [N,2])
print(q)

u = np.zeros([resolution,1])

for i in range(resolution):
    for j in range(N):
        for n in range(order):
            u[i] += h[j][n][i]*q[j][n]
  

# print(u)
##############################################
# PLOT RESULTS
##############################################
ax2 = plt.axes
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.title('Shape functions')

plt.plot(X,h[0,0,:], label='h_00')
# plt.plot(X,dh_dx[0,0,:]*dh_dx[1,0,:], label='h_00')
plt.plot(X,h[0,1,:], label='h_01')
plt.plot(X,h[1,0,:], label='h_10')
# plt.plot(X,dh_dx[3,1,:], label='h_10')

plt.plot(X,h[1,1,:], label='h_11')
plt.plot(X,h[2,0,:], label='h_20')
# plt.plot(X,h[2,1,:], label='h_21')
# plt.plot(X,h[3,0,:], label='h_30')
# plt.plot(X,h[3,1,:], label='h_31')
# plt.plot(X,h[4,0,:], label='h_40')
# plt.plot(X,h[5,0,:], label='h_50')
# plt.plot(X,sum_W_j, label='sum_W_j')


plt.legend()
plt.xlabel("x")
plt.ylabel("h")
plt.show()


ax3 = plt.axes
plt.rcParams["font.family"] = "serif"
plt.rcParams["mathtext.fontset"] = "dejavuserif"
plt.title('Axial bar deformation')
plt.plot(X,u_analytical,label='Analytical')
plt.scatter(x,u[0:resolution:int(resolution/(N-1))],color='r',label='MFS')
plt.plot(X,u,color='r')
# plt.plot(x,u,color='r')
plt.legend()
plt.xlabel("x")
plt.ylabel("u(x)")
plt.show()

error=[]
for a in range(resolution):
    error.append(u[a]-u_analytical[a])
# ax4 = plt.axes
# plt.rcParams["font.family"] = "serif"
# plt.rcParams["mathtext.fontset"] = "dejavuserif"
# plt.title('Deformation error')
# # plt.plot(X,u_analytical,label='Analytical')
# plt.plot(X,error,color='r',label='error')
# # plt.plot(x,u,color='r')
# plt.legend()
# plt.xlabel("x")
# plt.ylabel("absolute error")
# plt.show()

# ax4 = plt.axes
# plt.rcParams["font.family"] = "serif"
# plt.rcParams["mathtext.fontset"] = "dejavuserif"
# plt.title('Axial bar force')
# # plt.plot(X,u_analytical,label='Analytical')
# plt.scatter(x,f,color='r',label='MFS')
# plt.plot(x,f,color='r')
# plt.legend()
# plt.xlabel("x")
# plt.ylabel("f(x)")
# plt.show()
