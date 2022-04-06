# meshLess_2D v1

# Author(s): Thomas Aston
# Date updated: 09/02/2022

# Main file for solving 2D Poisson problems using MFS

# To do list:
# 
# 
# 
import math
from random import gauss
import numpy as np
import sympy as sym
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon

###############################################
# PROBLEM DEFINITION
##############################################

# We want to solve Poisson's equation on a 2D bi-unit square. The analytical
# displacement field is chosen such that tht following domain properties and
# values are valid. 

### Domain ###

Ndim = 2
order = 4
lx = 2  # Length in x-direction
ly = 2  # Length in y-direction

vertices = np.array([[-lx/2,-lx/2],[-lx/2,lx/2],[lx/2,lx/2],[lx/2,-lx/2]])
poly = Polygon(vertices)                                            

### Discretisation ###
Nx = 3  # number of nodes in x-direction 
Ny = 3 # number of nodes in y-direction
Ntot = Nx*Ny # total number of nodes

x = np.linspace(-lx/2, lx/2, Nx) # initialise x coords
y = np.linspace(-ly/2, ly/2, Ny) # initialise y coords

dx = lx/(Nx-1) # node spacing in x-direction
dy = ly/(Ny-1) # node spacing in y-direction

r = dx  # Sphere radius

# Generate coordinates
xv, yv = np.meshgrid(x, y)        
coor = np.array(list(zip(xv.ravel(), yv.ravel())))
# print('Nodal coordinates: ', coor)

### Analytical displacement field ###
NX = 30
NY = 30
X = np.linspace(-lx/2,lx/2,NX)
Y = np.linspace(-ly/2,ly/2,NY)

# Assemble analytical coordinates
coor_analytical = []
for i in X:
    for j in Y:
        coor_analytical.append([i,j])


def u_analytical(x,y):
    return (7*x+x**7) * np.cos(np.pi*y)


##############################################
# GENERATE SHAPE FUNCTIONS
##############################################

# W_i = np.zeros((Nx, Ny, NX, NY))

def shape_functions(node,point):
    xi = coor[node,0]
    yi = coor[node,1]

    point_x = point[0]
    point_y = point[1]

    s = math.hypot(point_x-xi, point_y-yi)/r

    if s<=1:
        W_i = 1-6*s**2+8*s**3-3*s**4
    else:
        W_i = 0

    sum_W_j = 0
    for j in range(Ntot):
        xj = coor[j,0]
        yj = coor[j,1]
        s_j = math.hypot(point_x-xj, point_y-yj)/r
        if s_j<=1:
            W_j=1-6*s_j**2+8*s_j**3-3*s_j**4
            sum_W_j+=W_j
        else:
            continue
    
    phi_i = W_i/sum_W_j 

    h = np.zeros([order,1])
    h[0] = phi_i
    h[1] = phi_i*(point_x-xi)/r
    h[2] = phi_i*(point_y-yi)/r
    h[3] = phi_i*(point_x-xi)/r*(point_y-yi)/r

    return W_i, h
# print(phi)

# def shape_functions(x,y):


##############################################
# GENERATE INTEGRATION POINTS
##############################################

gauss_grid = np.array(list(zip(xv.ravel(), yv.ravel())))
gauss_grid = Ntot*[gauss_grid]

def integration_points(I,J,int_degree):
    node_x = coor[I,0]
    node_y = coor[I,1]

    contrib_x = coor[J,0]
    contrib_y = coor[J,1]

    x_int, w_int = np.polynomial.legendre.leggauss(int_degree)
    xi, xj = np.meshgrid(x_int, x_int)
    wi, wj = np.meshgrid(w_int, w_int)


    gauss_grid = np.array(list(zip(xi.ravel(), xj.ravel())))
    weight_grid = np.array(list(zip(wi.ravel(), wj.ravel())))
    gauss_points = []
    gauss_weights = []

    if I==J: # If I = J we are integrating an interior/boundary sphere directly
        # Quadrant 1
        for a in range(int_degree**2):
            x = node_x - 0.5*r + 0.5*r*gauss_grid[a,0]
            y = node_y + 0.5*r + 0.5*r*gauss_grid[a,1]
            point = (x,y)
            pointInDomain = poly.contains_point(point)
            if pointInDomain == True and math.hypot(x-node_x, y-node_y) < r:
                gauss_points.append([x,y])
                gauss_weights.append([weight_grid[a,0], weight_grid[a,1]])
            else:   
                continue   
        # Quadrant 2
        for a in range(int_degree**2):
            x = node_x - 0.5*r + 0.5*r*gauss_grid[a,0]
            y = node_y - 0.5*r + 0.5*r*gauss_grid[a,1]
            point = (x,y)
            pointInDomain = poly.contains_point(point)
            if pointInDomain == True and math.hypot(x-node_x, y-node_y) < r:
                gauss_points.append([x,y])
                gauss_weights.append([weight_grid[a,0], weight_grid[a,1]])
            else:
                continue   
        # Quadrant 3
        for a in range(int_degree**2):
            x = node_x + 0.5*r + 0.5*r*gauss_grid[a,0]
            y = node_y - 0.5*r + 0.5*r*gauss_grid[a,1]
            point = (x,y)
            pointInDomain = poly.contains_point(point)
            if pointInDomain == True and math.hypot(x-node_x, y-node_y) < r:
                gauss_points.append([x,y])
                gauss_weights.append([weight_grid[a,0], weight_grid[a,1]])
            else:
                continue    
        # Quadrant 4
        for a in range(int_degree**2):
            x = node_x + 0.5*r + 0.5*r*gauss_grid[a,0]
            y = node_y + 0.5*r + 0.5*r*gauss_grid[a,1]
            point = (x,y)
            pointInDomain = poly.contains_point(point)
            if pointInDomain == True and math.hypot(x-node_x, y-node_y) < r:
                gauss_points.append([x,y])
                gauss_weights.append([weight_grid[a,0], weight_grid[a,1]])
            else:
                continue      
    else:  # Otherwise we are integrating a sphere overlap region
        if math.hypot(node_x-contrib_x, node_y-contrib_y) < 2*r:
        # Quadrant 1
            for a in range(int_degree**2):
                x = node_x - 0.5*r + 0.5*r*gauss_grid[a,0]
                y = node_y + 0.5*r + 0.5*r*gauss_grid[a,1]
                point = (x,y)
                pointInDomain = poly.contains_point(point)
                if pointInDomain == True and math.hypot(x-node_x, y-node_y) < r and math.hypot(x-contrib_x, y-contrib_y) < r:
                    gauss_points.append([x,y])
                    gauss_weights.append([weight_grid[a,0], weight_grid[a,1]])
                else:
                    continue   
            # Quadrant 2
            for a in range(int_degree**2):
                x = node_x - 0.5*r + 0.5*r*gauss_grid[a,0]
                y = node_y - 0.5*r + 0.5*r*gauss_grid[a,1]
                point = (x,y)
                pointInDomain = poly.contains_point(point)
                if pointInDomain == True and math.hypot(x-node_x, y-node_y) < r and math.hypot(x-contrib_x, y-contrib_y) < r:
                    gauss_points.append([x,y])
                    gauss_weights.append([weight_grid[a,0], weight_grid[a,1]])
                else:
                    continue   
            # Quadrant 3
            for a in range(int_degree**2):
                x = node_x + 0.5*r + 0.5*r*gauss_grid[a,0]
                y = node_y - 0.5*r + 0.5*r*gauss_grid[a,1]
                point = (x,y)
                pointInDomain = poly.contains_point(point)
                if pointInDomain == True and math.hypot(x-node_x, y-node_y) < r and math.hypot(x-contrib_x, y-contrib_y) < r:
                    gauss_points.append([x,y])
                    gauss_weights.append([weight_grid[a,0], weight_grid[a,1]])
                else:
                    continue    
            # Quadrant 4
            for a in range(int_degree**2):
                x = node_x + 0.5*r + 0.5*r*gauss_grid[a,0]
                y = node_y + 0.5*r + 0.5*r*gauss_grid[a,1]
                point = (x,y)
                pointInDomain = poly.contains_point(point)
                if pointInDomain == True and math.hypot(x-node_x, y-node_y) < r and math.hypot(x-contrib_x, y-contrib_y) < r:
                    gauss_points.append([x,y])
                    gauss_weights.append([weight_grid[a,0], weight_grid[a,1]])
                else:
                    continue    
    
    # Determine integral limits
    x1 = max(-lx/2,node_x-r, contrib_x-r)
    x2 = min(lx/2,node_x+r,contrib_x+r)
    y1 = max(-lx/2,node_y-r,contrib_y-r)
    y2 = min(lx/2,node_y+r,contrib_y+r)

    gauss_points = np.array(gauss_points)
    gauss_weights = np.array(gauss_weights)
    
    return gauss_points, gauss_weights, x1,x2,y1,y2

##############################################
# ASSEMBLE STIFFNESS MATRIX
##############################################
### Initialise empty arrays ###
# K_nested = np.zeros((Ntot,Ntot,order,order)) # Stiffness matrix
# KU_ImJn = np.zeros((Ntot,Ntot,order,order))

# K_open = np.zeros((Ntot*order,Ntot*order)) # Unravelled stiffness matrix

# for I in range(Ntot):
#     for J in range(Ntot):
#         gauss_points, gauss_weights,x1,x2,y1,y2 = integration_points(I,J,int_degree=6)
        
#         for m in range(order):
#             for n in range(order):
#                 if gauss_points.size==0:
#                     K_nested[I][J][m][n] = 0
#                 else:
#                     for counter in range(len(gauss_points)):
#                         # Generate shape functions for each point
#                         point = gauss_points[counter]
#                         h_I = shape_functions(I,point)
#                         h_J = shape_functions(J,point)
#                         delta = 0.0001
#                         xplus = point[0]+delta
#                         yplus = point[1]+delta
#                         point_xplus = (xplus,point[1])
#                         point_yplus = (point[0],yplus)
                        
#                         hI_xplus = shape_functions(I,point_xplus)
#                         hI_yplus = shape_functions(I,point_yplus)
#                         hJ_xplus = shape_functions(J,point_xplus)
#                         hJ_yplus = shape_functions(J,point_yplus)

#                         dhI_dx = (hI_xplus-h_I)/delta
#                         dhI_dy = (hI_yplus-h_I)/delta
#                         dhJ_dx = (hJ_xplus-h_J)/delta
#                         dhJ_dy = (hJ_yplus-h_J)/delta

#                         wi = gauss_weights[counter,0]
#                         wj = gauss_weights[counter,1]

#                         K_nested[I][J][m][n] += (dx/2)*(dy/2)*wj*wi*(dhI_dx[m]*dhJ_dx[n]+dhI_dy[m]*dhJ_dy[n])
#                         # K_nested[I][J][m][n] += ((x2-x1)/2)*((y2-y1)/2)*wj*wi*(dhI_dx[m]*dhJ_dx[n]+dhI_dy[m]*dhJ_dy[n])
#                         # K_nested[I][J][m][n] += wj*wi*(dhI_dx[m]*dhJ_dx[n]+dhI_dy[m]*dhJ_dy[n])
                        
#                         if coor[I,0] + r > lx/2: # sphere intersects dirichlet boundary
#                             dhImhJn_dx = (hI_xplus[m]*hJ_xplus[n] - h_I[m]*h_J[n])/delta
#                             # K_nested[I][J][m][n] -= ((x2-x1)/2)*((y2-y1)/2)*wj*wi*dhImhJn_dx
#                             K_nested[I][J][m][n] -= (dx/2)*(dy/2)*wj*wi*dhImhJn_dx

#                         K_open[order*I:order+order*I,order*J:order+order*J] = K_nested[I,J,:,:]
                    

##############################################
# IMPLEMENT BOUNDARY CONDITIONS
##############################################

# f_Im = np.zeros((Ntot, order))
# fHat_Im = np.zeros((Ntot, order))
# # KU_ImJn = np.zeros((Ntot,Ntot,order,order))
# # KU_ImJn_open = np.zeros((Ntot*order,Ntot*order))

# for I in range(Ntot):
#     for m in range(order):
       
#         gauss_points, gauss_weights,x1,x2,y1,y2 = integration_points(I,I,int_degree=6)
        
#         for counter in range(len(gauss_points)):
#             point = gauss_points[counter]
#             h_I = shape_functions(I,point)
#             wi = gauss_weights[counter,0]
#             wj = gauss_weights[counter,1]
            
#             f_xy = ((np.pi**2)*(7*point[0]+point[0]**7)-42*point[0]**5)*np.cos(np.pi*point[1])
#             f_Im[I][m] +=  (dx/2)*(dy/2)*wj*wi*h_I[m]*f_xy
#             # f_Im[I][m] +=  wj*wi*h_I[m]*f_xy
        
#         if coor[I,0] - r < -lx/2: # sphere intersects Neumann boundary
#             for counter in range(len(gauss_points)):
#                 point = gauss_points[counter]
#                 h_I = shape_functions(I,point)
#                 wi = gauss_weights[counter,0]
#                 wj = gauss_weights[counter,1]

#                 f_s = -14*np.cos(np.pi*point[1])
#                 # fHat_Im[I,m] += ((x2-x1)/2)*((y2-y1)/2)*wj*wi*h_I[m]*f_s
                
#                 fHat_Im[I,m] += (dx/2)*(dy/2)*wj*wi*h_I[m]*f_s
        
#         elif coor[I,0] + r > lx/2: # sphere intersects Dirichlet boundary
#             for counter in range(len(gauss_points)):
#                 point = gauss_points[counter]
#                 h_I = shape_functions(I,point)
#                 wi = gauss_weights[counter,0]
#                 wj = gauss_weights[counter,1]
                
#                 delta = 0.0001
#                 xplus = point[0]+delta
#                 yplus = point[1]+delta
#                 point_xplus = (xplus,point[1])
                
#                 hI_xplus = shape_functions(I,point_xplus)
                
#                 dhI_dx = (hI_xplus-h_I)/delta
#                 u_s = 8*np.cos(np.pi*point[1])
#                 # fHat_Im[I][m] -= ((x2-x1)/2)*((y2-y1)/2)*wj*wi*u_s*dhI_dx[m]
#                 fHat_Im[I][m] -= (dx/2)*(dy/2)*wj*wi*u_s*dhI_dx[m]

            # for J in range(Ntot):
            #     # if I==J:
            #     for n in range(order):
            #         for counter in range(len(gauss_points)):
            #             point = gauss_points[counter]
            #             h_I = shape_functions(I,point)
            #             wi = gauss_weights[counter,0]
            #             wj = gauss_weights[counter,1]
                        
            #             delta = 0.001
            #             xplus = point[0]+delta
            #             yplus = point[1]+delta
            #             point_xplus = (xplus,point[1])

            #             hI_xplus = shape_functions(I,point_xplus)
            #             dhI_dx = (hI_xplus-h_I)/delta
            #             h_J = shape_functions(J,point)
            #             hJ_xplus = shape_functions(J,point_xplus)
                        
            #             dhImhJn_dx = (hI_xplus[m]*hJ_xplus[n] - h_I[m]*h_J[n])/delta
            
            #             KU_ImJn[I][J][m][n] += ((x2-x1)/2)*((y2-y1)/2)*wj*wi*dhImhJn_dx
            #             KU_ImJn_open[order*I:order+order*I,order*J:order+order*J] = KU_ImJn[I,J,:,:]
        

##############################################
# SOLVE SYSTEM OF EQUATIONS
##############################################
# K_open = K_open-KU_ImJn_open
# print(f_Im)

# f = f_Im + fHat_Im
# f = np.ravel(f)

# q = np.linalg.solve(K_open,f)

##############################################
# REASSEMBLE FUNCTION FIELD
##############################################
# q = np.reshape(q, [Ntot,order])

# print(q)
# # u = np.zeros([Ntot,1])
# u=np.zeros([NX,NY])

# plotX, plotY = np.meshgrid(X,Y)
# for a in range(NX):
#     for b in range(NY):
#         point_x = X[a]
#         point_y = Y[b]
#         point = (point_x,point_y)
    
#         for J in range(Ntot):
#             h_J = shape_functions(J,point)
#             for n in range(order):
#                 u[b,a] += h_J[n]*q[J,n]

# for node in range(len(coor)):
#     point = coor[node]
#     for J in range(Ntot):
#         h_J = shape_functions(J,point)
#         for n in range(order):
#             u[node] += h_J[n]*q[J,n]

##############################################
# PLOT RESULTS
##############################################
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{cmbright}')


gauss_plot, weights,x1,x2,y1,y2 = integration_points(2,5,int_degree=6)
# print(weights)
ax10 = plt.axes()
plt.scatter(gauss_plot[:,0], gauss_plot[:,1])


# for i in range(len(gauss_plot[:,0])):
#     ax10.annotate(i, (gauss_plot[i,0], gauss_plot[i,1]))

plt.show()


# ax = plt.axes(projection='3d')
# plotX, plotY = np.meshgrid(X,Y)
# ax.contour(plotX,plotY,u,50,cmap='viridis')
# ax.scatter(coor[:,0],coor[:,1],s=10,color='k')
# ax.zaxis.set_rotate_label(False)
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('u(x,y', rotation=0)
# ax.set_title('MFS field')
# # ax.set_xlabel(r'$x$')
# # ax.set_ylabel(r'$y$')
# # ax.set_zlabel(r'$u(x,y)$', rotation = 0);
# plt.grid(False)
# plt.show()


# ax2 = plt.axes(projection='3d')
# plotX, plotY = np.meshgrid(X,Y)
# u_analytical = u_analytical(plotX,plotY)

# ax2.contour(plotX,plotY,u_analytical,50,cmap='viridis')
# ax2.zaxis.set_rotate_label(False)
# ax2.set_xlabel('x')
# ax2.set_ylabel('y')
# ax2.set_zlabel('u(x,y', rotation=0)
# ax2.set_title('Analytical field')
# # ax.set_xlabel(r'$x$')
# # ax.set_ylabel(r'$y$')
# # ax.set_zlabel(r'$u(x,y)$', rotation = 0);
# plt.grid(False)
# plt.show()

# u_error = u - u_analytical
# ax3 = plt.axes(projection='3d')
# ax3.contour(plotX,plotY,u_error,50,cmap='viridis')
# ax3.zaxis.set_rotate_label(False)
# ax3.set_xlabel('x')
# ax3.set_ylabel('y')
# ax3.set_zlabel('Error', rotation=0)
# ax3.set_title('Absolute error in MFS solution')
# # ax.set_xlabel(r'$x$')
# # ax.set_ylabel(r'$y$')
# # ax.set_zlabel(r'$u(x,y)$', rotation = 0);
# plt.grid(False)
# plt.show()

h_plot = np.zeros([NX,NY])
W_plot = np.zeros([NX,NY])

W_save = []
phi_save = []
phi1_save = []
phi2_save = []
phi3_save = []
coordinates_save = []
for a in range(NX):
    for b in range(NY):

        point = (X[a],Y[b])
        point_xplus = (point[0],point[1]+0.01)

        W, h = shape_functions(4,point)

        W_save.append((a,b,W))
        phi_save.append((a,b,h[0]))
        phi1_save.append((a,b,h[1]))
        phi2_save.append((a,b,h[2]))
        phi3_save.append((a,b,h[3]))
        
        W_plot[b,a] = W
        h_plot[b,a] = h[3]
        # h_plot[b,a] = (shape_functions(4,point_xplus)[0]-shape_functions(4,point)[0])/0.01

np.savetxt("W_plot.txt", W_save)
np.savetxt("phi_plot.txt", phi_save)
np.savetxt("phi1_plot.txt", phi1_save)
np.savetxt("phi2_plot.txt", phi2_save)
np.savetxt("phi3_plot.txt", phi3_save)

ax = plt.axes(projection='3d')
plotX, plotY = np.meshgrid(X,Y)
ax.contour(plotX,plotY,h_plot,50,cmap='viridis')
# ax.scatter(coor[:,0],coor[:,1],f_Im[:,0],s=10,color='k')
# ax.zaxis.set_rotate_label(False)
ax.set_xlabel('x')
ax.set_ylabel('y')
# ax.set_zlabel('u(x,y', rotation=0)
ax.set_title('Shape functions')
# ax.set_xlabel(r'$x$')
# ax.set_ylabel(r'$y$')
# ax.set_zlabel(r'$u(x,y)$', rotation = 0);
plt.grid(False)
plt.show()
