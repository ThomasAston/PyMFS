# meshLess_2Day v1
# Author(s): Thomas Aston
# Date updated: 02/12/2021 

# This is the central file for solving problems by calling the various functions required by the meshLess solver.

# To do list:
# Post-processing module
# Examine the choice of basis function
#

from functions import *
import time


start = time.time()
##########################################################################
# Set up the domain, select material
##########################################################################
my_domain = rectangle(point1=[0,0], point2=[2,2])
# my_material = material(E=200e9, v=0.3)
# my_Dmat = D_mat(my_material)

##########################################################################
# Discretise it by generating nodes
##########################################################################
my_nodes = nodes(my_domain, N=9, method='Regular')


##########################################################################
# Generate sphere elements and shape functions 
##########################################################################
my_elements = elements(my_nodes, my_domain, type='MFS', size=my_nodes.coor[1,0]-my_nodes.coor[0,0])


##########################################################################
# Impose boundary conditions on the nodes, apply prescribed forces
##########################################################################
my_DoFs = apply_DoFs(my_nodes, BC = 'bottom fixed', load = [0, 10000], location = [1,2])

##########################################################################
# Generate integration points 
##########################################################################
my_intPoints = integration_points(my_domain, my_nodes, my_elements, degree=6)


##########################################################################
# Integrate to find B and K
##########################################################################
my_basis = polynomial_basis(order=0)

my_Bmat = B_mat(my_domain, my_nodes, my_elements, x=1, y=0, current_node = 4)
my_Bmat.example_plot()

my_Kmat = K_mat(my_domain, my_nodes, my_elements, my_intPoints, my_DoFs)

##########################################################################
# Solve the system of equations 
##########################################################################
my_solution = solve(my_Kmat, my_DoFs)
print(my_solution.u)
print(my_solution.F)

end = time.time()
print('Time taken to solve: ',end-start)
##########################################################################
# Post processing 
##########################################################################
my_post = post(my_solution, my_domain, my_nodes, my_elements, my_intPoints)

