# Let's generate the K matrix, for the entire domain?

# To do list:
# Include the calculation of stiffness of overlap regions...
#
#

import numpy as np
from .polynomial_basis import *
from .B_mat import *
from .C_mat import *
from.apply_DoFs import *

class K_mat:
    def __init__(self, domain, nodes, elements, integration_points, DoFs):
        # print('Generating stiffness matrix...')
        self.domain = domain
        self.nodes = nodes
        self.elements = elements
        self.polynomial_basis = polynomial_basis(order=0)  
        self.integration_points = integration_points      
        self.DoFs = DoFs

        self.value = self.assemble_K()

    # Function for assembling the global stiffness matrix
    def assemble_K(self):

        K = np.zeros([2*self.nodes.num, 2*self.nodes.num])
        # K = []    
        
        # Generate the BC flag, used to determine if node is fixed.
        DoFs = self.DoFs.order
        DoFs = DoFs.astype(int)

        for I in range(0,2*self.nodes.num,2):
                for J in range(0, 2*self.nodes.num,2):
                    self.I = I
                    self.J = J

                    K[DoFs[I],DoFs[J]] = self.element_k()[0,0]
                    K[DoFs[I],DoFs[J+1]] = self.element_k()[0,1]
                    K[DoFs[I+1],DoFs[J]] = self.element_k()[1,0]
                    K[DoFs[I+1],DoFs[J+1]] = self.element_k()[1,1]

        return K
    
    # Function for calculating the stiffness matrix for each global coordinate 
    def element_k(self):
        I = self.I
        J = self.J
        
        # Work out the current node we're looking at based on I and J 
        x_node = self.nodes.coor[int(I/2)][0]
        y_node = self.nodes.coor[int(I/2)][1]
        pos_node = np.array((x_node,y_node))

        t = 0.001
        # Generate the D matrix
        my_material = material(E=70e9, v=0)
        my_Dmat = D_mat(my_material)
        D = my_Dmat.D

        k = np.zeros([2*len(self.polynomial_basis.value),2*len(self.polynomial_basis.value)])
        x_J = self.nodes.coor[int(J/2)][0]
        y_J = self.nodes.coor[int(J/2)][1]
        pos_J = np.array((x_J,y_J))
        if I == J:
            k = np.zeros([2*len(self.polynomial_basis.value),2*len(self.polynomial_basis.value)])
            for gp in range(0,len(self.integration_points.gauss_points)):
                
                # print(self.integration_points.gauss_points[gp,:])
                x_gauss = self.integration_points.gauss_points[gp,0]
                y_gauss = self.integration_points.gauss_points[gp,1]
                pos_gauss = np.array((x_gauss,y_gauss))

                w_i = self.integration_points.gauss_weights[gp,0]
                w_j = self.integration_points.gauss_weights[gp,1]
                # If distance between gauss point and node is inside sphere radius...
                if np.linalg.norm(pos_gauss-pos_node) < self.elements.size:
                    B_Im = B_mat(self.domain, self.nodes, self.elements, x_gauss, y_gauss, int(I/2))

                    for J in range(0, self.nodes.num):
                        B_Jn = B_mat(self.domain, self.nodes, self.elements, x_gauss, y_gauss, int(J))
              
                        k = k + w_j * w_i * t * np.transpose(B_Im.value) @ D @ B_Jn.value
                else: 
                    k = k
        else: #abs(np.linalg.norm(pos_J-pos_node)) <= 2*self.elements.size:
            for gp in range(0,len(self.integration_points.gauss_points)):
                x_gauss = self.integration_points.gauss_points[gp,0]
                y_gauss = self.integration_points.gauss_points[gp,1]
                pos_gauss = np.array((x_gauss,y_gauss))

                w_i = self.integration_points.gauss_weights[gp,0]
                w_j = self.integration_points.gauss_weights[gp,1]

                # If gauss point lies in an overlap region
                if np.linalg.norm(pos_gauss-pos_node) < self.elements.size and np.linalg.norm(pos_gauss-pos_J) < self.elements.size:
                    B_Im = B_mat(self.domain, self.nodes, self.elements, x_gauss, y_gauss, int(I/2))

                    for J in range(0, self.nodes.num):
                        B_Jn = B_mat(self.domain, self.nodes, self.elements, x_gauss, y_gauss, int(J))
              
                        k = k + w_j * w_i * t * np.transpose(B_Im.value) @ D @ B_Jn.value
                else: 
                    k = k      

        return k

