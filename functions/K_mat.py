'''
PyMFS stiffness matrix generation

The purpose of this script is to create an object containing the system
stiffness matrix to be passed to the PyMFS solver.
'''

'''
To-do list:

- Shape function generation
'''

import numpy as np

from functions.input_data import input_data
from .B_mat import *
from .C_mat import *
from .gauss_points import *

class K_mat:
    def __init__(self, input_data):
        self.input_data = input_data
        
        self.C = C(input_data)

        print('Generating stiffness matrix...')
        self.assemble()
        print('Stiffness matrix assembled. Time elapsed = ')
        

    def assemble(self):
        Ntot = len(self.input_data.node_coor)
        order = 3
        
        U_ImJn = np.zeros((Ntot,Ntot,order,order)) # BC stiffness matrix
        nested = np.zeros(((Ntot,Ntot,order,order))) # Nested stiffness matrix
        K_open = np.zeros((Ntot*order,Ntot*order)) # Unravelled stiffness matrix
        
        for I in range(Ntot):
            for J in range(Ntot):
                '''
                Generate gauss points and weights for current node
                and contributing node.
                '''
                I_coor = self.input_data.node_coor[I]
                J_coor = self.input_data.node_coor[J]
                gauss = gauss_points(I,J,self.input_data)
                
                '''Loop over gauss points to evaluate B matrix'''
                for m in range(order):
                    for n in range(order):
                        '''If no integration points exist in region, stiffness must be zero'''
                        if np.array(gauss.points).size==0:
                            nested[I][J][m][n] = 0
                        else:
                            '''Otherwise, evaluate the stiffness term by integrating over the
                            integration points.'''
                            for counter in range(len(gauss.points)):
                                B_Im = B(I,m,gauss.points[counter],self.input_data).mat
                                B_Jn = B(J,n,gauss.points[counter],self.input_data).mat

                                wi = np.array(gauss.weights)[counter,0]
                                wj = np.array(gauss.weights)[counter,1]

                                dx = self.input_data.radius
                                nested[I][J][m][n] += (dx/2)*(dx/2)*wj*wi*np.transpose(B_Im) @ C @ B_Jn


