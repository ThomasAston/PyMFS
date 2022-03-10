'''
PyMFS stiffness matrix generation
'''

'''
To-do list:

- Shape function generation
'''

import numpy as np
from .B_mat import *
from .C_mat import *
from.apply_DoFs import *

class K_mat:
    def __init__(self, input_data):
        self.input_data = input_data
        
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
                I_coor = self.input_data.node_coor[I]
                J_coor = self.input_data.node_coor[J]
                gauss = gauss_points(I,J,self.input_data)
        
