'''
PyMFS force vector generation
'''

'''
To-do list:

- Maybe move boundary force terms to their own functions?
'''

import numpy as np
from functions.input_data import input_data
from .gauss_points import *
from .boundary_intpoints import *
from .B_mat import *
from .C_mat import *
from .direction_cosines import*
from shapely.geometry import LineString, Point

class f_vec():
    def __init__(self, I, system_variables):
        '''
        Read inputs
        '''
        self.I = I
        self.input_data = system_variables.input_data
        self.C = system_variables.C
        '''
        Generate f_Im
        '''
        self.vec = self.generate()
    

    def generate(self):
        '''
        Read values, generate gauss points for current region, initialise
        empty element force vector.
        '''
        I = self.I
        C = self.C
        input_data = self.input_data
        surfaces = input_data.external_surfaces
        surfaces += input_data.internal_surfaces
        
        neumann_count = 0 # Initialise counter that counts number of neumann spheres encountered
        dirichlet_count = 0 # Initialise counter that counts number of dirichlet spheres encountered

        Ndof = 3
        Ndim = 2

        gauss = gauss_points(I,I,input_data)
        # f_Im_nested = np.zeros((Ndof,Ndof,Ndim,Ndim)) # Nested version of K
        f_Im = np.zeros(Ndim)       # Un-nested version of K
        fHat_Im = np.zeros(Ndim)    # Empty boundary force vector

        '''
        Loop over surfaces with prescribed loads.
        '''
        for surface_number in input_data.loads[4]:

            surf = surfaces[surface_number[0]]     # read current surface
            surface = LineString(surf)             # create surface object
            coor = Point(input_data.node_coor[I])  # create point object from node
            dist = coor.distance(surface)          # distance from node to surface

            f_sx = input_data.loads[0][neumann_count]
            f_sy = input_data.loads[1][neumann_count]
            f_s = np.array([f_sx, f_sy])
            
            '''
            Does current sphere intersect a Neumann boundary?
            '''
            if dist < input_data.radius:            
                '''
                If yes, generate integration points along boundary,
                then loop over DoFs.
                '''
                integration_points, integration_weights, y1,y2 = boundary_intpoints(I,I,input_data)
                for m in range(Ndof):
                    '''
                    Loop over boundary integration points
                    '''
                    for counter, point in enumerate(integration_points):
                        wj = integration_weights[counter]
                        point = (np.array(surf)[0,0],integration_points[counter])
                        h_Im = shape_functions(I,m,point,input_data)
                        H_Im = np.array([[h_Im, 0],\
                                         [0, h_Im]])
                        
                        '''
                        1D Gaussian product rule
                        '''
                        fHat_Im[:] += ((y2-y1)/2)*wj*H_Im @ f_s
                        neumann_count += 1
                    

        f_Im = f_Im + fHat_Im

        return f_Im
