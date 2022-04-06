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

        Ndof = 4
        Ndim = 2

        gauss = gauss_points(I,I,input_data)
        
        f_Im = np.zeros([Ndof,Ndim])       # Un-nested version of K
        fHat_Im = np.zeros([Ndof,Ndim])    # Empty boundary force vector

        '''
        Loop over surfaces with prescribed loads.
        '''
        for count,surface_number in enumerate(input_data.loads[4]):

            surf = surfaces[surface_number[0]]     # read current surface
            surface = LineString(surf)             # create surface object
            coor = Point(input_data.node_coor[I])  # create point object from node
            dist = coor.distance(surface)          # distance from node to surface

            f_sx = input_data.loads[0][count]
            f_sy = input_data.loads[1][count]
            f_s = np.array([f_sx, f_sy])
            
            '''
            Does current sphere intersect a Neumann boundary?
            '''
            if dist < input_data.radius:            
                '''
                If yes, loop over DoFs for current node,  
                then move along the surface.
                '''
                # integration_points, integration_weights, t1,t2 = boundary_intpoints(I,I,surf,input_data)
                integration_points,integration_weights, t1, t2, type = boundary_intpoints(I,I,surf,input_data)
                for m in range(Ndof):
                    '''
                    Loop over boundary integration points
                    '''
                    for counter, point in enumerate(integration_points):
                        wj = integration_weights[counter]
                        
                        # point = (np.array(surf)[0,0],integration_points[counter])
                        
                        h_Im = shape_functions(I,m,point,input_data)
                        H_Im = np.array([[h_Im, 0],\
                                         [0, h_Im]])

                        '''
                        1D Gaussian product rule
                        '''
                        fHat_Im[m,:] += ((t2-t1)/2)*wj*H_Im @ f_s
                        # integrand = H_Im @ f_s
                        # fHat_Im[m,:] = np.trapz(integrand,x=t_int)
            # neumann_count += 1

        '''
        Loop over surfaces with prescribed displacements
        '''
        for surface_number in input_data.BC[4]:
            surf = surfaces[surface_number[0]]     # read current surface
            surface = LineString(surf)             # create surface object
            coor = Point(input_data.node_coor[I])  # create point object from node
            dist = coor.distance(surface)          # distance from node to surface

            u_s = input_data.BC[0][dirichlet_count]
            v_s = input_data.BC[1][dirichlet_count]
            u_s = np.array([u_s, v_s])
            
            '''
            Does current sphere intersect a Dirichlet boundary?
            '''
            if dist < input_data.radius:            
                '''
                If yes, loop over DoFs for current node,  
                then move along the surface.
                '''
                # integration_points, integration_weights, t1,t2 = boundary_intpoints(I,I,surf,input_data)
                integration_points,integration_weights, t1, t2, type = boundary_intpoints(I,I,surf,input_data)
                for m in range(Ndof):
                    '''
                    Loop over boundary integration points
                    '''
                    for counter, point in enumerate(integration_points):
                        wj = integration_weights[counter]
                        # point = (np.array(surf)[0,0],integration_points[counter])
                        B_Im = B(I,m,point,input_data).mat
                        N = direction_cosines(surf)
                        '''
                        1D Gaussian product rule
                        '''
                        fHat_Im[m,:] -= ((t2-t1)/2)*wj*np.transpose(B_Im) @ self.C @ np.transpose(N) @ u_s
                        # integrand = np.transpose(B_Im) @ self.C @ np.transpose(N) @ u_s
                        # fHat_Im[m,:] = -np.trapz(integrand,x=t_int)
                        # dirichlet_count += 1

        f_Im = f_Im + fHat_Im

        return f_Im
