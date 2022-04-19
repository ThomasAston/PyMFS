'''
PyMFS stiffness matrix generation
'''

'''
To-do list:

- Add boundary stiffness term... maybe in form of another function?
'''

import numpy as np
from functions.input_data import input_data
from .gauss_points import *
from .B_mat import *
from .C_mat import *
from .direction_cosines import*
from .boundary_intpoints import*
from shapely.geometry import LineString, Point

class K_mat:
    def __init__(self, I, J, system_variables):
        '''
        Read inputs
        '''
        self.I = I
        self.J = J
        self.input_data = system_variables.input_data
        self.C = system_variables.C
        '''
        Generate K_ImJn
        '''
        self.mat = self.generate()


    def generate(self):
        '''
        Read values, generate gauss points for current region, initialise
        empty element stiffness matrices.
        '''
        I = self.I
        J = self.J
        C = self.C
        input_data = self.input_data
        surfaces = input_data.external_surfaces[:]
        surfaces += input_data.internal_surfaces[:]
        Ndof = 4
        Ndim = 2
        gauss = gauss_points(I,J,input_data)
        K_ImJn_nested = np.zeros((Ndof,Ndof,Ndim,Ndim)) # Nested version of K
        K_ImJn = np.zeros((Ndof*Ndim, Ndof*Ndim)) # Un-nested version of K
        KU_ImJn_nested = np.zeros((Ndof,Ndof,Ndim,Ndim)) # Nested version of DIRICHLET K
        KU_ImJn = np.zeros((Ndof*Ndim, Ndof*Ndim)) # Un-nested version of DIRICHLET K
        dirichlet_count = 0

        '''
        Loop over DoFs
        '''
        for m in range(Ndof):
            for n in range(Ndof):
                '''
                If no integration points exist in the current region, K_ImJn will
                be full of zeros and so we can skip calculations...
                '''
                if np.array(gauss.points).size==0:
                    type = ''
                    continue
                else:
                    '''
                    Otherwise, we must... ASSEMBLE STIFFNESS MATRIX.
                    Evaluate stiffness terms by integrating over the
                    integration points.
                    '''
                    for counter, point in enumerate(gauss.points):                        
                        '''
                        Evaluate B_Im and B_Jn and gauss coefficients.
                        '''
                        B_Im = B(I,m,point,input_data).mat
                        B_Jn = B(J,n,point,input_data).mat

                        wi = np.array(gauss.weights)[counter,0]
                        wj = np.array(gauss.weights)[counter,1]
                        dx = input_data.radius

                        K_ImJn_nested[m,n] += (dx/2)*(dx/2)*wj*wi*np.transpose(B_Im) @ C @ B_Jn
                        K_ImJn[m*2:m*2+2, n*2:n*2+2] = K_ImJn_nested[m,n]
                    '''
                    Loop over surfaces with prescribed displacements
                    '''
                    for surface_number in input_data.BC[4]:
                        surf = surfaces[surface_number[0]]     # read current surface
                        surface = LineString(surf)             # create surface object
                        coorI = Point(input_data.node_coor[I])  # create point object from node
                        distI = coorI.distance(surface)          # distance from node to surface
                        coorJ = Point(input_data.node_coor[J])  # create point object from node
                        distJ = coorJ.distance(surface)          # distance from node to surface

                        u_s = input_data.loads[0][dirichlet_count]
                        v_s = input_data.loads[1][dirichlet_count]
                        u_s = np.array([u_s, v_s])
                        
                        '''
                        Is current node on Dirichlet boundary?
                        '''
                        if distI < input_data.radius and distJ < input_data.radius:            
                            '''
                            If yes, loop over DoFs for current node,  
                            then move along the surface.
                            '''

                            '''
                            Loop over boundary integration points
                            '''
                            integration_points,integration_weights, t1, t2, type = boundary_intpoints(I,I,surf,input_data)
                            for counter, point in enumerate(integration_points):
                                wj = integration_weights[counter]

                                B_Im = B(I,m,point,input_data).mat
                                B_Jn = B(J,n,point,input_data).mat
                                h_Im = shape_functions(I,m,point,input_data)
                                H_Im = np.array([[h_Im, 0],\
                                                 [0, h_Im]])
                                h_Jn = shape_functions(J,n,point,input_data)
                                H_Jn = np.array([[h_Jn, 0],\
                                                 [0, h_Jn]])                                        
                                N = direction_cosines(surf)
                                '''
                                1D Gaussian product rule
                                ''' 
                                KU_ImJn_nested[m,n] += ((t2-t1)/2)*wj*(H_Im @ N @ self.C @ B_Jn + np.transpose(B_Im) @ self.C @ np.transpose(N) @ H_Jn)
                                KU_ImJn[m*2:m*2+2, n*2:n*2+2] = KU_ImJn_nested[m,n]
                        
                            type = 'Dirichlet'
                        else:
                            type = 'Interior'
        

        return K_ImJn, KU_ImJn, type