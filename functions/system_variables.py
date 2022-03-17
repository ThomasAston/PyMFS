'''
PyMFS generation of variables required for solving algebraic system
'''

'''
To-do list:

- Maybe combine K and f loops to make faster...
'''

import numpy as np
from functions.input_data import input_data
from .C_mat import *
from .K_mat import *
from .f_vec import *
from .gauss_points import *


class system_variables:
    def __init__(self, input_data):
        self.input_data = input_data
        self.C = C(input_data)

        print('Generating algebraic system...')
        self.K, self.f = self.assemble()
        print('Algebraic system assembled.')
        
    '''
    Function which loops over every sphere region to generate
    global stiffness matrix and force vector.
    '''
    def assemble(self):
        '''
        Read values and initialise empty stiffness matrix and force
        vector.
        '''
        Ntot = len(self.input_data.node_coor) # Number of nodes
        Ndof = 3 # Number of DoFs per direction of each node
        Ndim = 2 # Number of physical dimensions of problem domain.
        
        K = np.zeros((Ntot*Ndof*Ndim,Ntot*Ndof*Ndim))
        f = np.zeros((Ntot*Ndof*Ndim))
    
        '''
        Loop over nodes I and J (direct spheres and overlap regions)
        and evaluate stiffness matrix and force vector components.

        Force components are only evaluated for direct spheres whereas
        stiffness components include direct spheres AND overlap regions.
        '''
        for I in range(Ntot):
            f_Im = f_vec(I,self).vec
            f[Ndim*Ndof*I:Ndim*Ndof+Ndim*Ndof*I] = f_Im
            for J in range(Ntot):   
                K_ImJn = K_mat(I,J,self).mat
                K[Ndim*Ndof*I:Ndim*Ndof+Ndim*Ndof*I,Ndim*Ndof*J:Ndim*Ndof+Ndim*Ndof*J] = K_ImJn

        '''
        Loop over nodes I (direct spheres only) and evaluate force 
        vector components.   
        '''
        
        
        # for I in range(Ntot):
            
        
        
        #         '''
        #         Generate gauss points and weights for current node
        #         and contributing node.
        #         '''
        #         I_coor = self.input_data.node_coor[I]
        #         J_coor = self.input_data.node_coor[J]
        #         gauss = gauss_points(I,J,self.input_data)
                
        #         '''
        #         Generate list of all surfaces, to be indexed to determine which
        #         surfaces BCs are applied to.
        #         '''
        #         surfaces = self.input_data.external_surfaces
        #         surfaces += self.input_data.internal_surfaces
        #         # dirichlet_count=0
        #         # neumann_count=0
        #         '''Loop over gauss points to evaluate B matrix'''
        #         for m in range(order):
        #             for n in range(order):
        #                 '''
        #                 If no integration points exist in region, stiffness must be zero
        #                 '''
        #                 if np.array(gauss.points).size==0:
        #                     K_nested[I][J][m][n] = 0
        #                 else:
        #                     '''
        #                     STIFFNESS MATRIX ASSEMBLY
        #                     Evaluate the stiffness term by integrating over the
        #                     integration points.
        #                     '''
        #                     '''
        #                     FORCE VECTOR ASSEMBLY
        #                     If I=J and m=n, we can also evaluate f_Im inside this loop!
        #                     '''
        #                     for counter in range(len(gauss.points)):
        #                         B_Im = B(I,m,gauss.points[counter],self.input_data).mat
        #                         B_Jn = B(J,n,gauss.points[counter],self.input_data).mat

        #                         wi = np.array(gauss.weights)[counter,0]
        #                         wj = np.array(gauss.weights)[counter,1]

        #                         dx = self.input_data.radius
        #                         K_nested[I][J][m][n][:][:] += (dx/2)*(dx/2)*wj*wi*np.transpose(B_Im) @ self.C @ B_Jn
        #                         K_almostOpen[I,J,2*m:2+2*m,2*n:2+2*n] = K_nested[I,J,m,n]
        #                         K_open[2*order*I:2*order+2*order*I,2*order*J:2*order+2*order*J] = K_almostOpen[I,J]

        #                         dirichlet_count=0
        #                         for surface_number in self.input_data.BC[4]:
        #                             '''
        #                             Does current sphere intersect current Dirichlet boundary?
        #                             '''
        #                             surf = surfaces[surface_number[0]]
        #                             surface = LineString(surf)
        #                             coorI = Point(self.input_data.node_coor[I])
        #                             distI = coorI.distance(surface)
        #                             coorJ = Point(self.input_data.node_coor[J])
        #                             distJ = coorJ.distance(surface)

        #                             u_sx = self.input_data.BC[0][dirichlet_count]
        #                             u_sy = self.input_data.BC[1][dirichlet_count]
        #                             u_s = [u_sx, u_sy]
        #                             if distI < self.input_data.radius or distJ < self.input_data.radius:
        #                                 N = direction_cosines(surf)
        #                                 h_Im = shape_functions(I,m,gauss.points[counter],self.input_data)
        #                                 H_Im = np.array([[h_Im, 0],\
        #                                                  [0, h_Im]])
        #                                 h_Jn = shape_functions(J,n,gauss.points[counter],self.input_data)
        #                                 H_Jn = np.array([[h_Jn, 0],\
        #                                                  [0, h_Jn]])
                                        
        #                                 y2 = gauss.y2
        #                                 y1 = gauss.y1
        #                                 KU_ImJn[I][J][m][n] += ((y2-y1)/2)*wj*H_Im @ N @ self.C @ B_Jn + ((y2-y1)/2)*wj* np.transpose(B_Im) @ self.C @ np.transpose(N) @ H_Jn                                        
        #                                 KU_ImJn_almostOpen[I,J,2*m:2+2*m,2*n:2+2*n] = KU_ImJn[I,J,m,n]
        #                                 KU_ImJn_open[2*order*I:2*order+2*order*I,2*order*J:2*order+2*order*J] = KU_ImJn_almostOpen[I,J]
        #                                 dirichlet_count+=1


        #                         if I==J and m==n:
        #                             '''
        #                             Loop over surfaces with prescribed displacements
        #                             '''
        #                             dirichlet_count=0
        #                             neumann_count=0
        #                             for surface_number in self.input_data.BC[4]:
        #                                 '''
        #                                 Does current sphere intersect current Dirichlet boundary?
        #                                 '''
        #                                 surf = surfaces[surface_number[0]]
        #                                 surface = LineString(surf)
        #                                 coor = Point(self.input_data.node_coor[I])
        #                                 dist = coor.distance(surface)

        #                                 u_sx = self.input_data.BC[0][dirichlet_count]
        #                                 u_sy = self.input_data.BC[1][dirichlet_count]
        #                                 u_s = [u_sx, u_sy]
        #                                 if dist < self.input_data.radius:
        #                                     N = direction_cosines(surf)
        #                                     fHat_Im[I][m] -= ((y2-y1)/2)*wj* np.transpose(B_Im) @ self.C @ np.transpose(N) @ u_s
        #                                     dirichlet_count+=1
                                    
        #                             '''
        #                             Loop over surfaces with prescribed loads
        #                             '''
        #                             for surface_number in self.input_data.loads[4]:
        #                                 '''
        #                                 Does current sphere intersect current Neumann boundary
        #                                 '''
        #                                 surf = surfaces[surface_number[0]]
        #                                 surface = LineString(surf)
        #                                 coor = Point(self.input_data.node_coor[I])
        #                                 dist = coor.distance(surface)

        #                                 f_sx = self.input_data.loads[0][neumann_count]
        #                                 f_sy = self.input_data.loads[1][neumann_count]
        #                                 f_s = [f_sx, f_sy]
        #                                 if dist < self.input_data.radius:
        #                                     N = direction_cosines(surf)
        #                                     point = gauss.points[counter]
        #                                     h_Im = shape_functions(I,m,gauss.points[counter],self.input_data)
        #                                     H_Im = np.array([[h_Im, 0],\
        #                                                     [0, h_Im]])
        #                                     y2 = gauss.y2
        #                                     y1 = gauss.y1
        #                                     fHat_Im[I][m] += ((y2-y1)/2)*wj* H_Im @ f_s
        #                                     neumann_count+=1
    
        # f = f_Im + fHat_Im
        # f = np.ravel(f)
        
        # K_open = K_open - KU_ImJn_open

        return K, f


