'''
PyMFS stiffness matrix generation

The purpose of this script is to create an object containing the system
stiffness matrix to be passed to the PyMFS solver.
'''

'''
To-do list:

- Possible move K and f evaluation to their own files.
- Body forces not currently included
'''

import numpy as np
from functions.input_data import input_data
from .B_mat import *
from .C_mat import *
from .gauss_points import *
from .direction_cosines import*
from shapely.geometry import LineString, Point

class system_variables:
    def __init__(self, input_data):
        self.input_data = input_data
        self.C = C(input_data)

        print('Generating stiffness matrix...')
        self.K, self.f = self.assemble()
        print('Stiffness matrix assembled.')
        

    def assemble(self):
        Ntot = len(self.input_data.node_coor)
        order = 3
        
        KU_ImJn = np.zeros((Ntot,Ntot,order,order,2,2)) # BC stiffness matrix
        KU_ImJn_almostOpen = np.zeros(((Ntot,Ntot,order*2,order*2)))
        KU_ImJn_open = np.zeros((Ntot*order*2,Ntot*order*2))

        K_nested = np.zeros(((Ntot,Ntot,order,order,2,2))) # Nested stiffness matrix
        K_almostOpen = np.zeros(((Ntot,Ntot,order*2,order*2)))
        K_open = np.zeros((Ntot*order*2,Ntot*order*2)) # Unravelled stiffness matrix

        f_Im = np.zeros((Ntot,order,2))  # Force vector
        fHat_Im = np.zeros((Ntot,order,2))  # Boundary force vector

        for I in range(Ntot):
            for J in range(Ntot):
                '''
                Generate gauss points and weights for current node
                and contributing node.
                '''
                I_coor = self.input_data.node_coor[I]
                J_coor = self.input_data.node_coor[J]
                gauss = gauss_points(I,J,self.input_data)
                
                '''
                Generate list of all surfaces, to be indexed to determine which
                surfaces BCs are applied to.
                '''
                surfaces = self.input_data.external_surfaces
                surfaces += self.input_data.internal_surfaces
                # dirichlet_count=0
                # neumann_count=0
                '''Loop over gauss points to evaluate B matrix'''
                for m in range(order):
                    for n in range(order):
                        '''
                        If no integration points exist in region, stiffness must be zero
                        '''
                        if np.array(gauss.points).size==0:
                            K_nested[I][J][m][n] = 0
                        else:
                            '''
                            STIFFNESS MATRIX ASSEMBLY
                            Evaluate the stiffness term by integrating over the
                            integration points.
                            '''
                            '''
                            FORCE VECTOR ASSEMBLY
                            If I=J and m=n, we can also evaluate f_Im inside this loop!
                            '''
                            for counter in range(len(gauss.points)):
                                B_Im = B(I,m,gauss.points[counter],self.input_data).mat
                                B_Jn = B(J,n,gauss.points[counter],self.input_data).mat

                                wi = np.array(gauss.weights)[counter,0]
                                wj = np.array(gauss.weights)[counter,1]

                                dx = self.input_data.radius
                                K_nested[I][J][m][n][:][:] += (dx/2)*(dx/2)*wj*wi*np.transpose(B_Im) @ self.C @ B_Jn
                                K_almostOpen[I,J,2*m:2+2*m,2*n:2+2*n] = K_nested[I,J,m,n]
                                K_open[2*order*I:2*order+2*order*I,2*order*J:2*order+2*order*J] = K_almostOpen[I,J]

                                dirichlet_count=0
                                for surface_number in self.input_data.BC[4]:
                                    '''
                                    Does current sphere intersect current Dirichlet boundary?
                                    '''
                                    surf = surfaces[surface_number[0]]
                                    surface = LineString(surf)
                                    coorI = Point(self.input_data.node_coor[I])
                                    distI = coorI.distance(surface)
                                    coorJ = Point(self.input_data.node_coor[J])
                                    distJ = coorJ.distance(surface)

                                    u_sx = self.input_data.BC[0][dirichlet_count]
                                    u_sy = self.input_data.BC[1][dirichlet_count]
                                    u_s = [u_sx, u_sy]
                                    if distI < self.input_data.radius or distJ < self.input_data.radius:
                                        N = direction_cosines(surf)
                                        h_Im = shape_functions(I,m,gauss.points[counter],self.input_data)
                                        H_Im = np.array([[h_Im, 0],\
                                                         [0, h_Im]])
                                        h_Jn = shape_functions(J,n,gauss.points[counter],self.input_data)
                                        H_Jn = np.array([[h_Jn, 0],\
                                                         [0, h_Jn]])
                                        
                                        y2 = gauss.y2
                                        y1 = gauss.y1
                                        KU_ImJn[I][J][m][n] += ((y2-y1)/2)*wi*H_Im @ N @ self.C @ B_Jn + (dx/2)*(dx/2)*wj*wi* np.transpose(B_Im) @ self.C @ np.transpose(N) @ H_Jn                                        
                                        KU_ImJn_almostOpen[I,J,2*m:2+2*m,2*n:2+2*n] = KU_ImJn[I,J,m,n]
                                        KU_ImJn_open[2*order*I:2*order+2*order*I,2*order*J:2*order+2*order*J] = KU_ImJn_almostOpen[I,J]
                                        dirichlet_count+=1


                                if I==J and m==n:
                                    '''
                                    Loop over surfaces with prescribed displacements
                                    '''
                                    dirichlet_count=0
                                    neumann_count=0
                                    for surface_number in self.input_data.BC[4]:
                                        '''
                                        Does current sphere intersect current Dirichlet boundary?
                                        '''
                                        surf = surfaces[surface_number[0]]
                                        surface = LineString(surf)
                                        coor = Point(self.input_data.node_coor[I])
                                        dist = coor.distance(surface)

                                        u_sx = self.input_data.BC[0][dirichlet_count]
                                        u_sy = self.input_data.BC[1][dirichlet_count]
                                        u_s = [u_sx, u_sy]
                                        if dist < self.input_data.radius:
                                            N = direction_cosines(surf)
                                            fHat_Im[I][m] -= (dx/2)*(dx/2)*wj*wi* np.transpose(B_Im) @ self.C @ np.transpose(N) @ u_s
                                            dirichlet_count+=1
                                    
                                    '''
                                    Loop over surfaces with prescribed loads
                                    '''
                                    for surface_number in self.input_data.loads[4]:
                                        '''
                                        Does current sphere intersect current Neumann boundary
                                        '''
                                        surf = surfaces[surface_number[0]]
                                        surface = LineString(surf)
                                        coor = Point(self.input_data.node_coor[I])
                                        dist = coor.distance(surface)

                                        f_sx = self.input_data.loads[0][neumann_count]
                                        f_sy = self.input_data.loads[1][neumann_count]
                                        f_s = [f_sx, f_sy]
                                        if dist < self.input_data.radius:
                                            N = direction_cosines(surf)
                                            
                                            h_Im = shape_functions(I,m,gauss.points[counter],self.input_data)
                                            H_Im = np.array([[h_Im, 0],\
                                                            [0, h_Im]])
                                            y2 = gauss.y2
                                            y1 = gauss.y1
                                            fHat_Im[I][m] += ((y2-y1)/2)*wi* H_Im @ f_s
                                            neumann_count+=1
    
        f = f_Im + fHat_Im
        f = np.ravel(f)
        
        K_open = K_open - KU_ImJn_open

        return K_open, f


