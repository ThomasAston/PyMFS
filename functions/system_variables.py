'''
PyMFS generation of variables required for solving algebraic system
'''

'''
To-do list:

- Maybe combine K and f loops to make faster...
'''

from tabnanny import check
import numpy as np
# from functions.input_data import input_data
from .C_mat import *
from .K_mat import *
from .f_vec import *
from .gauss_points import *


class system_variables:
    def __init__(self, input_data):
        self.input_data = input_data
        self.C = C(input_data)

        print('Generating algebraic system...')
        self.K, self.f, self.delete = self.assemble()
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
        Ndof = 4 # Number of DoFs per direction of each node
        Ndim = 2 # Number of physical dimensions of problem domain.
        
        K = np.zeros((Ntot*Ndof*Ndim,Ntot*Ndof*Ndim))
        KU = np.zeros((Ntot*Ndof*Ndim,Ntot*Ndof*Ndim))
        f = np.zeros((Ntot*Ndof,Ndim))
        
        '''
        Loop over nodes I and J (direct spheres and overlap regions)
        and evaluate stiffness matrix and force vector components.

        Force components are only evaluated for direct spheres whereas
        stiffness components include direct spheres AND overlap regions.
        '''
        delete_rows=[]
        delete_cols=[]
        for I in range(Ntot):
            print('Assembling f for node: ',I,'..............')
            f_I = f_vec(I,self).vec
            f[Ndof*I:Ndof*I+Ndof] = f_I
            for J in range(Ntot):
                print('Assembling K for nodes: ',I, 'and', J,'..............')
                K_IJ, KU_IJ, type = K_mat(I,J,self).mat
                if type == 'Interior':
                    K[Ndim*Ndof*I:Ndim*Ndof+Ndim*Ndof*I,Ndim*Ndof*J:Ndim*Ndof+Ndim*Ndof*J] = K_IJ
                    KU[Ndim*Ndof*I:Ndim*Ndof+Ndim*Ndof*I,Ndim*Ndof*J:Ndim*Ndof+Ndim*Ndof*J] = KU_IJ
                elif type == 'Dirichlet':
                    K[Ndim*Ndof*I:Ndim*Ndof+Ndim*Ndof*I,Ndim*Ndof*J:Ndim*Ndof+Ndim*Ndof*J] = K_IJ
                    KU[Ndim*Ndof*I:Ndim*Ndof+Ndim*Ndof*I,Ndim*Ndof*J:Ndim*Ndof+Ndim*Ndof*J] = KU_IJ
                    
                    checkxI = Ndim*Ndof*I not in delete_rows
                    checkyI = Ndim*Ndof*I+1 not in delete_rows
                    checkxJ = Ndim*Ndof*J not in delete_cols
                    checkyJ = Ndim*Ndof*J+1 not in delete_cols

                    if checkxI == True:
                        delete_rows.append(Ndim*Ndof*I)
                    if checkyI==True:
                        delete_rows.append(Ndim*Ndof*I + 1)
                    if checkxJ == True:
                        delete_cols.append(Ndim*Ndof*J)
                    if checkyJ==True:
                        delete_cols.append(Ndim*Ndof*J + 1)

            
        f = np.reshape(f, [Ntot*Ndof*Ndim])
        K = K - KU
        
        f = np.delete(f,delete_rows,0)
        K = np.delete(K,delete_rows,0)
        K = np.delete(K,delete_cols,1)
       

        return K, f, delete_rows


