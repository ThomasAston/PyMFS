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
        Ndof = 3
        Ndim = 2
        gauss = gauss_points(I,J,input_data)
        K_ImJn_nested = np.zeros((Ndof,Ndof,Ndim,Ndim)) # Nested version of K
        K_ImJn = np.zeros((Ndof*Ndim, Ndof*Ndim)) # Un-nested version of K

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
                    continue
                else:
                    '''
                    Otherwise, we must... ASSEMBLE STIFFNESS MATRIX!
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

        
        return K_ImJn