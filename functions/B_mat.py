'''
PyMFS B matrix generation
'''

'''
To-do list:

- 
'''

from functions.shape_functions import shape_functions
import numpy as np

class B:
    def __init__(self, node_num, DOF, gauss_point, input_data):
        self.node_num = node_num
        self.DOF = DOF
        self.point = gauss_point
        self.input_data = input_data
        
        self.mat = self.generate()

        
    def generate(self):
        point = self.point
        
        h = shape_functions(self.node_num,self.DOF, point, self.input_data)

        delta = 0.0001
        xplus = point[0]+delta
        yplus = point[1]+delta
        point_xplus = (xplus,point[1])
        point_yplus = (point[0],yplus)

        h_xplus = shape_functions(self.node_num,self.DOF, point_xplus, self.input_data)
        h_yplus = shape_functions(self.node_num,self.DOF, point_yplus, self.input_data)
        
        dh_dx = (h_xplus-h)/delta
        dh_dy = (h_yplus-h)/delta

        B_mat = np.array([[dh_dx, 0], \
                          [0, dh_dy], \
                          [dh_dy, dh_dx]])

        return B_mat