'''
PyMFS elasticity matrix generation
'''

'''
To-do list:

- add flag for differentiating between plane stress and strain
'''

import numpy as np

def C(input_data):

    E = input_data.mat_prop[0]
    nu = input_data.mat_prop[1]

#     c_11 = E/(1-nu**2)
#     c_12 = (E*nu)/(1-nu**2)
#     c_33 = E/(2*(1+nu))

    c_11 = E*(1-nu)/((1+nu)*(1-2*nu))
    c_12 = E*nu/((1+nu)*(1-2*nu))
    c_33 = E/(2*(1+nu))


    C = np.array([[c_11, c_12, 0],\
         [c_12, c_11, 0],\
         [0   , 0   , c_33]])

    return C



