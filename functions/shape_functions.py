'''
PyMFS shape function generation.

The purpose of this script is to return a matrix containing a 
matrix of shape function values at a given point, corresponding 
to a specific point and degree of freedom.
'''

'''
To-do list:

- 
'''

import numpy as np
import math

def shape_functions(node_num, DOF, point, input_data):
    Ntot = len(input_data.node_coor)
    node_coor = np.array(input_data.node_coor)
    r = input_data.radius

    x_i = node_coor[node_num, 0]
    y_i = node_coor[node_num, 1]
    x_point = point[0]
    y_point = point[1]

    s = math.hypot(x_point-x_i, y_point-y_i)/r
    
    if s<=1:
        W_i = 1-6*s**2+8*s**3-3*s**4
    else:
        W_i = 0
    
    sum_W_j = 0
    for j in range(Ntot):
        xj = node_coor[j,0]
        yj = node_coor[j,1]
        s_j = math.hypot(x_point-xj, y_point-yj)/r
        if s_j<=1:
            W_j=1-6*s_j**2+8*s_j**3-3*s_j**4
            sum_W_j+=W_j
        else:
            continue
    
    phi_i = W_i/sum_W_j 

    if DOF == 0:
        h = phi_i
    elif DOF == 1:
        h = phi_i*(x_point-x_i)/r
    elif DOF == 2:
        h = phi_i*(y_point-y_i)/r
    elif DOF ==3: 
        h = phi_i*(x_point-x_i)/r*(y_point-y_i)/r

    # H = [h, 0,\
    #      0, h]

    return h