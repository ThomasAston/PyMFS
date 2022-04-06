'''
PyMFS boundary integration point generation
'''

'''
To-do list:

- Parametrise...
'''


import math
import numpy as np

def boundary_intpoints(I,J, surface, input_data):
    '''
    Read data: radius, current nodal coordinate, contributing
    nodal coordinate. Generate generalised gauss points.
    '''
    int_degree=6 # Choose degree of integration here
    r = input_data.radius
    input_data.node_coor = np.array(input_data.node_coor)
    node_x = input_data.node_coor[I,0]
    node_y = input_data.node_coor[I,1]
    contrib_x = input_data.node_coor[J,0]
    contrib_y = input_data.node_coor[J,1]

    x_int, w_int = np.polynomial.legendre.leggauss(int_degree)

    '''
    Determine what type of boundary we are on
    y2-y1 on surface=0, then surface is horizontal
    x2-x1 on surface=0, then surface is vertical
    otherwise, surface is curved.
    '''
    integration_points = []
    integration_weights = []

    if surface[1][1]-surface[0][1] == 0:
        t1 = max(0, node_x-r, contrib_x-r)
        t2 = min(2, node_x+r, contrib_x+r)
        # t_int = np.linspace(t1,t2,20)
        # for t in t_int:
        #     integration_points.append([t,node_y])
        if math.hypot(node_x-contrib_x, node_y-contrib_y)<2*r:
            for a in range(int_degree):
                integration_points.append([t1 + (t2-t1)/2 + x_int[a]*(t2-t1)/2,node_y])
                integration_weights.append(w_int[a])
        type = 'Horizontal'
    elif surface[1][0] - surface[0][0] == 0:
        t1 = max(0, node_y-r, contrib_y-r)
        t2 = min(2, node_y+r, contrib_y+r)
        # t_int = np.linspace(t1,t2,20)
        # for t in t_int:
        #     integration_points.append([node_x,t])
        if math.hypot(node_x-contrib_x, node_y-contrib_y)<2*r:
            for a in range(int_degree):
                integration_points.append([node_x,(t1 + (t2-t1)/2 + x_int[a]*(t2-t1)/2)])
                integration_weights.append(w_int[a])
        type = 'Vertical'
    
    integration_points = np.array(integration_points)
    integration_weights = np.array(integration_weights)

    return integration_points, integration_weights, t1, t2, type
    # return integration_points,t_int,type
