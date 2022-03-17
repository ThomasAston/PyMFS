'''
PyMFS boundary integration point generation
'''

'''
To-do list:

- Add comments explaining steps
- Add general functionality for ALL boundary types
'''


import math
import numpy as np

def boundary_intpoints(I,J, input_data):
    int_degree=6

    r = input_data.radius
    input_data.node_coor = np.array(input_data.node_coor)
    node_x = input_data.node_coor[I,0]
    node_y = input_data.node_coor[I,1]
    contrib_x = input_data.node_coor[J,0]
    contrib_y = input_data.node_coor[J,1]

    x_int, w_int = np.polynomial.legendre.leggauss(int_degree)

    y1 = max(0, node_y-r, contrib_y-r)
    y2 = min(2, node_y+r, contrib_y+r)

    integration_points = []
    integration_weights = []

    if math.hypot(node_x-contrib_x, node_y-contrib_y)<2*r:
        for a in range(int_degree):
            integration_points.append(y1 + (y2-y1)/2 + x_int[a]*(y2-y1)/2)
            integration_weights.append(w_int[a])
    
    integration_points = np.array(integration_points)
    integration_weights = np.array(integration_weights)

    return integration_points, integration_weights, y1, y2

