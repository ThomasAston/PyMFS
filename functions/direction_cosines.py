'''
PyMFS direction cosine generation

The purpose of this script is to return the direction
cosine components at a given point on a chosen boundary. 
'''

'''
To-do list:

- Currently only handles straight lines, adapt to curves
'''

import numpy as np

def direction_cosines(surface):

    '''
    y2-y1 on surface=0, then surface is horizontal
    x2-x1 on surface=0, then surface is vertical
    '''

    if surface[1][1]-surface[0][1] == 0:
        n_x = 0
        n_y = 1
    elif surface[1][0] - surface[0][0] == 0:
        n_x = 1
        n_y = 0
    else:
        n_x=0
        n_y=0    

    N = np.array([[n_x, 0, n_y],\
                  [0, n_y, n_x]])

    return N