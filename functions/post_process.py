'''
PyMFS post-processing module

Takes a solution object past into it from the solver and reassembles
function fields before plotting chosen results.
'''

from .shape_functions import *
import matplotlib
from matplotlib import pyplot as plt

class post_process:
    def __init__(self, solution):
        self.solution = solution
        self.u = self.assemble_displacement()
        print(self.u)

    '''
    Function for re-assembling the nodal displacements from the obtained
    solution degrees of freedom
    '''
    def assemble_displacement(self):
        Ntot = len(self.solution.input_data.node_coor)
        order = 3
        dimensions = 2
        u= np.zeros([Ntot,2])
        q = np.reshape(self.solution.q, [Ntot, order, 2])
        for J in range(Ntot):
            for n in range(order):
                for dim in range(dimensions):
                    h_Jn = shape_functions(J,n,self.solution.input_data.node_coor[J],self.solution.input_data)
                    H_Jn = np.array([[h_Jn, 0],\
                                     [0, h_Jn]])
                    u[J, dim] += h_Jn*q[J,n, dim]

        deformed_shape = self.solution.input_data.node_coor+u*10
        
        fig,ax = plt.subplots()
        fig.canvas.manager.set_window_title('PyMFS post-processing') 
        
        ax.scatter(np.array(self.solution.input_data.node_coor)[:,0], np.array(self.solution.input_data.node_coor)[:,1], color='k')
        ax.scatter(deformed_shape[:,0], deformed_shape[:,1], color='r')
        plt.show()
        
        return u
                