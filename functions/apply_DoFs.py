# Let's generate the system DoFs by considering boundary conditions 

# To do list:
# This whole thing can be much more robust.
# 
#

import numpy as np


##########################################################################
# Class for generating the system degrees of freedom based on applied
# boundary conditions and loads.
##########################################################################
class apply_DoFs:
    def __init__(self, nodes, BC, load, location):
        self.nodes = nodes
        self.BC = BC
        self.load = np.array(load)
        self.location = location


        self.assign_BCs()
        self.assign_loads()

    # Function which flags nodes subject to Dirichlet condition
    def assign_BCs(self):
        BC_flag = np.ones((self.nodes.num, 2))

        if self.BC == 'left fixed':
            for i in range(0, self.nodes.num):
                if self.nodes.coor[i,0] == 0:
                    BC_flag[i,:] = [0,0]

    # Function which assigns applied loads to relevant nodes
    def assign_loads(self):
        F_p = np.zeros((self.nodes.num, 2))
   
        for i in range(0, self.nodes.num):
            if self.nodes.coor[i,0] == self.location[0] and self.nodes.coor[i,1] == self.location[1]:
                F_p[i,0] = self.load[0] 
                F_p[i,1] = self.load[1]