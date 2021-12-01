# Let's generate the system DoFs by considering boundary conditions 

# To do list:
# Assemble the load vector in the same order as the degrees of freedom...
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

        self.order, self.q, self.DoF_num, self.DoC_num = self.assign_BCs()
        self.f = self.assign_loads()


    # Function which flags nodes subject to Dirichlet condition
    def assign_BCs(self):
        BC_flag = np.ones((self.nodes.num, 2))
    
        # Generate the BC flag
        if self.BC == 'left fixed':
            for i in range(0, self.nodes.num):
                if self.nodes.coor[i,0] == 0:
                    BC_flag[i,:] = [0,0]
                    # DOF+=1
                    # U_p.append(0)
                else:
                    BC_flag[i,:] = [1,1]
        
        BC_flag = np.ravel(BC_flag)

        # Number the system DoFs, starting with free DoFs... 
        DoFs = np.zeros((len(BC_flag),1))
        q = []
        
        DoF_num = 0
        # Loop over the nodes, if node is free then start numbering the DoFs
        for i in range(0,len(BC_flag)):
            if BC_flag[i]==1:
                DoFs[i] = DoF_num
                DoF_num += 1
                q.append(0)
            else: continue

        # Loop over the nodes again, this time numbering only the restricted DoFs
        DoC_num = DoF_num
        for i in range(0,len(BC_flag)):
            if BC_flag[i]==0:
                DoFs[i] = DoC_num
                DoC_num += 1
                q.append(0)
            else: continue

        # print(DoFs)
        q = np.vstack([q]).reshape(-1,1)
        
        return DoFs, q, DoF_num, DoC_num


    # Function which assigns applied loads to relevant nodes
    def assign_loads(self):
        f = np.zeros((self.nodes.num*2, 1))

        order =self.order.astype(int)
        
        for i in range(0, self.nodes.num*2, 2):
            if self.nodes.coor[int(i/2),0] == self.location[0] and self.nodes.coor[int(i/2),1] == self.location[1]:
                f[order[i]] = self.load[0] 
                f[order[i]+1] = self.load[1]
            else:
                f[order[i]] = 0
                f[order[i]+1] = 0

        # self.order
        return f