# Let's generate some nodes on the problem domain...

# To do list:
# 
#
#

import numpy as np
from scipy.spatial import cKDTree


class nodes:
    def __init__(self, domain, N, dc = None, method = 'Regular'):
        self.domain = domain

        # If user chooses desired number of nodes
        if N:
            self.N = N
            self.dc = np.sqrt(self.domain.area)/(np.sqrt(N)-1)
        # Else if they choose desired nodal spacing
        elif dc:
            self.dc = dc
            self.N = self.domain.area/dc**2 + 2*self.domain.area/dc + 1
        else:
            raise ZeroDivisionError
        
        self.coor = np.zeros([self.N,2])

        if method == 'Regular':
            self.distribute_nodes_reg()
        
        self.tree = cKDTree(self.coor)

    # Function for distributing nodes over domain with regular spacing
    def distribute_nodes_reg(self):
        bounds = self.domain.bounds()
        lx = bounds[1,0] - bounds[0,0]
        ly = bounds[1,1] - bounds[0,1]
        
        bounds_area = lx*ly
        
        area = self.domain.area
        N = self.N*bounds_area/area
        
        nx = np.sqrt(lx/ly*N)
        ny = np.sqrt(ly/lx*N)
        
        nx = int(nx)
        ny = int(ny)
        
        x = np.linspace(bounds[0,0], bounds[1,0], nx)
        y = np.linspace(bounds[0,1], bounds[1,1], ny)
        xv, yv = np.meshgrid(x, y)        
        self.coor = np.array(list(zip(xv.ravel(), yv.ravel())))
        self.num = len(self.coor)
        self.nx = nx
        self.ny= ny


    def plot(self):
        from matplotlib import pylab
        pylab.rcParams["font.family"] = "serif"
        pylab.rcParams["mathtext.fontset"] = "dejavuserif"
        self.domain.plot(True)

        pylab.rcParams["font.family"] = "serif"
        pylab.rcParams["mathtext.fontset"] = "dejavuserif"
        
        pylab.scatter(self.coor[:,0], 
                    self.coor[:,1],
                    s = 4,
                    color = 'k')
        pylab.show()