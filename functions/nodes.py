'''
PyMFS node generation
'''

'''
To-do list:

- Irregular node distributions
'''

import numpy as np
from scipy.spatial import cKDTree


class nodes:
    def __init__(self, domain, nx, ny, dc = None, method = 'Regular'):
        self.domain = domain

        # If user chooses desired number of nodes
        # if N:
        #     self.N = N
        #     self.dc = np.sqrt(self.domain.area)/(np.sqrt(N)-1)
        # # Else if they choose desired nodal spacing
        # elif dc:
        #     self.dc = dc
        #     self.N = self.domain.area/dc**2 + 2*self.domain.area/dc + 1
        # else:
        #     raise ZeroDivisionError
        
        # self.coor = np.zeros([self.N,2])
        self.nx = nx
        self.ny = ny
        if method == 'Regular':
            self.distribute_nodes_reg(nx,ny)
        
        # self.tree = cKDTree(self.coor)

    # Function for distributing nodes over domain with regular spacing
    def distribute_nodes_reg(self,nx,ny):
        bounds = self.domain.bounds()

        nx = int(nx)
        ny = int(ny)
        
        x = np.linspace(bounds[0,0], bounds[1,0]-0.001, nx)
        y = np.linspace(bounds[0,1]+0.001, bounds[1,1], ny)
        # xv, yv = np.meshgrid(x, y)
        self.coor=[]
        for i in range(len(x)):
            for j in range(len(y)):
                point = (x[i],y[j])
                pointInDomain = self.domain.polygon.contains_point(point)
                if pointInDomain == True:
                    pointInSubDomain=[]
                    for k in range(len(self.domain.subPolygon)):
                        if self.domain.subPolygon[k].contains_point(point) == True:
                            pointInSubDomain.append(1)
                    check = (1 not in pointInSubDomain)
                    if check == True:
                        self.coor.append([x[i],y[j]])
        
        self.coor = np.array(self.coor)
        
        # self.coor = np.array(list(zip(xv.ravel(), yv.ravel())))
        self.num = len(self.coor)
        self.nx = nx
        self.ny= ny


    def plot(self):
        
        from matplotlib import pylab
        from matplotlib.lines import Line2D
        from matplotlib.collections import LineCollection
        
        pylab.rc('text', usetex=True)
        pylab.rc('text.latex', preamble=r'\usepackage{cmbright}')
        # pylab.rcParams["font.family"] = "serif"
        # pylab.rcParams["mathtext.fontset"] = "dejavuserif"

        color = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'] 
        
        fig, ax = pylab.subplots()
        self.domain.polygon.set_alpha(0.3)
        x, y = zip(*self.domain.polygon.xy)

        line = Line2D(x, y, linestyle= '-')
        ax.add_line(line)
        ax.set_title(r'Computational Domain')
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$y$')
        ax.set_xlim((min(x)-0.1*max(x), 1.1*max(x)))
        ax.set_ylim((min(y)-0.1*max(y), 1.1*max(y)))
        line.set_color('k')
        ax.grid(True)
        ax.set_aspect(aspect= 1)
        ax.add_patch(self.domain.polygon)  
        # ax.scatter(x,y,color='r') 

        for i in range(len(self.domain.subPolygon)):
            x_sub, y_sub = zip(*self.domain.subPolygon[i].xy)
            line_sub = Line2D(x_sub, y_sub, linestyle= '-')
            line_sub.set_color('k')
            ax.add_line(line_sub)
            self.domain.subPolygon[i].set_color(color[7])
            # self.subPolygon.set_alpha(0.3)
            ax.add_patch(self.domain.subPolygon[i]) 

        pylab.scatter(self.coor[:,0], 
                    self.coor[:,1],
                    s = 4,
                    color = 'k')
       
        spheres=[]
        size = self.coor[0,1] - self.coor[1,1]
        for i in range(0, len(self.coor)):
            spheres.append(pylab.Circle((self.coor[i,0],self.coor[i,1]), size, fill=False, color='c'))
            ax.add_artist(spheres[i])
        
        
        pylab.show()