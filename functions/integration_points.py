# Let's generate integration points for each region...

# To do list:
# Remove integration points on left side of domain
# 
#

import numpy as np

##########################################################################
# Class for generating the coords of the system's integration points
##########################################################################
class integration_points:
    def __init__(self, domain, nodes, elements, degree):
        self.elements = elements
        self.domain = domain
        self.nodes = nodes
        self.degree = degree # degree x degree gauss points in a disk quadrant
        self.num_points = degree**2 # number of gauss points per quadrant 
        self.gauss_points = self.generate_points()
        
    
    # Function which generates gaussian integration points.
    def generate_points(self):

        x_i, w_i = np.polynomial.legendre.leggauss(self.degree)
        xv, yv = np.meshgrid(x_i, x_i)
        self.gauss_grid = np.array(list(zip(xv.ravel(), yv.ravel())))
        self.gauss_grid = self.nodes.num*[self.gauss_grid]
        self.gauss_points = []


        for i in range(0,self.nodes.num):
            current_x = self.nodes.coor[i,0]
            current_y = self.nodes.coor[i,1]
            bounds = self.domain.bounds()

            if current_x == bounds[0,0] or current_y == bounds[1,1]:
                continue
            else:
                for j in range(0,self.num_points):
                    x = current_x - (0.5*self.elements.size) + (0.5*self.elements.size*self.gauss_grid[i-1][j,0])
                    y = current_y + (0.5*self.elements.size) + (0.5*self.elements.size*self.gauss_grid[i-1][j,1])
                    gauss_coor = [x, y]
                    self.gauss_points.append(gauss_coor)
                

        self.gauss_points = np.array(self.gauss_points)


        # Plotting integration points on the domain
        from matplotlib import pylab
        from matplotlib.lines import Line2D
        from matplotlib.patches import Polygon

        color = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
        
        fig, ax = pylab.subplots()
        self.domain.polygon.set_alpha(0.3)
        
        x, y = zip(*self.domain.polygon.xy)
        line = Line2D(x, y, linestyle= '-')
        ax.add_line(line)
        ax.set_title('Computational Domain')
        ax.set_xlim((min(x)-0.1*max(x), 1.1*max(x)))
        ax.set_ylim((min(y)-0.1*max(y), 1.1*max(y)))
        line.set_color('k')
         
        pylab.scatter(self.nodes.coor[:,0], 
            self.nodes.coor[:,1],
            s = 10,
            color = 'k')
    
        spheres=[]
        for i in range(0, len(self.nodes.coor)):
            spheres.append(pylab.Circle((self.nodes.coor[i,0],self.nodes.coor[i,1]), self.elements.size, fill=False, color='red'))
            ax.add_artist(spheres[i])

        pylab.scatter(self.gauss_points[:,0],
            self.gauss_points[:,1],
            s=4,
            color = 'g')
        

        pylab.show()


        

        # return gauss_points

                
    # def plot_points(self):
 
                


        