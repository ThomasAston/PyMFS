# Let's generate the B matrix for a given node

# To do list:
# Should probably do some tidying up here
# Maybe move shepard_PU to its own class entirely???
# Make more efficient? Calculating derivatives currently takes ages...

import math
import numpy as np
from numpy.lib import math
from .domain import*

class B_mat:
    def __init__(self, domain, nodes, elements, x, y, current_node):
        self.domain = domain
        self.nodes = nodes
        self.elements = elements
        self.x = x
        self.y = y
        self.current_node = current_node

        self.phi = self.shepard_PU(x,y)
        self.value = self.generate_BMat(x,y)

    # Plotting an example...
    def example_plot(self):
        from matplotlib import pylab

        pylab.rcParams["font.family"] = "serif"
        pylab.rcParams["mathtext.fontset"] = "dejavuserif"

        ax = pylab.axes(projection='3d')
        i = self.current_node

        x = np.linspace(self.nodes.coor[i,0] - self.elements.size,self.nodes.coor[i,0] + self.elements.size,100)
        y = np.linspace(self.nodes.coor[i,1] - self.elements.size,self.nodes.coor[i,1] + self.elements.size,100)
        X, Y = np.meshgrid(x,y)

        phi = np.zeros([len(x), len(y)])
        for i in range(0,len(x)):
            for j in range(0,len(y)):
                phi[i,j] = self.shepard_PU(x[i], y[j])
        
        pylab.scatter(self.nodes.coor[:,0], 
                      self.nodes.coor[:,1],
                      s = 10,
                      color = 'k')
        # ax.plot_surface(X, Y, phi, rstride=1, cstride=1,
        #         cmap='viridis', edgecolor='none')
        ax.contour(X,Y,phi,50,cmap='viridis')
        ax.set_title(r'Weight function $\varphi$')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.zaxis.set_rotate_label(False)
        ax.set_zlabel(r'$\varphi$', rotation = 0);
        pylab.grid(False)
        # pylab.axis('off')

        pylab.show()
        ################################################################################################
        ax2 = pylab.axes(projection='3d')
        i = self.current_node

        dPHI_dx = np.zeros([len(x), len(y)])
        dPHI_dy = np.zeros([len(x), len(y)])
        for i in range(0,len(x)):
            for j in range(0,len(y)):
                dPHI_dx[i,j], dPHI_dy[i,j] = self.PHI(x[i], y[j])
        
        pylab.scatter(self.nodes.coor[:,0], 
                      self.nodes.coor[:,1],
                      s = 10,
                      color = 'k')
        # ax.plot_surface(X, Y, phi, rstride=1, cstride=1,
        #         cmap='viridis', edgecolor='none')
        ax2.contour(X,Y,dPHI_dy,50,cmap='viridis')
        ax2.set_title(r'Shape function derivatives')
        ax2.set_xlabel('x')
        ax2.set_ylabel('y')
        ax2.zaxis.set_rotate_label(False)
        ax2.set_zlabel(r'$\frac{d\phi}{dx}$', rotation = 0);
        # pylab.grid(False)
        # pylab.axis('off')

        pylab.show()

    # Function for generating the B matrix
    def generate_BMat(self, x, y):
        dPHI_dx, dPHI_dy = self.PHI(x,y)

        B_mat = np.array([[dPHI_dx, 0],
                          [0, dPHI_dy],
                          [dPHI_dx, dPHI_dy]])
        
        # print(B_mat)
        return B_mat


    # Function for evaluating shape function at a point
    def PHI(self, x, y):
        order = 0
        p = self.polynomial_basis(order, x, y)  
        PHI = self.shepard_PU(x,y)*p

        i = self.current_node
        x_i = self.nodes.coor[i][0]
        y_i = self.nodes.coor[i][1]

        dx = 0.1
        dy = 0.1
        PHI_right = self.shepard_PU(x+dx,y)*p 
        PHI_left = self.shepard_PU(x-dx,y)*p 
        PHI_up = self.shepard_PU(x,y+dy)*p 
        PHI_down = self.shepard_PU(x,y-dy)*p 

        dPHI_dx = (PHI_right - PHI_left)/(2*dx)
        dPHI_dy = (PHI_up-PHI_down)/(2*dy)
        
        return dPHI_dx, dPHI_dy


    # Function for calculating Shepard partition of unity weight at a given point 
    def shepard_PU(self, x, y): 
        i = self.current_node

        pos = np.array((x,y))
        x_i = self.nodes.coor[i][0]
        y_i = self.nodes.coor[i][1]
        pos_i = np.array((x_i, y_i))

        # s = abs(np.linalg.norm(pos-pos_i))/self.elements.size
        s = math.hypot(x-x_i, y-y_i)/self.elements.size

        # Quartic spline function
        if s<=1:
            W_i = 1 - 6*s**2 + 8*s**3 - 3*s**4
        else:
            W_i = 0
        
        W_j = []
        for j in range(0,self.nodes.num):
            x_j = self.nodes.coor[j][0]
            y_j = self.nodes.coor[j][1]
            pos_j = np.array((x_j,y_j))
            
            # s_j = np.linalg.norm(pos-pos_j)/self.elements.size
            
            s_j = math.hypot(x_i-x_j, y_i-y_j)/self.elements.size

            if s_j<=1:
                W_j.append(1 - 6*s_j**2 + 8*s_j**3 - 3*s_j**4)
            else:
                W_j.append(0)

        
        # Calculate value of shepard PU
        if sum(W_j) == 0:
            phi = 0
        else:
            phi = W_i/sum(W_j)

        return phi
    

    # Function for selecting the polynomial basis
    def polynomial_basis(self, order, x, y):
        if order==0:
            p = 1
        if order==1:
            p = np.array([1, x, y])
        if order==2:
            p = np.array([1, x, y, x**2, y**2, x*y])
        
        return p

