'''
PyMFS domain generation
'''

'''
To-do list:

- Coordinate driven curves
'''

import numpy as np
from matplotlib.patches import Polygon
from matplotlib.patches import PathPatch
from matplotlib.patches import Path


'''
Class for defining a 2D polygonal domain made up of n vertices following
a closed path...

'''
class domain:
    # Function which takes the vertices passed into the domain and
    # constructs a polygon from it. 
    def __init__(self, edges, subedges):
        self.edges = edges
        self.vertices=[]
        for edge in edges:
            self.vertices += edge
        
        self.subedges = subedges
        self.subvertices = subedges
        
        self.polygon = Polygon(self.vertices)
        self.path = self.polygon.get_path()
        
        self.subPolygon=[]        
        self.subPath=[]
        if self.subvertices:
            for i in range(len(self.subvertices)):
                self.subPolygon.append(Polygon(self.subvertices[i]))
                self.subPath.append(self.subPolygon[i].get_path())

    # Function for calculating the area of the domain using Shoelace
    # formula: https://en.wikipedia.org/wiki/Shoelace_formula
    @property
    def area(self):
        x = self.vertices[:,0]
        y = self.vertices[:,1]

        return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

    # Returns a (minx, miny, maxx, maxy) tuple (float values) that bounds the object.
    def bounds(self):
        return np.array([np.min(self.vertices,axis=0), np.max(self.vertices,axis=0)])

    # Return true for points inside domain
    def contains(self, points):
        return np.array(self.path.contains_points(points))
    
    # Return true if all points inside domain
    def contains_all(self,points):
        return self.contains(points).all()

    # Function for plotting the domain
    def plot(self, show=True):
        
        from matplotlib import pylab
        from matplotlib.lines import Line2D
        pylab.rc('text', usetex=True)
        pylab.rc('text.latex', preamble=r'\usepackage{cmbright}')
        # pylab.rcParams["font.family"] = "serif"
        # pylab.rcParams["mathtext.fontset"] = "dejavuserif"

        color = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'] 
        
        fig, ax = pylab.subplots()
        self.polygon.set_alpha(0.3)
        x, y = zip(*self.polygon.xy)

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
        ax.add_patch(self.polygon)   

        for i in range(len(self.subPolygon)):
            x_sub, y_sub = zip(*self.subPolygon[i].xy)
            line_sub = Line2D(x_sub, y_sub, linestyle= '-')
            line_sub.set_color('k')
            ax.add_line(line_sub)
            self.subPolygon[i].set_color(color[7])
            # self.subPolygon.set_alpha(0.3)
            ax.add_patch(self.subPolygon[i])  
        if show:
            pylab.show()


'''
Function for defining a straight edge in the domain, using start and end
points.
'''
def straight_line(point1=[0, 0], point2=[1, 0]):
    point1_x = point1[0]
    point1_y = point1[1]
    point2_x = point2[0]
    point2_y = point2[1]

    if point2_x-point1_x ==0 : 
        gradient = 'undefined' # Line is vertical
        y_intercept = 'undefined'
    else:
        gradient = (point2_y-point1_y)/(point2_x-point1_x)
        y_intercept = point2_y-gradient*point2_x
    
    # vertices = []
    # if gradient == 'undefined':
    #     x = point1_x
    #     for y in np.linspace(point1_y,point2_y,100):
    #         vertices.append([x,y])
    # else: 
    #     for x in np.linspace(point1_x,point2_x,100):
    #         y = gradient*x + y_intercept
    #         vertices.append([x,y])
    # # vertices = []
    # # vertices.append(point1)
    # # # vertices = [[x,y] for x in np.linspace(point1_x,point2_x,100)[:-1] for y in np.linspace(point1_y,point2_y,100)[:-1]]
    # # for x in np.linspace(point1_x,point2_x,100):
    # #     for y in np.linspace(point1_y,point2_y,100):
    # #         vertices.append([x,y])
    
    # # vertices.append(point2)
    vertices = [point1,point2]


    return vertices

'''
Function for defining a circle in the domain, using centre point
and radius.
'''
def circle(center=np.array([0,0]), radius=1.0):
    resolution=100
    vertices = [[radius*np.cos(x)+center[0],radius*np.sin(x)+center[1]]
                for x in np.linspace(0,2*np.pi,resolution)[:-1]]
    
    return vertices

'''
Function for defining a circular segment in the domain, using centre
point, radius, and start and end angles relative to the horizontal axis
passing through the centre point. Anticlockwise direction is positive.
'''
def circular_segment(center=np.array([0,0]), radius=1.0, start=0, end=np.pi/2):
    resolution=100
    vertices = [[radius*np.cos(x)+center[0],radius*np.sin(x)+center[1]]
                for x in np.linspace(start,end,resolution)[:-1]]
    vertices.append([radius*np.cos(end)+center[0],radius*np.sin(end)+center[1]])
    return vertices
