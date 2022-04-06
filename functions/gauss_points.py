'''
PyMFS gauss point generation
'''

'''
To-do list:

-
'''

import numpy as np
import math
from matplotlib.patches import Polygon

import matplotlib
from matplotlib import pylab
from matplotlib.lines import Line2D


class gauss_points:
    def __init__(self, I, J,input_data):
        self.I = I
        self.J = J
        self.input_data = input_data

        '''External polygon'''
        vertices=[]
        for surf in self.input_data.external_surfaces:
            vertices += surf
        self.polygon = Polygon(vertices)
        self.path = self.polygon.get_path()

        '''Internal polygons'''
        self.subPolygon=[]
        self.subPath=[]
        if self.input_data.internal_surfaces:
            for i in range(len(self.input_data.internal_surfaces)):
                self.subPolygon.append(Polygon(self.input_data.internal_surfaces[i]))
                self.subPath.append(self.subPolygon[i].get_path())

        '''Degree of gaussian quadrature per quadrant'''
        self.int_degree=6
        
        '''Generate points for current node I, contributing node J'''
        self.points, self.weights,self.y1,self.y2 = self.generate()
        
        pylab.rc('text', usetex=True)
        pylab.rc('text.latex', preamble=r'\usepackage{cmbright}')
        
        



    '''
    Function which returns gaussian points and weights which lie inside
    the domain for the current node I and contributing node J
    '''
    def generate(self): 
        coor = np.array(self.input_data.node_coor) 
        node_x = coor[self.I,0]
        node_y = coor[self.I,1]
        contrib_x = coor[self.J,0]
        contrib_y = coor[self.J,1]
        r = self.input_data.radius
        '''Generate gauss points and weights'''
        x_int, w_int = np.polynomial.legendre.leggauss(self.int_degree)
        xi, xj = np.meshgrid(x_int, x_int)
        wi, wj = np.meshgrid(w_int, w_int)

        gauss_grid = np.array(list(zip(xi.ravel(), xj.ravel())))
        weight_grid = np.array(list(zip(wi.ravel(), wj.ravel())))
        points = []
        weights = []
        xpoints = []
        
        y1 = max(0,node_y-r, contrib_y-r)
        y2 = min(2,node_y+r, contrib_y+r)
        
        '''If I==J we are directly integrating a sphere'''
        if self.I==self.J:
            '''Populate quadrant 1'''
            for a in range(self.int_degree**2):
                x = node_x - 0.5*r + 0.5*r*gauss_grid[a,0]
                y = node_y + 0.5*r + 0.5*r*gauss_grid[a,1]
                point = (x,y)
                pointInDomain = self.polygon.contains_point(point)
                if pointInDomain==True and math.hypot(x-node_x, y-node_y) < r:
                    pointInSubDomain = []
                    for b in range(len(self.subPolygon)):
                        if self.subPolygon[b].contains_point(point)==True:
                            pointInSubDomain.append(1)
                    check = (1 not in pointInSubDomain)
                    if check==True:
                        points.append([x,y])
                        weights.append([weight_grid[a,0], weight_grid[a,1]])
                    else:
                        xpoints.append([x,y])
                        continue
                elif pointInDomain==True and math.hypot(x-node_x, y-node_y) > r:
                    xpoints.append([x,y])
            '''Populate quadrant 2'''
            for a in range(self.int_degree**2):
                x = node_x - 0.5*r + 0.5*r*gauss_grid[a,0]
                y = node_y - 0.5*r + 0.5*r*gauss_grid[a,1]
                point = (x,y)
                pointInDomain = self.polygon.contains_point(point)
                if pointInDomain==True and math.hypot(x-node_x, y-node_y) < r:
                    pointInSubDomain = []
                    for b in range(len(self.subPolygon)):
                        if self.subPolygon[b].contains_point(point)==True:
                            pointInSubDomain.append(1)
                    check = (1 not in pointInSubDomain)
                    if check==True:
                        points.append([x,y])
                        weights.append([weight_grid[a,0], weight_grid[a,1]])
                    else:
                        xpoints.append([x,y])
                        continue
                elif pointInDomain==True and math.hypot(x-node_x, y-node_y) > r:
                    xpoints.append([x,y])
            '''Populate quadrant 3'''
            for a in range(self.int_degree**2):
                x = node_x + 0.5*r + 0.5*r*gauss_grid[a,0]
                y = node_y - 0.5*r + 0.5*r*gauss_grid[a,1]
                point = (x,y)
                pointInDomain = self.polygon.contains_point(point)
                if pointInDomain==True and math.hypot(x-node_x, y-node_y) < r:
                    pointInSubDomain = []
                    for b in range(len(self.subPolygon)):
                        if self.subPolygon[b].contains_point(point)==True:
                            pointInSubDomain.append(1)
                    check = (1 not in pointInSubDomain)
                    if check==True:
                        points.append([x,y])
                        weights.append([weight_grid[a,0], weight_grid[a,1]])
                    else:
                        xpoints.append([x,y])
                        continue
                elif pointInDomain==True and math.hypot(x-node_x, y-node_y) > r:
                    xpoints.append([x,y])
            '''Populate quadrant 4'''
            for a in range(self.int_degree**2):
                x = node_x + 0.5*r + 0.5*r*gauss_grid[a,0]
                y = node_y + 0.5*r + 0.5*r*gauss_grid[a,1]
                point = (x,y)
                pointInDomain = self.polygon.contains_point(point)
                if pointInDomain==True and math.hypot(x-node_x, y-node_y) < r:
                    pointInSubDomain = []
                    for b in range(len(self.subPolygon)):
                        if self.subPolygon[b].contains_point(point)==True:
                            pointInSubDomain.append(1)
                    check = (1 not in pointInSubDomain)
                    if check==True:
                        points.append([x,y])
                        weights.append([weight_grid[a,0], weight_grid[a,1]])
                    else:
                        xpoints.append([x,y])
                        continue
                elif pointInDomain==True and math.hypot(x-node_x, y-node_y) > r:
                    xpoints.append([x,y])
        else:
            '''
            Otherwise I=\=J, so we are attempting to integrate a sphere overlap.
            If spheres overlap, populate overlap region.
            '''
            if math.hypot(node_x-contrib_x, node_y-contrib_y) < 2*r:
                '''Populate quadrant 1'''
                for a in range(self.int_degree**2):
                    x = node_x - 0.5*r + 0.5*r*gauss_grid[a,0]
                    y = node_y + 0.5*r + 0.5*r*gauss_grid[a,1]
                    point = (x,y)
                    pointInDomain = self.polygon.contains_point(point)
                    if pointInDomain==True and math.hypot(x-node_x, y-node_y) < r and math.hypot(x-contrib_x, y-contrib_y) < r:
                        pointInSubDomain = []
                        for b in range(len(self.subPolygon)):
                            if self.subPolygon[b].contains_point(point)==True:
                                pointInSubDomain.append(1)
                        check = (1 not in pointInSubDomain)
                        if check==True:
                            points.append([x,y])
                            weights.append([weight_grid[a,0], weight_grid[a,1]])
                        else:
                            xpoints.append([x,y])
                            continue
                    elif pointInDomain==True and math.hypot(x-node_x, y-node_y) > r or pointInDomain==True and math.hypot(x-contrib_x,y-contrib_y)>r:
                        xpoints.append([x,y])
                '''Populate quadrant 2'''
                for a in range(self.int_degree**2):
                    x = node_x - 0.5*r + 0.5*r*gauss_grid[a,0]
                    y = node_y - 0.5*r + 0.5*r*gauss_grid[a,1]
                    point = (x,y)
                    pointInDomain = self.polygon.contains_point(point)
                    if pointInDomain==True and math.hypot(x-node_x, y-node_y) < r and math.hypot(x-contrib_x, y-contrib_y) < r:
                        pointInSubDomain = []
                        for b in range(len(self.subPolygon)):
                            if self.subPolygon[b].contains_point(point)==True:
                                pointInSubDomain.append(1)
                        check = (1 not in pointInSubDomain)
                        if check==True:
                            points.append([x,y])
                            weights.append([weight_grid[a,0], weight_grid[a,1]])
                        else:
                            xpoints.append([x,y])
                            continue
                    elif pointInDomain==True and math.hypot(x-node_x, y-node_y) > r or pointInDomain==True and math.hypot(x-contrib_x,y-contrib_y)>r:
                        xpoints.append([x,y])
                
                '''Populate quadrant 3'''
                for a in range(self.int_degree**2):
                    x = node_x + 0.5*r + 0.5*r*gauss_grid[a,0]
                    y = node_y - 0.5*r + 0.5*r*gauss_grid[a,1]
                    point = (x,y)
                    pointInDomain = self.polygon.contains_point(point)
                    if pointInDomain==True and math.hypot(x-node_x, y-node_y) < r and math.hypot(x-contrib_x, y-contrib_y) < r:
                        pointInSubDomain = []
                        for b in range(len(self.subPolygon)):
                            if self.subPolygon[b].contains_point(point)==True:
                                pointInSubDomain.append(1)
                        check = (1 not in pointInSubDomain)
                        if check==True:
                            points.append([x,y])
                            weights.append([weight_grid[a,0], weight_grid[a,1]])
                        else:
                            continue
                    elif pointInDomain==True and math.hypot(x-node_x, y-node_y) > r or pointInDomain==True and math.hypot(x-contrib_x,y-contrib_y)>r:
                        xpoints.append([x,y])
                '''Populate quadrant 4'''
                for a in range(self.int_degree**2):
                    x = node_x + 0.5*r + 0.5*r*gauss_grid[a,0]
                    y = node_y + 0.5*r + 0.5*r*gauss_grid[a,1]
                    point = (x,y)
                    pointInDomain = self.polygon.contains_point(point)
                    if pointInDomain==True and math.hypot(x-node_x, y-node_y) < r and math.hypot(x-contrib_x, y-contrib_y) < r:
                        pointInSubDomain = []
                        for b in range(len(self.subPolygon)):
                            if self.subPolygon[b].contains_point(point)==True:
                                pointInSubDomain.append(1)
                        check = (1 not in pointInSubDomain)
                        if check==True:
                            points.append([x,y])
                            weights.append([weight_grid[a,0], weight_grid[a,1]])
                        else:
                            continue
                    elif pointInDomain==True and math.hypot(x-node_x, y-node_y) > r or pointInDomain==True and math.hypot(x-contrib_x,y-contrib_y)>r:
                        xpoints.append([x,y])
                        
                        
        '''Test plot'''
        # fig,ax=pylab.subplots()
        
        # radius = 1
        # center = [1,1]
        # start = 0
        # end = np.pi*2
        # resolution=100
        # vertices = [[radius*np.cos(x)+center[0],radius*np.sin(x)+center[1]]
        #         for x in np.linspace(start,end,resolution)[:-1]]
        # # radius = 1
        # # center = [0,1]
        # # start = -np.pi/6
        # # end = -np.pi/2
        # # resolution=100
        # # vertices += [[radius*np.cos(x)+center[0],radius*np.sin(x)+center[1]]
        # #         for x in np.linspace(start,end,resolution)[:-1]]

        # # vertices.append([[radius*np.cos(x)+center[0],radius*np.sin(x)+center[1]] for x in np.linspace(start,end,resolution)[:-1]])
        # # vertices.append([0,0])
        # # vertices.append([0,1])
        # polygon = Polygon(vertices, facecolor=(0,0.541176,0.541176,0.3))
        # ax.add_patch(polygon) 
        
        # # if points:
        # #     pylab.scatter(np.array(points)[:,0], np.array(points)[:,1], s=10, color = [(0,0.541176,0.541176)])
        # # if xpoints:    
        # #     pylab.scatter(np.array(xpoints)[:,0], np.array(xpoints)[:,1], s=10, color = 'r')

        # pts = pylab.scatter(np.array(self.input_data.node_coor)[:,0], 
        #                     np.array(self.input_data.node_coor)[:,1],
        #                     s = 25,
        #                     color = 'k')

        # x,y = zip(*self.polygon.xy)
        # line = Line2D(x,y,linestyle='-', color='k', linewidth=1)
        # ax.set_xlim((min(x)-0.1*max(x), 1.1*max(x)))
        # ax.set_ylim((min(y)-0.1*max(y), 1.1*max(y)))
        # ax.add_line(line)
        # ax.set_aspect(aspect= 1)
        
        # sphereI = pylab.Circle((np.array(self.input_data.node_coor)[self.I,0],np.array(self.input_data.node_coor)[self.I,1]), self.input_data.radius, fill=False, color=(0.6627,0.6627,0.6627))
        # sphereJ = pylab.Circle((np.array(self.input_data.node_coor)[self.J,0],np.array(self.input_data.node_coor)[self.J,1]), self.input_data.radius, fill=False, color=(0.6627,0.6627,0.6627))
        # ax.add_patch(sphereI)
        # ax.add_patch(sphereJ)

        # pylab.show()

        return points, weights, y1, y2
