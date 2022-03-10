'''
PyMFS gauss point generation
'''

'''
To-do list:

- This is currently extremely slow, figure out how to make faster.
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
        self.points, self.weights = self.generate()

        '''Test plot'''
        # fig,ax=pylab.subplots()
        # if self.points:
        #     pylab.scatter(np.array(self.points)[:,0], np.array(self.points)[:,1])

        # x,y = zip(*self.polygon.xy)
        # line = Line2D(x,y,linestyle='-')
        # ax.add_line(line)
        # pylab.show()


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
                        continue
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
                        continue
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
                        continue
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
                        continue
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
                            continue
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
                            continue
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

        return points, weights
