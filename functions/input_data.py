'''
PyMFS input data unpacking
'''

'''
To-do list:

- Create input data object which can be passed wherever necessary
'''

import json
import numpy as np
from shapely.geometry import LineString, Point

class input_data:
    def __init__(self, job_ID):
        self.input_file = job_ID
        self.external_surfaces, self.internal_surfaces, self.node_coor, \
            self.radius, self.mat_prop, self.BC, self.loads = self.unpack()
        sphere_types = self.sphere_types()
        # print(sphere_types)

    '''
    Function which unpacks the input .mfs file. 
    '''
    def unpack(self):
        with open(self.input_file,'r') as f:
            previous_line=''
            for line in f:
                if 'External surfaces' in previous_line:
                    external_surfaces = json.loads(line)
                elif 'Internal surfaces' in previous_line:
                    internal_surfaces = json.loads(line)    
                elif 'Nodal coordinates' in previous_line:
                    node_coor = json.loads(line)
                elif 'Sphere radius' in previous_line:
                    radius = json.loads(line)
                elif 'Physical properties' in previous_line:
                    mat_prop = json.loads(line)
                elif 'Prescribed displacements' in previous_line:
                    BC = json.loads(line)
                elif 'Prescribed loads' in previous_line:
                    loads = json.loads(line)
                
                previous_line=line
    
        return external_surfaces, internal_surfaces, node_coor, radius, mat_prop, BC, loads


    '''
    Function which flags each sphere, to determine
    whether it is INTERIOR, NEUMANN or DIRICHLET
    '''
    def sphere_types(self):
        
        '''
        Loop over every nodal coordinate.
        For each coordinate, loop over domain surfaces.
        If distance between node and surface < sphere radius, assign sphere type.
        '''
        sphere_types = []

        for coor in self.node_coor:
            coor = Point(coor)
            sphere_type=[]
            intersects = []
            for surf in range(len(self.internal_surfaces)):
                surface = self.internal_surfaces[surf]
                surface = LineString(surface)
                dist = coor.distance(surface)
                
                '''Does the current sphere intersect an external surface?'''
                if dist < self.radius:
                    intersects.append('YES')
                else:
                    intersects.append('NO')
            
    
            if 'YES' in intersects:
                sphere_type = 1
            else:
                sphere_type = 0

            sphere_types.append(sphere_type)

        return sphere_types
