'''
PyMFS input data unpacking
'''

'''
To-do list:

- Create input data object which can be passed wherever necessary
'''

import json

class input_data:

    def __init__(self, job_ID):
        self.surfaces, self.node_coor, self.radius, self.mat_prop, self.BC, self.loads= self.unpack()

    '''
    Function which unpacks the input .mfs file. 
    '''
    def unpack(self):
        with open(self.input_file,'r') as f:
            previous_line=''
            for line in f:
                if 'Surfaces' in previous_line:
                    surfaces = json.loads(line)
                elif 'Nodal coordinates' in previous_line:
                    node_coor = json.loads(line)
                elif 'Sphere radius' in previous_line:
                    radius = json.loads(line)
                elif 'Physical properties' in previous_line:
                    mat_prop = json.loads(line)
                elif 'Boundary conditions' in previous_line:
                    BC = json.loads(line)
                elif 'Loads' in previous_line:
                    loads = json.loads(line)
                
                previous_line=line
    
        return surfaces, node_coor, radius, mat_prop, BC, loads