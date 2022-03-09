'''
PyMFS solver module
'''


'''
To-do list:

- Unpack input file to return properties
- Shape function generation
- Integration point generation
- C matrix
- B matrix
- K matrix
- solve system
'''

import numpy as np

class solve:
    def __init__(self, job_ID):
        self.input_file = job_ID
        self.surfaces, self.node_coor, self.radius, self.mat_prop, self.BC, self.loads = self.unpack(self.input_file)


    def unpack(self):
        open_file = open(self.input_file,'r')

        
        return 



        


