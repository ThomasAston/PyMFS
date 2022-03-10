'''
PyMFS gauss point generation
'''

'''
To-do list:

- 
'''

import numpy as np

class gauss_points:
    def __init__(self, I, J,):
        self.I = I
        self.J = J

        self.int_degree=6 # Control the number of integration points using this variable
        
        self.points, self.weights = self.generate()

    def generate(self): 
