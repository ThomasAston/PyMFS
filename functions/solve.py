'''
PyMFS solver module
'''

'''
To-do list:

- Shape function generation
- Integration point generation
- C matrix
- B matrix
- K matrix
- solve system
'''

import numpy as np
from .input_data import*
from .K_mat import*
from .f_vec import*

class solve:
    def __init__(self, job_ID):
    
        self.input_data = input_data(job_ID)

        self.K = K_mat(self.input_data)
        self.f = f_vec(self.input_data)

    def solve_system(self):
        q = np.linalg.solve(self.K.open, self.f.open)






        


