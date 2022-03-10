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
from input_data import*
from K_mat import*
from f_vec import*

class solve:
    def __init__(self, job_ID):
    
        input_data = input_data(job_ID)

        self.K = K_mat(input_data)
        self.f = f_vec(input_data)






        


