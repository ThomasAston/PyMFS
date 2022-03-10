'''
PyMFS elasticity matrix generation
'''

'''
To-do list:

- Shape function generation
'''

import numpy as np

class D_mat:
    def __init__(self, material):
        self.material = material
        self.D = self.generate_D()

        # print(self.D)

    def generate_D(self):
        E = self.material.E
        v = self.material.v

        d11 = (E)/(1-v**2)
        d12 = (E*v)/(1-v**2)
        d33 = E/(2*(1+v))

        D = np.array([[d11, d12, 0],
                     [d12, d11, 0],
                     [0, 0, d33]])

        return D

