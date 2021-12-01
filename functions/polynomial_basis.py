# Let's generate the polynomial basis

# To do list:
# 

import numpy as np

class polynomial_basis:
    def __init__(self,order):
        self.order = order

        self.value = self.generate_basis(order,1,1)

    # Function for selecting the polynomial basis
    def generate_basis(self, order, x, y):
        if order==0:
            p = np.array([1])
        if order==1:
            p = np.array([1, x, y])
        if order==2:
            p = np.array([1, x, y, x**2, y**2, x*y])
        
        return p
