# Let's define the material used in the problem domain

# To do list:
# 
# 
#


import numpy as np

##########################################################################
# Class for defining the material of the domain
##########################################################################
class material:
    def __init__(self, E, v):
        self.properties = [E,v]
        self.E = E
        self.v = v

    def display_properties(self):
        print("Young's modulus: ",self.E)
        print("Poisson's ratio:  ",self.v)

##########################################################################
# Class for defining an elastic material
##########################################################################
# class elastic(material):
#     def __init__(self):
#         E, v = self.properties()
#         material.__init__(self, E, v)
        
#     def properties(self):
#         E = input("Young's modulus: ")
#         v = input("Poisson's ratio: ")
#         return E, v