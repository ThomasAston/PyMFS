''' 
PyMFS v1
Author(s): Thomas Aston
Date updated: 08/03/2022 

PyMFS has been developed as part of the MEng final year individual project:
'Develping a code for the method of finite spheres'.  

This is the central file for PyMFS... 

This file should provide users with sufficient functionality to generate PyMFS 
input files, then solve and analyse the results of the problem posed.  
'''

from functions import *
import time

'''
---------------------------------------------------------------------------------
------------------------ 1. Pre-processing --------------------------------------
---------------------------------------------------------------------------------
'''
'''
In the case that a .mfs job file has already been generated, this section can be 
commented out, and you may skip straight to section 2 (2. Solving).

The purpose of this section is to generate a .mfs job file which can be read by 
the PyMFS solver, with the purpose of solving solid mechanics problems using 
the method of finite spheres:
https://link.springer.com/article/10.1007/s004660050481
'''

'''
First, Set up the domain geometry with scripting commands.
For a detailed explanation please consult the PyMFS user guide.
'''
'''
Material addition:
'''
edge1 = straight_line(point1=[5, 0], point2=[10, 0])
edge2 = circular_segment(center=np.array([0,0]), radius=10, start=0, end=np.pi/2)
edge3 = straight_line(point1=[0, 10], point2=[0, 5])
edge4 = circular_segment(center=np.array([0,0]), radius=5, start=np.pi/2, end=0)
my_edges = (edge1, edge2, edge3, edge4)

'''
Material removal:
'''
subedge1 = circle(center=np.array([7,2]),radius=1) 
subedge2 = circle(center=np.array([2,7]),radius=1) 
my_subedges = (subedge1, subedge2)
# my_subedges = ()

'''
Generate a domain object from the given edges:
'''
my_domain = domain(my_edges, my_subedges)

'''
Discretise the domain by selecting number of nodes.
Note, at present sphere sizing is uniform only.
'''
my_nodes = nodes(my_domain, nx=11, ny=11, method='Regular')

'''
Enter the pre-processing UI to review geometry, set boundary conditions and
submit the job for solving:
'''
pre_process(my_domain, my_nodes, job_ID='job1')


'''
---------------------------------------------------------------------------------
--------------------------- 2. Solving ------------------------------------------
---------------------------------------------------------------------------------
'''
start = time.time()

'''
Select the input file to be solved and send it to PyMFS solver:
'''
solve(job_ID='job1.mfs')

end = time.time()
time_taken = start-end
print(time_taken)

'''
---------------------------------------------------------------------------------
------------------------ 3. Post-processing -------------------------------------
---------------------------------------------------------------------------------
'''
'''
Enter post-processing module. 
'''
