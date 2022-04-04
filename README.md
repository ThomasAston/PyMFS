[![forthebadge made-with-python](http://ForTheBadge.com/images/badges/made-with-python.svg)](https://www.python.org/)

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) 

# PyMFS

<a href="https://github.com/ThomasAston/PyMFS/graphs/contributors">
  <img src="https://contrib.rocks/image?repo=ThomasAston/PyMFS" />
</a>

Made with [contrib.rocks](https://contrib.rocks).

## Table of contents
* [About](#about)
* [Technologies](#technologies)
* [Setup](#setup)
* [To Do](#to_do)
* [Getting Involved](#getting_involved)


## About
PyMFS is a Python implementation of the method of finite spheres. The method of finite spheres (MFS) is a _truly_ meshless numerical method developed by Suvran De and Klaus-Jürgen Bathe [[1]]. The MFS has here been implemented in the first open-source code which adopts the method.
	
## Technologies
Project is created with:
* Python 3.10.0
	
## Setup
To run this project, enter

```
$ cd ../PyMFS
$ python3 PyMFS.py
```

The file that is run and the materials which are used are currently set in PyMFS.py by the following lines

```
edge1 = straight_line(point1=[0,0], point2=[2,0])
edge2 = straight_line(point1=[2,0], point2=[2,2])
edge3 = straight_line(point1=[2,2], point2=[0,2])
edge4 = straight_line(point1=[0,2], point2=[0,0])
my_edges = [edge1, edge2, edge3, edge4]
# '''
# Material removal:
# '''
# # subedge1 = circle(center=np.array([7,2]),radius=1) 
# # subedge2 = circle(center=np.array([2,7]),radius=1) 
# # my_subedges = [subedge1, subedge2]
my_subedges = []

'''
Generate a domain object from the given edges:
'''
my_domain = domain(my_edges, my_subedges)

'''
Discretise the domain by selecting number of nodes.
Note, at present sphere sizing is uniform only.
'''
my_nodes = nodes(my_domain, nx=3, ny=3, method='Regular')

'''
Enter the pre-processing UI to view geometry, set boundary conditions and
submit the job for solving:
'''
pre_process(my_domain, my_nodes, job_ID='job1')
```
## To do
- [ ] Complete solver for linear elasticity problems.
- [ ] Add tests

## Getting Involved
For any suggestions, please [create a new issue](https://github.com/ThomasAston/PyMFS/issues).

## Citations
1. De, S., Bathe, K. The method of finite spheres. Computational Mechanics 25, 329–345 (2000). 
 
[1]: https://doi.org/10.1007/s004660050481

