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
PyMFS is a Python implementation of the method of finite spheres. The method of finite spheres (MFS) is a _truly_ meshless numerical method developed by Suvran De and Klaus-Jürgen Bathe [1]. The MFS has here been implemented in the first open-source code which adopts the method.
	
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
# Input file and materials file
filename = 'Tinput.tro'
material = 'materials.json'
```
## To do
- [ ] Complete solver for linear elasticity problems.

## Getting Involved
For any suggestions, please [create a new issue](https://github.com/ThomasAston/PyMFS/issues).

## Citations
1. De, S., Bathe, K. The method of finite spheres. Computational Mechanics 25, 329–345 (2000). 
 
[1]: https://doi.org/10.1007/s004660050481

