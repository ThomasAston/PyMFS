# Let's post-process the obtained results...

# To do list:
# 
#
#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

class post:
    def __init__(self, solution, domain, nodes, elements, integration_points):
        self.nodes = nodes
        self.solution = solution

        self.plot_disp()

    def plot_disp(self):
        fig_1 = plt.figure(1) 
        plt.rcParams["font.family"] = "serif"
        plt.rcParams["mathtext.fontset"] = "dejavuserif"
        
        plt.title('Deformed shape')
        axdef = fig_1.add_subplot(111)
        
        x_undeformed = self.nodes.coor[:,0]
        y_undeformed = self.nodes.coor[:,1]

        x_deformed = self.nodes.coor[:,0] + 1000*self.solution.u[:,0]
        y_deformed = self.nodes.coor[:,1] + 1000*self.solution.u[:,1]

        plt.scatter(x_deformed, y_deformed, color = 'r')
        plt.scatter(x_undeformed,y_undeformed,color='k')

        count = 1 # Annotate node numbers
        for i in range(0, self.nodes.num):
            plt.annotate(count, xy = (x_deformed[i], y_deformed[i]), color='r')
            count +=1
        count = 1
        for i in range(0, self.nodes.num):
            plt.annotate(count, xy = (x_undeformed[i], y_undeformed[i]), color='k')
            count +=1

        plt.show()


    def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=-1):
        if n==-1:
            n=cmap.N
        new_cmap = mcolors.LinearSegmentedColormap.from_list('trunc({name},{a:2f},{b:2f})'.format(name=cmap.name, a=minval, b=maxval),cmap(np.linspace(minval,maxval,n)))
        return new_cmap