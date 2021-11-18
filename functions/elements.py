# Let's generate some elements from the node list...

# To do list:
# 
#
#


##########################################################################
# Class for generating elements for the list of nodes
##########################################################################
class elements:
    def __init__(self, nodes, domain, type = 'MFS', size = 0.1):
        self.nodes = nodes
        self.domain = domain
        self.size = size
        self.type = type
        if type=='MFS':
            self.generate_MFS()

    # Function for generating MFS nodes if MFS node type chosen...
    def generate_MFS(self):

        functional_form = 'quartic spline'

        
        # Plotting spheres on the domain
        from matplotlib import pylab
        from matplotlib.lines import Line2D
        from matplotlib.patches import Polygon

        color = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
        
        fig, ax = pylab.subplots()
        self.domain.polygon.set_alpha(0.3)
        
        x, y = zip(*self.domain.polygon.xy)
        line = Line2D(x, y, linestyle= '-')
        ax.add_line(line)
        ax.set_title('Computational Domain')
        ax.set_xlim((min(x)-0.1*max(x), 1.1*max(x)))
        ax.set_ylim((min(y)-0.1*max(y), 1.1*max(y)))
        line.set_color('k')
         
        pylab.scatter(self.nodes.coor[:,0], 
            self.nodes.coor[:,1],
            s = 4,
            color = 'k')
    
        spheres=[]
        for i in range(0, len(self.nodes.coor)):
            spheres.append(pylab.Circle((self.nodes.coor[i,0],self.nodes.coor[i,1]), self.size, fill=False, color='red'))
            ax.add_artist(spheres[i])
        

        pylab.show()