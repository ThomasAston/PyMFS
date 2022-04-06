'''
PyMFS pre-processing module
'''

'''
To-do list:

- Speed up using blitting
- Add identifier for if problem is plane stress or plane strain?
'''

import numpy as np
from matplotlib.widgets import CheckButtons, Button, TextBox

class pre_process:
    def __init__(self, domain, nodes, job_ID):
        self.domain = domain
        self.nodes = nodes
        self.job_ID = job_ID
        self.u_s, self.f_s, self.PhysicalProperties = self.init_window()
        self.compile_MFS()

    '''
    Open the pre-processing window and generate user-interface
    '''
    def init_window(self):
        ### Plot the geometry ###
        import matplotlib
        matplotlib.use("Qt5Agg")
        from matplotlib import pylab
        from matplotlib.lines import Line2D
        
        pylab.rc('text', usetex=True)
        pylab.rc('text.latex', preamble=r'\usepackage{cmbright}')

        color = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w'] 

        fig, ax = pylab.subplots()
        fig.canvas.manager.set_window_title('PyMFS pre-processing') 
        self.domain.polygon.set_alpha(0.3)
        self.domain.polygon.set_color((0,0.541176,0.541176))
        x, y = zip(*self.domain.polygon.xy)
        
        ### Plot outlines and add patches ###
        line = Line2D(x, y, linestyle= '-')
        ax.add_line(line)
        line.set_color('k')
        
        ax.add_patch(self.domain.polygon)  
        for i in range(len(self.domain.subPolygon)):
            x_sub, y_sub = zip(*self.domain.subPolygon[i].xy)
            line_sub = Line2D(x_sub, y_sub, linestyle= '-')
            line_sub.set_color('k')
            ax.add_line(line_sub)
            self.domain.subPolygon[i].set_color(color[7])
            ax.add_patch(self.domain.subPolygon[i]) 

        ### Initialise highlighted edges as hidden lines ###
        plotted_edges = []
        for edge in self.domain.edges:
            edge = np.array(edge)
            plotted_edges.append(ax.plot(edge[:,0], edge[:,1], color='r', visible=False))
        plotted_subedges = []
        for subedge in self.domain.subedges:
            subedge = np.array(subedge)
            plotted_subedges.append(ax.plot(subedge[:,0], subedge[:,1], color='r', visible=False))

        ### Computational domain plotting ### 
        ax.set_title(r'Computational Domain')
        ax.set_xlabel(r'$x$')
        ax.set_ylabel(r'$y$')
        ax.set_xlim((min(x)-0.1*max(x), 1.1*max(x)))
        ax.set_ylim((min(y)-0.1*max(y), 1.1*max(y)))
        ax.grid(True)
        ax.set_aspect(aspect= 1)
        
        pts = pylab.scatter(self.nodes.coor[:,0], 
                    self.nodes.coor[:,1],
                    s = 4,
                    color = 'k')

        spheres=[]
        self.size = abs(self.nodes.coor[0,1] - self.nodes.coor[1,1])
        for i in range(0, len(self.nodes.coor)):
            spheres.append(pylab.Circle((self.nodes.coor[i,0],self.nodes.coor[i,1]), self.size, fill=False, color=(0.6627,0.6627,0.6627)))
            ax.add_artist(spheres[i])

        ### Maximise window ###
        ### with Tk
        # mng = plt.get_current_fig_manager()
        # mng.resize(*mng.window.maxsize())
        ### with Qt
        figManager = pylab.get_current_fig_manager()
        figManager.window.showMaximized()
        
        ### Button formatting ###
        EDGESax = pylab.axes([0.05, 0.58, 0.2, 0.3])
        
        BCax = pylab.axes([0.05, 0.45, 0.2, 0.1])
        BCax_inner1 = pylab.axes([0.15, 0.46, 0.04, 0.03])
        BCax_inner2 = pylab.axes([0.15, 0.51, 0.04, 0.03])
        BCax_inner3 = pylab.axes([0.2, 0.48, 0.04, 0.04])
        
        BCBox1 = TextBox(BCax_inner1, r'$u^s_y(x,y)=$ ')
        BCBox1.set_val("0")
        BCBox2 = TextBox(BCax_inner2, r'$u^s_x(x,y)=$ ')
        BCBox2.set_val("0")

        BCax.text(0.03,0.45, 'Prescribed \n displacements:') 
        BCax.get_yaxis().set_visible(False)
        BCax.get_xaxis().set_visible(False)
        
        LOADax = pylab.axes([0.05,0.3,0.2,0.1])
        LOADax_inner1 = pylab.axes([0.15, 0.31, 0.04, 0.03])
        LOADax_inner2 = pylab.axes([0.15, 0.36, 0.04, 0.03])
        LOADax_inner3 = pylab.axes([0.2, 0.33, 0.04, 0.04])
        
        LOADax.get_yaxis().set_visible(False)
        LOADax.get_xaxis().set_visible(False)
        LOADax.text(0.03,0.45, 'Traction \n vector:')
        LOADBox1 = TextBox(LOADax_inner1, r'$f^s_y(x,y)=$ ')
        LOADBox1.set_val("0")
        LOADBox2 = TextBox(LOADax_inner2, r'$f^s_x(x,y)=$ ')
        LOADBox2.set_val("0")

        MATax = pylab.axes([0.05,0.12,0.2,0.15])
        MATax.get_yaxis().set_visible(False)
        MATax.get_xaxis().set_visible(False)
        MATax.text(0.03,0.45, 'Physical \n properties:')
        MATax_inner1 = pylab.axes([0.15, 0.13, 0.04, 0.03])
        MATax_inner2 = pylab.axes([0.15, 0.18, 0.04, 0.03])
        MATax_inner3 = pylab.axes([0.15, 0.23, 0.04, 0.03])
        MATBox3 = TextBox(MATax_inner3, r'$E=$ ')
        MATBox3.set_val("100")
        MATBox2 = TextBox(MATax_inner2, r'$\nu=$ ')
        MATBox2.set_val("0.3")
        MATBox1 = TextBox(MATax_inner1, r'$t=$ ')
        MATBox1.set_val("1")

        SUBMITax = pylab.axes([0.76, 0.45, 0.15, 0.1])
        
        ### Generate check buttons ###
        EDGE_labels = []
        EDGE_flags = []
        for i in range(1,len(self.domain.edges)+1):
            EDGE_labels.append(str('Edge %3.0f' % i))
            EDGE_flags.append(False)
        for i in range(1,len(self.domain.subedges)+1):
            EDGE_labels.append(str('Annulus %3.0f' % i))
            EDGE_flags.append(False)

        EDGEcheck = CheckButtons(EDGESax, EDGE_labels, EDGE_flags)
        
        ### Toggle edge visibility when toggling check buttons ###
        def EDGEcheck_click(label):
            if 'Edge' in label:
                edge_number = [int(s) for s in label.split() if s.isdigit()]
                edge_number = int(edge_number[0]) - 1
                plotted_edges[edge_number][0].set_visible(not plotted_edges[edge_number][0].get_visible())
                pylab.draw()
            elif 'Annulus' in label:
                subedge_number = [int(s) for s in label.split() if s.isdigit()]
                subedge_number = int(subedge_number[0]) - 1
                plotted_subedges[subedge_number][0].set_visible(not plotted_subedges[subedge_number][0].get_visible())
                pylab.draw()

        
        EDGEcheck.on_clicked(EDGEcheck_click)

        ### Generate buttons ###
        BCbutton = Button(BCax_inner3, 'Apply')
        LOADbutton = Button(LOADax_inner3, 'Apply')
        SUBMITbutton = Button(SUBMITax, 'Submit job')

        ### Interactivity ###
        self.u_sx = []
        self.u_sx_flag = []
        self.u_sy = []
        self.u_sy_flag = []
        self.u_s_edges = []

        def applyBC(expression):
            if BCBox2.text:
                self.u_sx.append(float(BCBox2.text))
                self.u_sx_flag.append(0)
            else:
                self.u_sx.append(0)
                self.u_sx_flag.append(1)

            if BCBox1.text:
                self.u_sy.append(float(BCBox1.text))
                self.u_sy_flag.append(0)
            else:
                self.u_sy.append(0)
                self.u_sy_flag.append(1)           

            
            active_edges = []
            status = EDGEcheck.get_status()
            for i in range(len(status)):
                if status[i]==True:
                    active_edges.append(i)
            self.u_s_edges.append(active_edges)

            return self.u_sx, self.u_sy, self.u_sx_flag, self.u_sy_flag, self.u_s_edges
        
        BCbutton.on_clicked(applyBC)

        self.f_sx = []
        self.f_sx_flag = []
        self.f_sy = []
        self.f_sy_flag = []
        self.f_s_edges = []
        def applyLOAD(expression):
            if LOADBox2.text:
                self.f_sx.append(float(LOADBox2.text))
                self.f_sx_flag.append(0)
            else:
                self.f_sx.append(0)
                self.f_sx_flag.append(1)

            if LOADBox1.text:
                self.f_sy.append(float(LOADBox1.text))
                self.f_sy_flag.append(0)
            else:
                self.f_sy.append(0)
                self.f_sy_flag.append(1)   
            
            active_edges = []
            status = EDGEcheck.get_status()
            for i in range(len(status)):
                if status[i]==True:
                    active_edges.append(i)
            self.f_s_edges.append(active_edges)

            return self.f_sx, self.f_sy, self.f_s_edges

        LOADbutton.on_clicked(applyLOAD)


        def submitJob(expression):
            pylab.close()
        SUBMITbutton.on_clicked(submitJob)

        pylab.show()

        self.u_s = [self.u_sx, self.u_sy, self.u_sx_flag, self.u_sy_flag, self.u_s_edges]
        self.f_s = [self.f_sx, self.f_sy, self.f_sx_flag, self.f_sy_flag, self.f_s_edges]
        self.E = MATBox3.text
        self.nu = MATBox2.text
        self.t = MATBox1.text
        self.PhysicalProperties = [float(self.E),float(self.nu),float(self.t)]

        return self.u_s, self.f_s, self.PhysicalProperties

    '''
    Compile the .mfs input file for given pre-processing inputs
    '''
    def compile_MFS(self):
        lines = [
                 f"# Job: {self.job_ID}.mfs",\
                 '',\
                 "# External surfaces", \
                 str(self.domain.edges), \
                 '',\
                 "# Internal surfaces", \
                 str(self.domain.subedges), \
                 '',\
                 "# Nodal coordinates", \
                 str(self.nodes.coor.tolist()), \
                 '',\
                 "# Sphere radius", \
                 str(self.size), \
                 '',\
                 "# Physical properties", \
                 str(self.PhysicalProperties), \
                 '',\
                 "# Prescribed displacements [[u_sx], [u_sy], [u_sx_flag], [u_sy_flag], [surfaces]]", \
                 str(self.u_s), \
                 '', \
                 "# Prescribed loads [[f_sx], [f_sy], [f_sx_flag], [f_sy_flag], [surfaces]]", \
                 str(self.f_s)]

        job = f"{self.job_ID}.mfs"
        with open(job, 'w') as f:
            for line in lines:
                f.write(line)
                f.write('\n')

        return 


        
        

    


        


        






   