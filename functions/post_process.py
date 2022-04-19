'''
PyMFS post-processing module

Takes a solution object past into it from the solver and reassembles
function fields before plotting chosen results.
'''
from .B_mat import *
from .C_mat import *
from .shape_functions import *
from matplotlib import pyplot as plt
from matplotlib import ticker
import csv
from matplotlib.patches import Polygon
from matplotlib.cm import ScalarMappable
from matplotlib.lines import Line2D
from copy import copy

class post_process:
    def __init__(self, solution):
        self.solution = solution
        self.reassemble()

    '''
    Function for reassembling field quantities from obtained DoFs, and dumping into
    output file. Results are then also plotted.
    '''
    def reassemble(self):
        
        '''
        Initialise output file, write first row
        '''
        # output_file = (self.solution.job_ID[:-3],'csv')
        output_file = 'examples\\2d_elastostatics\\hole\\hole_17x33.csv'
        header = ['x', 'y', 'u_x', 'u_y', 'strain_xx', 'strain_yy', 'strain_xy', 'stress_xx','stress_yy', 'stress_xy']

        '''External polygon'''
        vertices=[]
        for surf in self.solution.input_data.external_surfaces:
            vertices += surf
        self.polygon = Polygon(vertices)
        # self.polygon.set_color('none')
        self.polygon.set_alpha(0)
        self.polygon.set_color((0,0.541176,0.541176))
        self.path = self.polygon.get_path()

        '''Internal polygons'''
        self.subPolygon=[]
        self.subPath=[]
        if self.solution.input_data.internal_surfaces:
            for i in range(len(self.solution.input_data.internal_surfaces)):
                self.subPolygon.append(Polygon(self.solution.input_data.internal_surfaces[i]))
                self.subPath.append(self.subPolygon[i].get_path())

        '''Domain bounds'''
        bounds = np.array([np.min(vertices,axis=0), np.max(vertices,axis=0)])

        Ntot = len(self.solution.input_data.node_coor)
        order = 4
        dimensions = 2
        
        nx = 11
        ny = 11
        X = np.linspace(bounds[0,0]+0.0000001, bounds[1,0]-0.0000001, nx)
        Y = np.linspace(bounds[0,1]+0.0000001, bounds[1,1]-0.0000001, ny)
        plotX, plotY = np.meshgrid(X,Y)
        
        Ntot = len(self.solution.input_data.node_coor)
        dimensions = 2
        u = np.zeros([Ntot,2])

        q = np.reshape(self.solution.q, [Ntot, order, 2])
        # q = np.reshape(self.solution.q, [Ntot, 2, order])
        for J in range(Ntot):
            for n in range(order):
                h_Jn = shape_functions(J,n,self.solution.input_data.node_coor[J],self.solution.input_data)
                H_Jn = np.array([[h_Jn, 0],\
                                    [0, h_Jn]])
                u[J, :] += H_Jn@q[J,n]


        deformed_shape = self.solution.input_data.node_coor + u*2
        fig,ax = plt.subplots()
        fig.canvas.manager.set_window_title('Deformed shape') 
        # cp = plt.contourf(plotX, plotY, strain[:,:,0])
        # plt.colorbar(cp)
        # ax.scatter(np.array(self.solution.input_data.node_coor)[:,0], np.array(self.solution.input_data.node_coor)[:,1], color='k')
        ax.scatter(deformed_shape[:,0], deformed_shape[:,1], color = [(0,0.541176,0.541176)])
        ax.set_aspect(aspect= 1)
        plt.show()


        u = np.zeros([nx,ny,2])
        strain=np.zeros([nx,ny,3])
        stress=np.zeros([nx,ny,3]) 
        with open(output_file,'w',newline='') as f:
            writer = csv.writer(f)
            writer.writerow(header)
            for a in range(nx):
                for b in range(ny):
                    print('Reassembling for x: ', a,', y: ', b)
                    point_x = X[a]
                    point_y = Y[b]
                    point = (point_x,point_y)
                    pointInDomain = self.polygon.contains_point(point)
                    if pointInDomain==True:
                        pointInSubDomain = []
                        for poly in range(len(self.subPolygon)):
                            if self.subPolygon[poly].contains_point(point)==True:
                                pointInSubDomain.append(1)
                        check = (1 not in pointInSubDomain)
                        if check==True:
                            for J in range(Ntot):
                                for n in range(order):
                                    # if J*order*dimensions not in self.solution.delete[:] and J*order*dimensions+1 not in self.solution.delete[:]:
                                    B_Jn = B(J,n,point,self.solution.input_data).mat
                                    h_Jn = shape_functions(J,n,point,self.solution.input_data)
                                    H_Jn = np.array([[h_Jn, 0],\
                                                    [0, h_Jn]])
                                    u[b,a,:] += H_Jn @ q[J,n]
                                    strain[b,a,:] += B_Jn @ q[J,n]
                                    stress[b,a,:] += C(self.solution.input_data) @ B_Jn @ q[J,n]

                
                    writer.writerow([point_x,point_y,u[b,a,0],u[b,a,1],strain[b,a,0],strain[b,a,1],strain[b,a,2],stress[b,a,0],stress[b,a,1],stress[b,a,2]])



        plt.rc('text', usetex=True)
        plt.rc('text.latex', preamble=r'\usepackage{cmbright}')
        plt.rcParams.update({'font.size': 18})


#########################################################################################################
        '''
        PLOTTING u_x
        '''
        # vmin=0
        # vmax=0.25
        fig,ax = plt.subplots()
        fig.canvas.manager.set_window_title('Ux') 
                
        '''External polygon'''
        vertices=[]
        for surf in self.solution.input_data.external_surfaces:
            vertices += surf
        self.polygon = Polygon(vertices)
        ax.add_patch(self.polygon)
        self.path = self.polygon.get_path()
        x, y = zip(*self.polygon.xy)
        line = Line2D(x, y, linestyle= '-')
        ax.add_line(line)
        line.set_color('k')

        cp = ax.contourf(plotX, plotY, u[:,:,0],levels=400)# levels=np.linspace(vmin,vmax,400))
        clb = fig.colorbar(ScalarMappable(norm=cp.norm, cmap=cp.cmap))
        clb.ax.set_title(r'$u_x~\mathrm{[m]}$')
        tick_locator = ticker.MaxNLocator(nbins=5)
        clb.locator = tick_locator
        clb.update_ticks()
        ax.set_xlabel(r'$x~\mathrm{[m]}$')
        ax.set_ylabel(r'$y~\mathrm{[m]}$')
        ax.scatter(np.array(self.solution.input_data.node_coor)[:,0], np.array(self.solution.input_data.node_coor)[:,1], color='k',s = 4)
        ax.set_aspect(aspect= 1)
        plt.xlim([bounds[0,0]-0.1, bounds[1,0]+0.1])
        plt.ylim([bounds[0,1]-0.1, bounds[1,1]+0.1])
        
        for col in cp.collections:
            col.set_clip_path(self.polygon)

        for i in range(len(self.subPolygon)):
            x_sub, y_sub = zip(*self.subPolygon[i].xy)
            line_sub = Line2D(x_sub, y_sub, linestyle= '-')
            line_sub.set_color('k')
            ax.add_line(line_sub)
            self.subPolygon[i].set_color('w')
            subPoly = copy(self.subPolygon[i])
            ax.add_patch(subPoly)
        
        plt.tight_layout()
        plt.show()


##############################################################################################################

        '''
        PLOTTING u_y
        '''
        # vmin=-0.05
        # vmax=0.2
        fig,ax = plt.subplots()
        fig.canvas.manager.set_window_title('Uy') 
                
        '''External polygon'''
        vertices=[]
        for surf in self.solution.input_data.external_surfaces:
            vertices += surf
        self.polygon = Polygon(vertices)
        ax.add_patch(self.polygon)
        self.path = self.polygon.get_path()
        x, y = zip(*self.polygon.xy)
        line = Line2D(x, y, linestyle= '-')
        ax.add_line(line)
        line.set_color('k')

        cp = ax.contourf(plotX, plotY, u[:,:,1],levels=100)
        clb = fig.colorbar(ScalarMappable(norm=cp.norm, cmap=cp.cmap))
        clb.ax.set_title(r'$u_y~\mathrm{[m]}$')
        tick_locator = ticker.MaxNLocator(nbins=5)
        clb.locator = tick_locator
        clb.update_ticks()
        ax.set_xlabel(r'$x~\mathrm{[m]}$')
        ax.set_ylabel(r'$y~\mathrm{[m]}$')
        ax.scatter(np.array(self.solution.input_data.node_coor)[:,0], np.array(self.solution.input_data.node_coor)[:,1], color='k',s = 4)
        ax.set_aspect(aspect= 1)
        plt.xlim([bounds[0,0]-0.1, bounds[1,0]+0.1])
        plt.ylim([bounds[0,1]-0.1, bounds[1,1]+0.1])
        
        for col in cp.collections:
            col.set_clip_path(self.polygon)

        for i in range(len(self.subPolygon)):
            x_sub, y_sub = zip(*self.subPolygon[i].xy)
            line_sub = Line2D(x_sub, y_sub, linestyle= '-')
            line_sub.set_color('k')
            ax.add_line(line_sub)
            self.subPolygon[i].set_color('w')
            subPoly = copy(self.subPolygon[i])
            ax.add_patch(subPoly)
        
        plt.tight_layout()
        plt.show()

##############################################################################################################

        '''
        PLOTTING Strain_xx
        '''
        # vmin=-0.05
        # vmax=0.2
        fig,ax = plt.subplots()
        fig.canvas.manager.set_window_title('Strain_xx') 
                
        '''External polygon'''
        vertices=[]
        for surf in self.solution.input_data.external_surfaces:
            vertices += surf
        self.polygon = Polygon(vertices)
        ax.add_patch(self.polygon)
        self.path = self.polygon.get_path()
        x, y = zip(*self.polygon.xy)
        line = Line2D(x, y, linestyle= '-')
        ax.add_line(line)
        line.set_color('k')

        cp = plt.contourf(plotX, plotY, strain[:,:,0],levels=100)
        clb = plt.colorbar(cp)
        clb.ax.set_title(r'$\epsilon_{xx}~\mathrm{[-]}$')
        tick_locator = ticker.MaxNLocator(nbins=5)
        clb.locator = tick_locator
        clb.update_ticks()
        ax.set_xlabel(r'$x~\mathrm{[m]}$')
        ax.set_ylabel(r'$y~\mathrm{[m]}$')
        ax.scatter(np.array(self.solution.input_data.node_coor)[:,0], np.array(self.solution.input_data.node_coor)[:,1], color='k',s = 4)
        ax.set_aspect(aspect= 1)
        plt.xlim([bounds[0,0]-0.1, bounds[1,0]+0.1])
        plt.ylim([bounds[0,1]-0.1, bounds[1,1]+0.1])
        
        for col in cp.collections:
            col.set_clip_path(self.polygon)

        for i in range(len(self.subPolygon)):
            x_sub, y_sub = zip(*self.subPolygon[i].xy)
            line_sub = Line2D(x_sub, y_sub, linestyle= '-')
            line_sub.set_color('k')
            ax.add_line(line_sub)
            self.subPolygon[i].set_color('w')
            ax.add_patch(copy(self.subPolygon[i]))
        
        plt.tight_layout()
        plt.show()
        

        ##############################################################################################################

        '''
        PLOTTING Strain_yy
        '''
        # vmin=-0.05
        # vmax=0.2
        fig,ax = plt.subplots()
        fig.canvas.manager.set_window_title('Strain_yy') 
                
        '''External polygon'''
        vertices=[]
        for surf in self.solution.input_data.external_surfaces:
            vertices += surf
        self.polygon = Polygon(vertices)
        ax.add_patch(self.polygon)
        self.path = self.polygon.get_path()
        x, y = zip(*self.polygon.xy)
        line = Line2D(x, y, linestyle= '-')
        ax.add_line(line)
        line.set_color('k')

        cp = plt.contourf(plotX, plotY, strain[:,:,1],levels=100)
        clb = plt.colorbar(cp)
        clb.ax.set_title(r'$\epsilon_{yy}~\mathrm{[-]}$')
        tick_locator = ticker.MaxNLocator(nbins=5)
        clb.locator = tick_locator
        clb.update_ticks()
        ax.set_xlabel(r'$x~\mathrm{[m]}$')
        ax.set_ylabel(r'$y~\mathrm{[m]}$')
        ax.scatter(np.array(self.solution.input_data.node_coor)[:,0], np.array(self.solution.input_data.node_coor)[:,1], color='k',s = 4)
        ax.set_aspect(aspect= 1)
        plt.xlim([bounds[0,0]-0.1, bounds[1,0]+0.1])
        plt.ylim([bounds[0,1]-0.1, bounds[1,1]+0.1])
        
        for col in cp.collections:
            col.set_clip_path(self.polygon)

        for i in range(len(self.subPolygon)):
            x_sub, y_sub = zip(*self.subPolygon[i].xy)
            line_sub = Line2D(x_sub, y_sub, linestyle= '-')
            line_sub.set_color('k')
            ax.add_line(line_sub)
            self.subPolygon[i].set_color('w')
            ax.add_patch(copy(self.subPolygon[i]))
        
        plt.tight_layout()
        plt.show()

        ##############################################################################################################

        '''
        PLOTTING Strain_xy
        '''
        # vmin=-0.05
        # vmax=0.2
        fig,ax = plt.subplots()
        fig.canvas.manager.set_window_title('Strain_xy') 
                
        '''External polygon'''
        vertices=[]
        for surf in self.solution.input_data.external_surfaces:
            vertices += surf
        self.polygon = Polygon(vertices)
        ax.add_patch(self.polygon)
        self.path = self.polygon.get_path()
        x, y = zip(*self.polygon.xy)
        line = Line2D(x, y, linestyle= '-')
        ax.add_line(line)
        line.set_color('k')

        cp = plt.contourf(plotX, plotY, strain[:,:,2],levels=100)
        clb = plt.colorbar(cp)
        clb.ax.set_title(r'$\epsilon_{xy}~\mathrm{[-]}$')
        tick_locator = ticker.MaxNLocator(nbins=5)
        clb.locator = tick_locator
        clb.update_ticks()
        ax.set_xlabel(r'$x~\mathrm{[m]}$')
        ax.set_ylabel(r'$y~\mathrm{[m]}$')
        ax.scatter(np.array(self.solution.input_data.node_coor)[:,0], np.array(self.solution.input_data.node_coor)[:,1], color='k',s = 4)
        ax.set_aspect(aspect= 1)
        plt.xlim([bounds[0,0]-0.1, bounds[1,0]+0.1])
        plt.ylim([bounds[0,1]-0.1, bounds[1,1]+0.1])
        
        for col in cp.collections:
            col.set_clip_path(self.polygon)

        for i in range(len(self.subPolygon)):
            x_sub, y_sub = zip(*self.subPolygon[i].xy)
            line_sub = Line2D(x_sub, y_sub, linestyle= '-')
            line_sub.set_color('k')
            ax.add_line(line_sub)
            self.subPolygon[i].set_color('w')
            ax.add_patch(copy(self.subPolygon[i]))
        
        plt.tight_layout()
        plt.show()
        
        
        ##############################################################################################################

        '''
        PLOTTING Stress_xx
        '''
        # vmin=-0.05
        # vmax=0.2
        fig,ax = plt.subplots()
        fig.canvas.manager.set_window_title('Stress_xx') 
                
        '''External polygon'''
        vertices=[]
        for surf in self.solution.input_data.external_surfaces:
            vertices += surf
        self.polygon = Polygon(vertices)
        ax.add_patch(self.polygon)
        self.path = self.polygon.get_path()
        x, y = zip(*self.polygon.xy)
        line = Line2D(x, y, linestyle= '-')
        ax.add_line(line)
        line.set_color('k')

        cp = plt.contourf(plotX, plotY, stress[:,:,0],levels=100)
        clb = plt.colorbar(cp)
        clb.ax.set_title(r'$\sigma_{xx}~\mathrm{[-]}$')
        tick_locator = ticker.MaxNLocator(nbins=5)
        clb.locator = tick_locator
        clb.update_ticks()
        ax.set_xlabel(r'$x~\mathrm{[m]}$')
        ax.set_ylabel(r'$y~\mathrm{[m]}$')
        ax.scatter(np.array(self.solution.input_data.node_coor)[:,0], np.array(self.solution.input_data.node_coor)[:,1], color='k',s = 4)
        ax.set_aspect(aspect= 1)
        plt.xlim([bounds[0,0]-0.1, bounds[1,0]+0.1])
        plt.ylim([bounds[0,1]-0.1, bounds[1,1]+0.1])
        
        for col in cp.collections:
            col.set_clip_path(self.polygon)

        for i in range(len(self.subPolygon)):
            x_sub, y_sub = zip(*self.subPolygon[i].xy)
            line_sub = Line2D(x_sub, y_sub, linestyle= '-')
            line_sub.set_color('k')
            ax.add_line(line_sub)
            self.subPolygon[i].set_color('w')
            ax.add_patch(copy(self.subPolygon[i]))
        
        plt.tight_layout()
        plt.show()

 ##############################################################################################################

        '''
        PLOTTING Stress_yy
        '''
        # vmin=-0.05
        # vmax=0.2
        fig,ax = plt.subplots()
        fig.canvas.manager.set_window_title('Stress_yy') 
                
        '''External polygon'''
        vertices=[]
        for surf in self.solution.input_data.external_surfaces:
            vertices += surf
        self.polygon = Polygon(vertices)
        ax.add_patch(self.polygon)
        self.path = self.polygon.get_path()
        x, y = zip(*self.polygon.xy)
        line = Line2D(x, y, linestyle= '-')
        ax.add_line(line)
        line.set_color('k')

        cp = plt.contourf(plotX, plotY, stress[:,:,1],levels=100)
        clb = plt.colorbar(cp)
        clb.ax.set_title(r'$\sigma_{yy}~\mathrm{[-]}$')
        tick_locator = ticker.MaxNLocator(nbins=5)
        clb.locator = tick_locator
        clb.update_ticks()
        ax.set_xlabel(r'$x~\mathrm{[m]}$')
        ax.set_ylabel(r'$y~\mathrm{[m]}$')
        ax.scatter(np.array(self.solution.input_data.node_coor)[:,0], np.array(self.solution.input_data.node_coor)[:,1], color='k',s = 4)
        ax.set_aspect(aspect= 1)
        plt.xlim([bounds[0,0]-0.1, bounds[1,0]+0.1])
        plt.ylim([bounds[0,1]-0.1, bounds[1,1]+0.1])
        
        for col in cp.collections:
            col.set_clip_path(self.polygon)

        for i in range(len(self.subPolygon)):
            x_sub, y_sub = zip(*self.subPolygon[i].xy)
            line_sub = Line2D(x_sub, y_sub, linestyle= '-')
            line_sub.set_color('k')
            ax.add_line(line_sub)
            self.subPolygon[i].set_color('w')
            ax.add_patch(copy(self.subPolygon[i]))
        
        plt.tight_layout()
        plt.show()

        return u
        
        
                