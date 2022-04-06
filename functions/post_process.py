'''
PyMFS post-processing module

Takes a solution object past into it from the solver and reassembles
function fields before plotting chosen results.
'''
from .B_mat import *
from .C_mat import *
from .shape_functions import *
import matplotlib
from matplotlib import pyplot as plt

class post_process:
    def __init__(self, solution):
        self.solution = solution
        self.u, self.u_x, self.u_y = self.deformed_shape()
        # self.strain = self.strain()
        # print(self.u)

    '''
    Function for re-assembling the nodal displacements from the obtained
    solution degrees of freedom
    '''
    def deformed_shape(self):
        Ntot = len(self.solution.input_data.node_coor)
        order = 4
        dimensions = 2
        u= np.zeros([Ntot,2])
        q = np.reshape(self.solution.q, [Ntot, order, 2])
        for J in range(Ntot):
            for n in range(order):
                for dim in range(dimensions):
                    h_Jn = shape_functions(J,n,self.solution.input_data.node_coor[J],self.solution.input_data)
                    H_Jn = np.array([[h_Jn, 0],\
                                     [0, h_Jn]])
                    u[J, dim] += h_Jn*q[J,n, dim]
        
        deformed_shape = self.solution.input_data.node_coor + u*0.3
        fig,ax = plt.subplots()
        fig.canvas.manager.set_window_title('Deformed shape') 
        # cp = plt.contourf(plotX, plotY, strain[:,:,0])
        # plt.colorbar(cp)
        ax.scatter(np.array(self.solution.input_data.node_coor)[:,0], np.array(self.solution.input_data.node_coor)[:,1], color='k')
        ax.scatter(deformed_shape[:,0], deformed_shape[:,1], color = [(0,0.541176,0.541176)])
        ax.set_aspect(aspect= 1)
        plt.show()

        u_x = np.zeros([30,30])
        u_y = np.zeros([30,30])
        X = np.linspace(0,2,30)
        Y = np.linspace(0,2,30)
        plotX, plotY = np.meshgrid(X,Y)
        for a in range(30):
            for b in range(30):
                point_x = X[a]
                point_y = Y[b]
                point = (point_x,point_y)
                for J in range(Ntot):
                    for n in range(order):
                        for dim in range(dimensions):  
                            h_Jn = shape_functions(J,n,point,self.solution.input_data)
                            H_Jn = np.array([[h_Jn, 0],\
                                            [0, h_Jn]])
                            u[J, dim] += h_Jn*q[J,n, dim]

                            if dim == 0:
                                u_x[b,a] += h_Jn*q[J,n, dim]
                            elif dim==1:
                                u_y[b,a] += h_Jn*q[J,n, dim]

        # deformed_shape = self.solution.input_data.node_coor + u*0.3
        fig,ax = plt.subplots()
        fig.canvas.manager.set_window_title('Ux') 
        cp = plt.contourf(plotX, plotY, u_x, levels=100)
        plt.colorbar(cp)
        ax.scatter(np.array(self.solution.input_data.node_coor)[:,0], np.array(self.solution.input_data.node_coor)[:,1], color='k')
        # ax.scatter(deformed_shape[:,0], deformed_shape[:,1], color = [(0,0.541176,0.541176)])
        ax.set_aspect(aspect= 1)
        plt.show()

        fig,ax = plt.subplots()
        fig.canvas.manager.set_window_title('Uy') 
        cp = plt.contourf(plotX, plotY, u_y, levels=100)
        plt.colorbar(cp)
        ax.scatter(np.array(self.solution.input_data.node_coor)[:,0], np.array(self.solution.input_data.node_coor)[:,1], color='k')
        # ax.scatter(deformed_shape[:,0], deformed_shape[:,1], color = [(0,0.541176,0.541176)])
        ax.set_aspect(aspect= 1)
        plt.show()

        return u, u_x, u_y

    '''
    Function for calculating strain field
    '''
    def strain(self):
        Ntot = len(self.solution.input_data.node_coor)
        order = 4
        dimensions = 2
        q = np.reshape(self.solution.q, [Ntot, order, 2])
        X = np.linspace(0,2,30)
        Y = np.linspace(0,2,30)
        plotX, plotY = np.meshgrid(X,Y)
        
        strain=np.zeros([30,30,3])
        strain_xx = np.zeros([30,30])
        strain_yy = np.zeros([30,30])
        strain_xy = np.zeros([30,30])
        
        stress=np.zeros([30,30,3]) 
        stress_xx = np.zeros([30,30])
        stress_yy = np.zeros([30,30])
        stress_xy = np.zeros([30,30])
        for a in range(30):
            for b in range(30):
                point_x = X[a]
                point_y = Y[b]
                point = (point_x,point_y)
                for J in range(Ntot):
                    for n in range(order):
                        # for dim in range(3):
                        B_Jn = B(J,n,point,self.solution.input_data).mat
                        u = np.array([self.u_x[b,a],self.u_y[b,a]])
                        dh_dx = B_Jn[0,0]
                        dh_dy = B_Jn[1,1]
                        strain[b,a,:] += B_Jn @ q[J,n]
                        stress[b,a,:] += C(self.solution.input_data) @ B_Jn @ q[J,n]

                        for dim in range(dimensions):
                            if dim==0:
                                strain_xx[b,a] += dh_dx*q[J,n, dim]
                            if dim==1:
                                strain_yy[b,a] += dh_dy*q[J,n, dim]
                # strain_xx[b,a] = strain[b,a,0]
                # strain_yy[b,a] = strain[b,a,1]
                strain_xy[b,a] = strain[b,a,2]
       
        fig,ax = plt.subplots()
        fig.canvas.manager.set_window_title('Strain_xx') 
        cp = plt.contourf(plotX, plotY, strain_xx,levels=100)
        plt.colorbar(cp)
        ax.scatter(np.array(self.solution.input_data.node_coor)[:,0], np.array(self.solution.input_data.node_coor)[:,1], color='k')
        ax.set_aspect(aspect= 1)
        plt.show()

        fig,ax = plt.subplots()
        fig.canvas.manager.set_window_title('Strain_yy') 
        cp = plt.contourf(plotX, plotY, strain_yy,levels=100)
        plt.colorbar(cp)
        ax.scatter(np.array(self.solution.input_data.node_coor)[:,0], np.array(self.solution.input_data.node_coor)[:,1], color='k')
        ax.set_aspect(aspect= 1)
        plt.show()

        fig,ax = plt.subplots()
        fig.canvas.manager.set_window_title('Strain_xy') 
        cp = plt.contourf(plotX, plotY, strain_xy,levels=100)
        plt.colorbar(cp)
        ax.scatter(np.array(self.solution.input_data.node_coor)[:,0], np.array(self.solution.input_data.node_coor)[:,1], color='k')
        ax.set_aspect(aspect= 1)
        plt.show()

        
        fig,ax = plt.subplots()
        fig.canvas.manager.set_window_title('Stress_xx') 
        cp = plt.contourf(plotX, plotY, stress[:,:,0],levels=100)
        plt.colorbar(cp)
        ax.scatter(np.array(self.solution.input_data.node_coor)[:,0], np.array(self.solution.input_data.node_coor)[:,1], color='k')
        ax.set_aspect(aspect= 1)
        plt.show()

        
        return strain
        
        
                