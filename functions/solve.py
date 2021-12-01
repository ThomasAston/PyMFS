# Let's solve the system of matrices for unknown forces and displacements

# To do list:
# 
#
#

import numpy as np

class solve:
    def __init__(self, K_mat, DoFs):
        self.K_mat = K_mat
        self.DoFs = DoFs

        self.u, self.F = self.solve()
    

    def solve(self):
        DoF_num = self.DoFs.DoF_num
        DoC_num = self.DoFs.DoC_num

        # Slice the stiffness matrix for free and restricted degrees of freedom.
        K_ff = self.K_mat.value[0:DoF_num, 0:DoF_num]
        K_fr = self.K_mat.value[0:DoF_num, DoF_num:DoC_num+DoF_num]
        K_rf = self.K_mat.value[DoF_num:DoC_num+DoF_num, 0:DoF_num]
        K_rr = self.K_mat.value[DoF_num:DoC_num+DoF_num, DoF_num:DoC_num+DoF_num]

        # Slice force vector for the known forces
        F_f = self.DoFs.f[0:DoF_num]

        # Read the known displacements
        q_r = self.DoFs.q[DoF_num:DoC_num+DoF_num]

        # print(q_r)
        # Calculate the unknowns    
        F = F_f - (K_fr @ q_r)
        q_f = np.linalg.solve(K_ff, F)
        F_r = (K_rf @ q_f) + (K_rr @ q_r)

        # Reassemble values in nodal order
        q = np.zeros((len(self.DoFs.order), 1))
        q[0:DoF_num] = q_f
        q[DoF_num:DoF_num+DoC_num] = q_r
        
        q_ordered = np.zeros((len(q),1))
        for i in range(0, len(self.DoFs.order)):
            # print(self.DoFs.order[i])
            q_ordered[i] = q[int(self.DoFs.order[i])]

        u = np.zeros((int(len(q_ordered)/2),2))
        for i in range(0, len(q_ordered), 2):
            u[int(i/2),0] = q_ordered[i]
            u[int(i/2),1] = q_ordered[i+1]
        
    
        F = np.zeros((len(self.DoFs.order), 1))
        F[0:DoF_num] = F_f
        F[DoF_num:DoF_num+DoC_num] = F_r
        
        F_ordered = np.zeros((len(q),1))
        for i in range(0, len(self.DoFs.order)):
            # print(self.DoFs.order[i])
            F_ordered[i] = F[int(self.DoFs.order[i])]

        F = np.zeros((int(len(F_ordered)/2),2))
        for i in range(0, len(F_ordered), 2):
            F[int(i/2),0] = F_ordered[i]
            F[int(i/2),1] = F_ordered[i+1]

        return u, F