# -*- coding: utf-8 -*-
"""
Created on Sat Oct  1 18:33:11 2022

@author: Yuanshun Xie
"""

from functions import *
import numpy as np
import time
import csv
import matplotlib.pyplot as plt

class kinematic:
    def __init__(self,robot,T_total,q_e):
        self.q_a=np.zeros((robot.dimension,1))
        self.q_e=q_e*np.ones((robot.dimension,1))
        self.T_total=T_total
        self.robot=robot
        
    def calculate_jointwinkel(self):
        T_joint = self.T_total/self.robot.dimension
        N_joint = int(np.floor(T_joint/self.robot.dt))
        T=self.robot.dt*np.array(list(range(self.robot.dimension*N_joint)))
        
        Q=np.zeros((self.robot.dimension,len(T)))
        W=np.zeros((3,len(T)))
        dot_W= np.zeros((3,len(T)))
        V=np.zeros((3,4,self.robot.dimension,len(T)))
    
        a=time.time()  
        for i in range(self.robot.dimension):
            delta_q=(self.q_e[i]-self.q_a[i])/N_joint
            
            self.robot.dot_q=np.zeros((self.robot.dimension,1))
            self.robot.dot_q[i]= (self.q_e[i]-self.q_a[i])/T_joint
            
            for k in range(N_joint):
                self.robot.q[i]+=delta_q
                self.robot.time+=self.robot.dt
                
                self.dk_position_vectorchain()
                
                self.calculate_dk_velocity()
            
                Q[:,N_joint*i+k] = self.robot.q
                W[:,N_joint*i+k] = np.transpose(self.robot.w)[0]
                dot_W[:,N_joint*i+k] = np.transpose(self.robot.dot_w)[0]
                
                for l in range(self.robot.dimension):
                    V[:,0,l,N_joint*i+k] =  np.transpose(self.robot.arms[l].B0_r_i)[0]
                    V[:,1:4,l,N_joint*i+k] = self.robot.arms[l].A_i0
                    
        self.T, self.V, self.W = T, V, W
        print(time.time()-a)            
        
        
    def visualize(self, filename):   
     
        self.write_data(filename)
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        W = self.W
        ax.plot3D(W[0, :], W[1, :], W[2, :], 'gray')
        ax.set_xlabel('x /m')
        ax.set_ylabel('y /m')
        ax.set_zlabel('z /m')
        ax.set_title('Bahn im Arbeitsraum')
        fig.savefig('./temp.png')
        
    def write_data(self, filename):      
        T, V, N_q = self.T, self.V, self.robot.dimension
        with open(filename,'w',newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['# time\tB0_r_i(1)\tB0_r_i(2)\tB0_r_i(3)\t0_A_i(1,1)\t0_A_i(1,2)\t0_A_i(1,3)\t0_A_i(2,1)\t...'])
            for t in range(len(T)):
                curr_row = str(T[t])+'\t'
                for q in range(V.shape[2]):
                    for i in range(V.shape[1]):
                        for j in range(V.shape[0]):
                            curr_row = curr_row+str(V[j,i,q,t])+'\t'
                writer.writerow([curr_row])
            
                
    def dk_position_vectorchain(self):
    
        for i in range(self.robot.dimension):
        
       
            self.robot.arms[i].Bv_r_vi = np.array([[                      self.robot.arms[i].a],
                             [-np.sin(self.robot.arms[i].alpha)*self.robot.arms[i].d],
                             [np.cos(self.robot.arms[i].alpha)*self.robot.arms[i].d]])

            self.robot.arms[i].A_iv = np.matmul(Az(self.robot.q[self.robot.arms[i].number]), \
                                      Ax(self.robot.arms[i].alpha))
            
                
            if self.robot.arms[i].last == -1:
                self.robot.arms[i].A_i0 = self.robot.arms[i].A_iv
            else:
                self.robot.arms[i].A_i0 = np.matmul(self.robot.arms[i].A_iv, self.robot.arms[self.robot.arms[i].last].A_i0)
            
            if self.robot.arms[i].last == -1:
                self.robot.arms[i].Bi_r_i = np.matmul(self.robot.arms[i].A_iv, self.robot.arms[i].Bv_r_vi)
            else:
                self.robot.arms[i].Bi_r_i = np.matmul(self.robot.arms[i].A_iv, \
                                            (self.robot.arms[self.robot.arms[i].last].Bi_r_i + \
                                             self.robot.arms[i].Bv_r_vi))
            # print(self.robot.arms[i].Bi_r_i, self.robot.arms[i].Bv_r_vi)
            self.robot.arms[i].B0_r_i = np.matmul(np.transpose(self.robot.arms[i].A_i0), self.robot.arms[i].Bi_r_i)
        
        self.robot.w = np.matmul(np.transpose(self.robot.arms[len(self.robot.arms)-1].A_i0), \
                        (self.robot.arms[len(self.robot.arms)-1].Bi_r_i + self.robot.BN_r_N_tcp))
        # print(self.robot.arms[len(self.robot.arms)-1].Bi_r_i, self.robot.BN_r_N_tcp, self.robot.w)
        
    def calculate_dk_velocity(self):
        for i in range(self.robot.dimension):
            self.robot.arms[i].Bi_omega_rel =[[0],[0],[self.robot.dot_q[self.robot.arms[i].number]]]
            
            if self.robot.arms[i].last == -1:
                self.robot.arms[i].Bi_omega = self.robot.arms[i].Bi_omega_rel
            else:
                self.robot.arms[i].Bi_omega = np.matmul(self.robot.arms[i].A_iv, \
                                               self.robot.arms[self.robot.arms[i].last].Bi_omega)+\
                                               self.robot.arms[i].Bi_omega_rel
                                               
            if self.robot.arms[i].last == -1:
                self.robot.arms[i].Bi_dot_r_i = np.zeros((3,1))
            else:
                self.robot.arms[i].Bi_dot_r_i = np.matmul(self.robot.arms[i].A_iv, \
                                                (self.robot.arms[self.robot.arms[i].last].Bi_dot_r_i - \
                                                (np.matmul(tilde(self.robot.arms[i].Bv_r_vi), \
                                                 self.robot.arms[self.robot.arms[i].last].Bi_omega))))    

        self.robot.dot_w = np.matmul(np.transpose(self.robot.arms[self.robot.dimension-1].A_i0), \
                          (self.robot.arms[self.robot.dimension-1].Bi_dot_r_i - \
                          (np.matmul(tilde(self.robot.BN_r_N_tcp), \
                          self.robot.arms[self.robot.dimension-1].Bi_omega))))
            
            
            
            
            