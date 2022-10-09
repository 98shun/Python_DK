# -*- coding: utf-8 -*-
"""
Created on Sat Oct  8 00:26:26 2022

@author: Yuanshun Xie
"""
from robots import robot, arm
from kinematics import *
import numpy as np


#define a robot with 6 arms
n_arms=6

        
DH = np.array([[0, np.pi/2, np.pi, -np.pi/2, np.pi/2, -np.pi/2 ], #here are the alphas
      [0, 0, 0.184, 0.103, 0, 0],#here are the a
      [0, 0, 0, 0.196, 0, 0.058]])#here are the d

DH_matrix=[]
for i in range(n_arms):
    DH_matrix.append(arm(DH[0,i], DH[1,i], DH[2,i],i))
    
    
dt=0.001    
T_total=2.0
BN_r_N_tcp=np.array([[0],[0],[0.035]])

my_rob=robot(n_arms, DH_matrix, dt,BN_r_N_tcp)

q_e = np.pi/4

my_kinematic=kinematic(my_rob, T_total, q_e)
my_kinematic.calculate_jointwinkel()
my_kinematic.visualize('trajectory_D_python.csv')
