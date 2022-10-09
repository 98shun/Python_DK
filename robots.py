# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 17:04:11 2022

@author: Yuanshun Xie
"""

import numpy as np

class arm:
    def __init__(self, alpha, a, d, i):
        self.alpha=alpha
        self.a=a
        self.d=d
        self.Bv_r_vi=0
        self.A_iv=0
        self.A_i0=0
        self.Bi_r_i=0
        self.B0_r_i=0
        self.number=i
        self.last=i-1
        
        self.Bi_omega_rel=0
        self.Bi_omega=0
        self.Bi_dot_r_i=0


    
class robot:
    """
    In this class, the model of the robot will be established and all the structural variables are listed.
    """
    
    def __init__(self, n_arms,DH,dt,BN_r_N_tcp):
          
# Structural information about the robot    
        self.dimension = n_arms # the dimension of the work space is given and stored in the robot
        self.time = 0
        self.dt=dt #Time step size

# Length, width and thickness of the single arm
        self.length = np.zeros((1,n_arms))
        self.width = np.zeros((1,n_arms))
        self.thickness = np.zeros((1,n_arms))

# DH parameters
        self.arms=DH
        
#joint angle description
        self.q=np.zeros((self.dimension))
        self.dot_q=0
        
#position of the TCP in B0 system
        self.w=0
        self.dot_w=0
        self.BN_r_N_tcp=BN_r_N_tcp
        