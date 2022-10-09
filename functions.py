# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 23:47:39 2022

@author: Yuanshun Xie
"""

import numpy as np

def Ax(phi):
    A = np.array([[1,         0,        0],
        [0,  np.cos(phi), np.sin(phi)],
        [0, -np.sin(phi), np.cos(phi)]])
    return A

def Ay(phi):
    A =  np.array([[np.cos(phi), 0, -np.sin(phi)],
         [       0,    1,         0],
         [np.sin(phi), 0,  np.cos(phi)]])
    return A

def Az(phi):
    A = np.array([[np.cos(phi), np.sin(phi), 0],
         [-np.sin(phi), np.cos(phi), 0],
         [       0,        0, 1]])
    return A

def dot_Ax(phi, dot_phi):
    dotA=np.array([[0,                0,                 0      ],
          [0,-np.sin(phi)*dot_phi,  np.cos(phi)*dot_phi],
          [0,-np.cos(phi)*dot_phi, -np.sin(phi)*dot_phi]])
    return dotA

def dot_Ay(phi, dot_phi):  
    dotA = np.array([[ -np.sin(phi)*dot_phi, 0, -np.cos(phi)*dot_phi],
            [  0,                   0,                 0   ],
            [  np.cos(phi)*dot_phi, 0, -np.sin(phi)*dot_phi]])
    return dotA

def dot_Az(phi, dot_phi):   
    dotA = np.array([[ -np.sin(phi)*dot_phi,  np.cos(phi)*dot_phi, 0],
            [ -np.cos(phi)*dot_phi, -np.sin(phi)*dot_phi, 0],
            [  0,                   0,                    0]])
    return dotA

def dh_trafo(alpha, a, d, theta):
    # DH-Matrix
    D_vi = np.array([[np.cos(theta),            -np.sin(theta),          0,                a                    ],
            [np.sin(theta)*np.cos(alpha), np.cos(theta)*np.cos(alpha), -np.sin(alpha), -np.sin(alpha)*d],
            [np.sin(theta)*np.sin(alpha), np.cos(theta)*np.sin(alpha), np.cos(alpha), np.cos(alpha)*d  ],
            [0,                           0,                           0,          1                   ]])
    return D_vi

def tilde(v):
    T = np.array([[    0, -v[2],  v[1]],
          [v[2],     0, -v[0]],
         [-v[1],  v[0],    0]])
    return T

