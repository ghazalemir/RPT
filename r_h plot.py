# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 16:16:09 2020

@author: ghaza
"""
import math
import numpy as np
from numpy import linalg as LA
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



def find_h_rho(position_particle,position_detector_face,position_detector_middle):
    #position_particle=p
    #position_detector_face=d
    #position_detector_middle=dd
    #Xpd
    distance_vector=[position_particle[0]-position_detector_face[0],position_particle[1]-position_detector_face[1],position_particle[2]-position_detector_face[2]]
    #ez'
    detector_orientation=[(position_detector_face[0]-position_detector_middle[0]),position_detector_face[1]-position_detector_middle[1],position_detector_face[2]-position_detector_middle[2]]
    norm1=LA.norm(detector_orientation)
    detector_orientation= detector_orientation/norm1
    
    
    #z'=Z
    Z=np.dot(distance_vector,detector_orientation)
    rho=np.fabs(Z)
    #ey'
    yprime=distance_vector-Z*detector_orientation
    
    norm2=LA.norm(yprime)
    eyprime=yprime/norm2
    Y=np.linalg.norm(yprime)
    h=np.fabs(Y)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    x=[position_detector_face[0],position_detector_face[0],position_detector_face[0],position_detector_face[0]]
    y=[position_detector_face[1],position_detector_face[1],position_detector_face[1],position_detector_face[1]]
    z=[position_detector_face[2],position_detector_face[2],position_detector_face[2],position_detector_face[2]]
    u=[distance_vector[0],detector_orientation[0],yprime[0],eyprime[0]]
    v=[distance_vector[1],detector_orientation[1],yprime[1],eyprime[1]]
    w=[distance_vector[2],detector_orientation[2],yprime[2],eyprime[2]]
    ax.quiver(x, y, z, u, v, w)
    ax.set_xlim([-4, 4])
    ax.set_ylim([-4, 4])
    ax.set_zlim([-4, 4])
    
    plt.show()
    
    print(distance_vector,yprime,detector_orientation,eyprime)

    return Z,Y,distance_vector,yprime,detector_orientation,eyprime



find_h_rho([0,0,1],[2,0,3],[4,0,3])



