# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 16:36:07 2020

@author: ghaza
"""

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
    #print(detector_orientation)
    
    #z'=Z
    Z=np.dot(distance_vector,detector_orientation)
    h=np.fabs(Z)
    M=h*detector_orientation
    #print(detector_orientation)
    #print(M)
    #print(Z)
    
    #ex'
    xprime=Z*detector_orientation-distance_vector
    
    X=LA.norm(xprime)
    exprime=xprime/X
    #print(exprime)
    #ey'
    eyprime=np.cross(detector_orientation,exprime)
    rho=np.fabs(X)
    #New coordinate calculation
    x_new=position_particle[0]-M[0]
    y_new=position_particle[1]-M[1]
    z_new=position_particle[2]-M[2]
    new=[x_new,y_new,z_new]
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    x=[position_detector_face[0],new[0],new[0],new[0]]
    y=[position_detector_face[1],new[1],new[1],new[1]]
    z=[position_detector_face[2],new[2],new[2],new[2]]
    u=[distance_vector[0],detector_orientation[0],exprime[0],eyprime[0]]
    v=[distance_vector[1],detector_orientation[1],exprime[1],eyprime[1]]
    w=[distance_vector[2],detector_orientation[2],exprime[2],eyprime[2]]
    ax.quiver(x, y, z, u, v, w)
    ax.set_xlim([-4, 4])
    ax.set_ylim([-4, 4])
    ax.set_zlim([-4, 4])
    
    plt.show()
    
    print(detector_orientation,exprime,eyprime,new,rho,h)

    return Z,X,distance_vector,xprime,detector_orientation,exprime



find_h_rho([1,2,3],[2,2,4],[3,2,4])




