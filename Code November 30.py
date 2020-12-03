# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 12:01:59 2020

@author: ghaza
"""

import math
import numpy as np
from numpy import linalg as LA
from sympy.abc import x, y, z, X, Y, Z, t, r
from sympy.solvers import solve
from scipy.optimize import minimize

"""
"""
#The function "find_h_rho" calculate the rho,h,M,exprime,eyprime,detector_orientation

def find_h_rho(position_particle,position_detector_face,position_detector_middle):
    
    distance_vector=[]
    detector_orientation=[]
    
    for i in range(0,3):
        distance_vector.append(position_particle[i]-position_detector_face[i])
        
    
    for i in range (0,3):
        detector_orientation.append(position_detector_face[i]-position_detector_middle[i])
        
    #ez'
    detector_orientation= detector_orientation/LA.norm(detector_orientation)
    K=np.dot(distance_vector,detector_orientation)
    h=np.fabs(K)   
    rho=np.fabs(LA.norm(K*detector_orientation-distance_vector))
    #vector of h
    M=h*detector_orientation
    #ex'
    xprime=K*detector_orientation-distance_vector
    exprime=xprime/LA.norm(xprime)
    #ey'
    eyprime=np.cross(detector_orientation,exprime)
    
    return rho,h,M,exprime,eyprime,detector_orientation

"""
"""

def solid_angle (position_particle,position_detector_face,
                 position_detector_middle):
    
    
    V_tetha=[]
    V_alpha=[]
    V_alpha_max=[]
    V_OA=[]
    V_OB=[]
    W_tetha=[]
    W_alpha=[]
    V_depth_detector=[]
    V_depth=[]
    
    r_cristal=1
    l_cristal=6
    r_reactor=8
    mu_d=0.21477
    mu_r=4
    #R=0.01
    #tau=10
    

    rho,h,M,exprime,eyprime,detector_orientation=find_h_rho(position_particle,
                                        position_detector_face,position_detector_middle)
    
    #efficiency=0
    a=0
    b=0
    for i in range (10):
        
        if rho < r_cristal:
           tetha_max = np.arctan((r_cristal+rho)/h)
           tetha_cri = np.arctan((r_cristal-rho)/h)
           tetha_min = 0.0
           tetha=np.arccos(math.cos(tetha_min)- np.random.random_sample() * 
                             (math.cos(tetha_min)-math.cos(tetha_max)))
           V_tetha.append(tetha)            
           W_tetha.append((math.cos(tetha_min) - math.cos(tetha_max))/2)
           
           if tetha < tetha_cri :
                alpha_max = math.pi
                V_alpha_max.append(alpha_max)
                alpha=alpha_max*(2* np.random.random_sample()-1)
                V_alpha.append(alpha)
                W_alpha.append(1)
                
                OA = rho * math.cos(alpha) + pow(pow(r_cristal,2) - 
            (pow(rho,2)*pow(math.sin(alpha),2)),0.5)
                V_OA.append(OA)
           else:
                 alpha_max = math.acos((pow(rho, 2) + pow((h * math.tan(tetha)), 2) -
                                        pow(r_cristal, 2)) / (2 * h * rho * math.tan(tetha)))
                 V_alpha_max.append(alpha_max)
                 alpha=alpha_max*(2* np.random.random_sample()-1)
                 V_alpha.append(alpha)
                 W_alpha.append(alpha_max/math.pi)
                 OA = rho * math.cos(alpha) + pow(pow(r_cristal,2) - 
            (pow(rho,2)*pow(math.sin(alpha),2)),0.5)
                 V_OA.append(OA)
                 
           depth_detector=path_length_detector_one(OA,h,rho,r_cristal,l_cristal,alpha,tetha)
           V_depth_detector.append(depth_detector)
            

        else:
            alpha_max = np.arcsin(r_cristal/rho)
            alpha=alpha_max*(2*np.random.random_sample()-1)
            V_alpha.append(alpha)
            W_alpha.append(alpha_max/math.pi)
            
            OB = rho * math.cos(alpha) - pow(pow(r_cristal,2) - 
            (pow(rho,2)*pow(math.sin(alpha),2)),0.5)
            V_OB.append(OB)
            
            
            OA = rho * math.cos(alpha) + pow(pow(r_cristal,2) - 
            (pow(rho,2)*pow(math.sin(alpha),2)),0.5)
            V_OA.append(OA)
            
            tetha_min = np.arctan(OB/(h+l_cristal))          
            tetha_max = np.arctan(OA/h)  
            tetha_cri = np.arctan(OB/(h+l_cristal))
            tetha = np.arccos(math.cos(tetha_min) - np.random.random_sample() * (math.cos(tetha_min)
            - math.cos(tetha_max)))
            V_tetha.append(tetha)
            W_tetha.append (math.cos(tetha_min) - math.cos(tetha_max)/2)
            depth_detector=path_length_detector_two(OA,OB,rho,r_cristal,l_cristal,alpha,tetha,h,tetha_cri)
            V_depth_detector.append(depth_detector)
        
        depth=path_length_reactor(alpha,tetha,r_reactor,position_particle,position_detector_face,
                                  position_detector_middle)
        V_depth.append(depth)
        
        a+=W_alpha[i]*W_tetha[i]*(1-np.exp(-1*mu_d*depth_detector))
        b+=W_alpha[i]*W_tetha[i]*(1-np.exp(-1*mu_d*depth_detector))*np.exp(-1*depth)
    a=a/i
    b=b/i
    print (depth)

    return (a,b,i)
    
    

    
        
"""        
"""

def path_length_detector_one(OA,h,rho,r_cristal,l_cristal,alpha,tetha):
    #Point source viewing only the top of the detector
    tetha_one= np.arctan(OA/(h+l_cristal))
    tetha_two= np.arctan(OA/h)
    if tetha < tetha_one:
        depth_detector=l_cristal/np.cos(tetha)
    elif tetha>tetha_one and tetha<tetha_two:
        depth_detector=(OA/np.sin(tetha))-(h/np.cos(tetha))
    return depth_detector

"""
"""
 
def path_length_detector_two(OA,OB,rho,r_cristal,l_cristal,alpha,tetha,h,tetha_cri):
    #Point source viewing the top and the lateral surface of the detector
    if tetha<tetha_cri:
        if (h+l_cristal)*np.tan(tetha)<OA:
            depth_detector=((h+l_cristal)/np.cos(tetha))-OB/np.sin(tetha)
        else:
            depth_detector=(OA-OB)/np.sin(tetha)
    else:
         if (h+l_cristal)*np.tan(tetha)<OA:
             depth_detector=l_cristal/np.cos(tetha)
         else:
             depth_detector=(OA/np.sin(tetha))-h/np.cos(tetha)
    return depth_detector
                         

"""
"""

def path_length_reactor(alpha,tetha,r_reactor,position_particle,position_detector_face,
                        position_detector_middle):
    
      
        
        rho,h,M,exprime,eyprime,detector_orientation=find_h_rho(position_particle,
                                        position_detector_face,position_detector_middle)
        
        newOrigin=[]
        position_particle_translation=[]
        intersection_point=[]
        intersection_point_1=[]
        intersection_point_2=[]
        distance_vector_1=[]
        distance_vector_2=[]
        depth_vector=[]
        
        for i in range (0,3):
            newOrigin.append(position_particle[i]-M[i])
        #print(newOrigin)     
   
    #Translation
        for i in range (0,3):
            position_particle_translation.append(position_particle[i]-newOrigin[i])
        
    
        A=np.array((exprime,eyprime,detector_orientation))
        A_inv=np.linalg.inv(A)
        #print(A_inv)
        B=np.array([position_particle_translation[0],position_particle_translation[1]
        ,position_particle_translation[2]])
       
    
    #Transformation
        particle_new=np.matmul(A,B)
        #print(particle_new)
        
        X=particle_new[0]+t*math.sin(tetha)*math.cos(alpha)
        Y=particle_new[1]+t*math.sin(alpha)*math.sin(tetha)
        Z=particle_new[2]+t*math.cos(math.pi-tetha)
        
        
        
        x=newOrigin[0]
        y=newOrigin[1]
        z=newOrigin[2]
        #print(newOrigin)

        
        B=[X,Y,Z]
        
        #print(B)
        
        for i in range (0,3):
                x+=A_inv[0][i]*B[i]
                y+=A_inv[1][i]*B[i]
                z+=A_inv[2][i]*B[i]
                
        intersection_point.append(x)
        intersection_point.append(y)
        intersection_point.append(z)
        
        
        func=pow(x,2)+pow(y,2)-pow(r_reactor,2)
        M=solve(func,t)

      
        for i in range(0,3):
            
            point_1=intersection_point[i].subs(t,M[0])
            intersection_point_1.append(point_1)
            point_2=intersection_point[i].subs(t,M[1])
            intersection_point_2.append(point_2)
        
        
        for i in range(0,3):
            distance_vector_1.append(position_detector_face[i]-intersection_point_1[i])
            distance_vector_2.append(position_detector_face[i]-intersection_point_2[i])
            
        
         #Convert sympy float to normal float type  
        distance_vector_1=str(distance_vector_1)
        distance_vector_1=eval(distance_vector_1)
        distance_vector_2=str(distance_vector_2)
        distance_vector_2=eval(distance_vector_2)
        
        distance_1=np.fabs(LA.norm(distance_vector_1))
        distance_2=np.fabs(LA.norm(distance_vector_2))
        
        
        if distance_2>distance_1:
            for i in range(0,3):
                depth_vector.append(position_particle[i]-intersection_point_1[i])
        else:
            for i in range(0,3):
                depth_vector.append(position_particle[i]-intersection_point_2[i])
                
        depth_vector=str(depth_vector)
        depth_vector=eval(depth_vector)
        depth=np.fabs(LA.norm(depth_vector))       

           
        return depth

"""
"""
#Minimization by Nelder-Mead
#I should continue from here
"""
def fun(x):
    T=1
    nui=1
    phi=1
    position_particle=np.array([0,1,2])
    position_detector_face=np.array([8.4,1,5])
    position_detector_middle=np.array([12,1,5])
    #here instead of gining the values how python can read it from the mesh or how I can just call the a,b from solidangle 
    #without giving the position of particle and detector
    a,b,i=solid_angle (position_particle,position_detector_face,position_detector_middle)
    print(a,b)
    g=T*nui*x[0]*phi*((x[1]/1)*a+b)/((1+x[2]*nui)*phi*(x[1]/1)*a+b)
    return (g)


x0=np.array([0,0,0])
res=minimize(fun,x0,method='nelder-mead',options={'xatol': 1e-8, 'disp': True})
print(res.x)

#x0=np.array([1,1,1])
#res = minimize(minimization, x0, method='nelder-mead',
               #options={'xatol': 1e-8, 'disp': True})
"""
#find_h_rho([0,1,7],[8.4,1,5],[12,1,5])
#transformation([-1,2,3],[-4,3,1],[-5,3,2])
#path_length_reactor(0,0,8,[7.9,1,0],[10,1,5],[12,1,5])
solid_angle([-1,-1,-1],[8.4,1,5],[12,1,5])
#T*nui*x[0]*phi*((x[1]/1)*a+b)/(1+x[2]*nui)*phi*((x[1]/1)*a+b
