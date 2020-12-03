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
#Open the file the that includes particle positions
position_particle=np.loadtxt("./Testdata.txt")

#The function "find_h_rho" calculate the rho and h for all the positions
def find_h_rho(position_detector_face,position_detector_middle):
    
    detector_orientation=[]
    V_rho=[]
    V_h=[]
    
    #This "for" loop calculates the orientation of the detector using two points on the detector axis
    for i in range (0,3):
        detector_orientation.append(position_detector_face[i]-position_detector_middle[i])
        
    detector_orientation= detector_orientation/LA.norm(detector_orientation)
    
    
    #Loop over all the positions to calculate rho and h  
    for j in range (0,len(position_particle)):
        
        #Distance vector is a vector from a point on the detector face to the particle position
        distance_vector=[]
        for i in range (0,3):
            distance_vector.append(position_particle[j][i]-position_detector_face[i])

        #Projection of the distance vector on the detector axis which is parallel to the z' axis
        K=np.dot(distance_vector,detector_orientation)
        h=np.fabs(K) 
        V_h.append(h)
        
        rho=np.fabs(LA.norm(K*detector_orientation-distance_vector))
        V_rho.append(rho)

    return V_rho,V_h
#rho and h as output of this function will be used in solid angle calculation   


"""
"""

def solid_angle (position_detector_face,
                 position_detector_middle):
    r_cristal=1
    l_cristal=6
    r_reactor=8
    Monte_Carlo_iteration=5
    
    """
    The four following vectors (V_...) store weight factors for alpha and theta, and the path length
    in detector and reactor for all the position. They contain vectors that stores 
    the mentioned parameters for each iteration of the Monte Carlo
    """
    V_W_tetha=[]
    V_W_alpha=[]
    V_V_depth_detector=[]
    V_V_depth=[] 
    
    
    #Call the function that calculates rho and h for each position
    V_rho,V_h=find_h_rho(position_detector_face,position_detector_middle)
    
    #Loop over all the position to do the Monte Carlo calculation
    for j in range (0,len(position_particle)):
        
        """
        W_tetha and W_alpha store the weight factor for each iteration of the Monte Carlo
        for a specific position
        """
        W_tetha=[]
        W_alpha=[]
        
        """
        The following vectors starting with (V_...) store related parameter for each iteration
        of the Monte Carlo for a specific position
        """
        V_tetha=[]
        V_alpha=[]
        V_alpha_max=[]
        V_OA=[]
        V_OB=[]
        V_depth_detector=[]
        V_depth=[]
    
        #Start the Monte Carlo
        
        for i in range (Monte_Carlo_iteration):
            
            """
            First "if" defines that the particle position is in the detector's face 
            or it can view the detector from both lateral side and top.
            """
            if V_rho[j] < r_cristal:
                #Check from here
               tetha_max = np.arctan((r_cristal+V_rho[j])/V_h[j])
               tetha_cri = np.arctan((r_cristal-V_rho[j])/V_h[j])
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
                
                  OA = V_rho[j] * math.cos(alpha) + pow(pow(r_cristal,2) - 
                  (pow(V_rho[j],2)*pow(math.sin(alpha),2)),0.5)
                  V_OA.append(OA)
               else:
                    alpha_max = math.acos((pow(V_rho[j], 2) + pow((V_h[j] * math.tan(tetha)), 2) -
                                        pow(r_cristal, 2)) / (2 * V_h[j] * V_rho[j] * math.tan(tetha)))
                    V_alpha_max.append(alpha_max)
                    alpha=alpha_max*(2* np.random.random_sample()-1)
                    V_alpha.append(alpha)
                    W_alpha.append(alpha_max/math.pi)
                    OA = V_rho[j] * math.cos(alpha) + pow(pow(r_cristal,2) - 
                        (pow(V_rho[j],2)*pow(math.sin(alpha),2)),0.5)
                    V_OA.append(OA)
                 
               depth_detector=path_length_detector_one(OA,V_h[j],V_rho[j],r_cristal,l_cristal,alpha,tetha)
               V_depth_detector.append(depth_detector)
               depth=path_length_reactor(alpha,tetha,r_reactor,position_detector_face,
                                  position_detector_middle,position_particle[j])
               V_depth.append(depth)
               
            

            else:
                alpha_max = np.arcsin(r_cristal/V_rho[j])
                alpha=alpha_max*(2*np.random.random_sample()-1)
                V_alpha.append(alpha)
                W_alpha.append(alpha_max/math.pi)
            
                OB = V_rho[j] * math.cos(alpha) - pow(pow(r_cristal,2) - 
                 (pow(V_rho[j],2)*pow(math.sin(alpha),2)),0.5)
                V_OB.append(OB)
            
            
                OA = V_rho[j] * math.cos(alpha) + pow(pow(r_cristal,2) - 
             (pow(V_rho[j],2)*pow(math.sin(alpha),2)),0.5)
                V_OA.append(OA)
            
                tetha_min = np.arctan(OB/(V_h[j]+l_cristal))          
                tetha_max = np.arctan(OA/V_h[j])  
                tetha_cri = np.arctan(OB/(V_h[j]+l_cristal))
                tetha = np.arccos(math.cos(tetha_min) - np.random.random_sample() * (math.cos(tetha_min)
                   - math.cos(tetha_max)))
                V_tetha.append(tetha)
                W_tetha.append (math.cos(tetha_min) - math.cos(tetha_max)/2)
                depth_detector=path_length_detector_two(OA,OB,V_rho[j],r_cristal
                                                        ,l_cristal,alpha,tetha,V_h[j],tetha_cri)
                V_depth_detector.append(depth_detector)
                
        
            depth=path_length_reactor(alpha,tetha,r_reactor,position_detector_face,
                                  position_detector_middle,position_particle[j])
            V_depth.append(depth)
            
        V_W_tetha.append(W_tetha)
        V_W_alpha.append(W_alpha)
        V_V_depth_detector.append(V_depth_detector)
        V_V_depth.append(V_depth)
      
          
     
        #Moshakhasan code in ja moshkel dare halesh kon
    #print(V_W_tetha)
    #print(V_W_tetha)
    return (V_W_tetha)
    
    

    
        
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

def path_length_reactor(alpha,tetha,r_reactor,position_detector_face,
                        position_detector_middle,position_particle):
    
      
        #print(position_particle)
        #rho,h,detector_orientation=find_h_rho(position_detector_face,position_detector_middle)
        detector_orientation=[]
        newOrigin=[]
        position_particle_translation=[]
        intersection_point=[]
        intersection_point_1=[]
        intersection_point_2=[]
        distance_vector_1=[]
        distance_vector_2=[]
        depth_vector=[]
        detector_orientation.append(position_detector_face[0]-position_detector_middle[0])
        detector_orientation.append(position_detector_face[1]-position_detector_middle[1])
        detector_orientation.append(position_detector_face[2]-position_detector_middle[2])
        detector_orientation= detector_orientation/LA.norm(detector_orientation)
        
        distance_vector=[]
        distance_vector.append(position_particle[0]-position_detector_face[0])
        distance_vector.append(position_particle[1]-position_detector_face[1])
        distance_vector.append(position_particle[2]-position_detector_face[2])
        K=np.dot(distance_vector,detector_orientation)
        #K=np.dot(distance_vector,detector_orientation)
        h=np.fabs(K) 
        rho=np.fabs(LA.norm(K*detector_orientation-distance_vector))
        
    #vector of h
    
        M=h*detector_orientation
        
    #ex'
        xprime=K*detector_orientation-distance_vector
        exprime=xprime/LA.norm(xprime)
        
    #ey'
        eyprime=np.cross(detector_orientation,exprime)
        
        
      
        
        for i in range (0,3):
           
            #print(position_particle[i])
            
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
find_h_rho([4,1,5],[6,1,5])
#transformation([-1,2,3],[-4,3,1],[-5,3,2])
#path_length_reactor(0,0,8,[7.9,1,0],[10,1,5],[12,1,5])
solid_angle([4,1,5],[6,1,5])
#T*nui*x[0]*phi*((x[1]/1)*a+b)/(1+x[2]*nui)*phi*((x[1]/1)*a+b
