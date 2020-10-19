# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 10:51:43 2020

@author: ghaza
"""

import math
#import utils as utils
import numpy as np
rand_alpha = np.random.random(5)
rand_tetha = np.random.random(5)
#print(rand_alpha)
#print(rand_tetha)
V_tetha = []
V_alphamax = []
V_alpha = []
V_OB = []
V_OA = []
V_tethamax = []
V_tethamin = []
W_tetha = []
W_alpha = []
def solid_angle (rho,h,r_cristal,l_cristal):
    
    if rho < r_cristal:
        
        tetha_max = np.arctan((r_cristal+rho)/h)
        tetha_cri = np.arctan((r_cristal-rho)/h)
        tetha_min = 0.0
    #else:
        #alpha_max = np.arcsin(r_cristal/rho)
        #V_alphamax.append(alpha_max)
    
    if rho < r_cristal:
        for i in range (5):
            tetha = np.arccos(math.cos(tetha_min)-rand_tetha[i]*(math.cos(tetha_min)-math.cos(tetha_max)))
            V_tetha.append(tetha)
            W_tetha.append((math.cos(tetha_min) - math.cos(tetha_max))/2)
            if tetha < tetha_cri :
                alpha_max = 2 * math.pi
                V_alphamax.append(alpha_max)
                W_alpha.append(1)
                # I am not sure that is 2pi or pi
            else:
                alpha_max = np.arccos(rho**2 - r_cristal**2 + (h**2) * (math.tan(tetha)**2) / 2*h*rho*math.tan(tetha))
                V_alphamax.append(alpha_max)
                W_alpha.append(V_alphamax[i]/math.pi)
    else:
        
        #print(V_alphamax)
        for i in range (5):
            alpha_max = np.arcsin(r_cristal/rho)
            V_alphamax.append(alpha_max)
            W_alpha.append(V_alphamax[i]/math.pi)
            
            alpha = alpha_max * (2 * rand_alpha[i] - 1)
            V_alpha.append(alpha)
            
            OB = rho * math.cos(V_alpha[i]) - pow(pow(r_cristal,2) - (pow(rho,2)*pow(math.sin(V_alpha[i]),2)),0.5)
            V_OB.append(OB)
            
            OA = rho * math.cos(V_alpha[i]) + pow(pow(r_cristal,2) - (pow(rho,2)*pow(math.sin(V_alpha[i]),2)),0.5)
            V_OA.append(OA)
            
            tetha_min = np.arctan(V_OB[i]/(h+l_cristal))
            V_tethamin.append(tetha_min)
            
            tetha_max = np.arctan(V_OA[i]/h)
            V_tethamax.append(tetha_max)
            
            tetha = np.arccos(math.cos(V_tethamin[i]) - rand_tetha * (math.cos(V_tethamin[i])- math.cos(V_tethamax[i])))
            V_tetha.append(tetha)
            W_tetha.append((math.cos(tetha_min) - math.cos(tetha_max))/2)
    print(W_tetha)
    print(W_alpha)
            

            
            
            
solid_angle(1,1,0.5,1)

    

        