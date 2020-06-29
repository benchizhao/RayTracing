# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 15:52:11 2020

@author: Benchi Zhao
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 16:42:19 2020

@author: Benchi Zhao
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

class RayTracing:

    def __init__(self,x,y,theta):
        '''
        __init__ (self,x,y,theta)
            Gives the initial state of the ray.
            
        Parameters
        ------------
        self.x: folat
            Initial x-position of the ray.
        self.y: float
            Initial y-position of the ray.
        self.theta: float
            The angle between the horizontal and the ray path (in degree). 
        self.state: list
            When ray interacting with the optical equipments, the ray state will change,
            all states are recorded in self.state.
        '''
        self.x = x
        self.y = y
        self.theta = theta
        self.state = []
        self.n_1 = 1
        self.n_2 = 1.5
        self.M = False
        self.Gauss = True
        
    def slope2rad(self,slope):
        return np.arctan(slope)
    
    def angle2slope(self,angle):
        return np.tan(np.deg2rad(angle))
    
    def gaussian(self,position):
        return 1/(np.sqrt(2*math.pi)) * math.exp(-position**2/2) 
    
    def R(self):
        theta_i = self.slope2rad(self.state[-1][2])
        theta_t = np.sqrt(1-(self.n_1/self.n_2 * np.sin(theta_i))**2)
        result = (abs((self.n_1*np.cos(theta_t)-self.n_2*np.cos(theta_i))/(self.n_1*np.cos(theta_t)+self.n_2*np.cos(theta_i))))**2
        return result
    
    def T(self):
        return 1-self.R()
    
    def hit_point(self,r):
        theta = self.slope2rad(self.state[-1][2])
        y = self.state[-1][1]
        m = y*np.cos(math.pi/2-theta)+np.sqrt(y*y*np.cos(math.pi/2-theta)**2 - (y**2-r**2))
        x = m*np.sin(math.pi/2-theta)
        y = m*np.sin(theta)
        
        self.state[-1][0] -= x
        self.state[-1][1] -= y
        self.state[-1][4] -= x
        self.state[-1][5] -= y
        pass
    
        
    def ray(self):
        '''
        ray(self)
            Append the initial ray state into the total ray state.
            8 elements are added into array, firt 4 element describe the refraction
            second 4 describe the reflection.
        '''
        if self.Gauss == True:
            ray_state = np.array([self.x,self.y,self.theta,self.gaussian(self.y),0,0,0,0])
            self.state.append(ray_state)
        else:
            ray_state = np.array([self.x,self.y,self.theta,1,0,0,0,0])
            self.state.append(ray_state)
        
    def free_propagate(self,distance):
        if self.M == False:        
            ray_state_1 = np.dot(np.array([[1,distance],[0,1]]),self.state[-1][1:3])
            ray_state_1 = np.append(self.state[-1][0]+distance,ray_state_1)
            ray_state_1 = np.append(ray_state_1,self.state[-1][3])
            
            ray_state_2 = np.dot(np.array([[1,distance],[0,1]]),self.state[-1][5:7])
            ray_state_2 = np.append(self.state[-1][4]+distance,ray_state_2)
            ray_state_2 = np.append(ray_state_2,self.state[-1][7])
            
            ray_state = np.append(ray_state_1,ray_state_2)         
        else:           
#            ray_state_2 = np.dot(np.array([[1,-distance],[0,1]]),self.state[-1][5:7])
#            ray_state_2 = np.append(self.state[-1][4]-distance,ray_state_2)
#            ray_state_2 = np.append(ray_state_2,self.state[-1][-1])
            
            ray_state_1 = np.dot(np.array([[1,distance],[0,1]]),self.state[-1][1:3])
            ray_state_1 = np.append(self.state[-1][0]+distance,ray_state_1)
            ray_state_1 = np.append(ray_state_1,self.state[-1][3])
            
            ray_state_2 = np.dot(np.array([[1,-distance],[0,1]]),self.state[-1][5:7])
            ray_state_2 = np.append(self.state[-1][4]-distance,ray_state_2)
            ray_state_2 = np.append(ray_state_2,self.state[-1][-1])
        
            ray_state = np.append(ray_state_1,ray_state_2)
        self.state.append(ray_state)
        return ray_state
    
    
    def full_reflection(self):
#        self.M = True
        ray_state = np.dot(np.array([[1,0],[0,-1]]),self.state[-1][1:3])
        ray_state = np.append(self.state[-1][0],ray_state)
        ray_state = np.append(ray_state,self.state[-1][-1])
        return ray_state
        
    def curved_interface(self,r):
        self.M = False
        ray_state_1 = np.dot(np.array([[1,0],[((self.n_1-self.n_2)/(r*self.n_2)),self.n_1/self.n_2]]),self.state[-1][1:3])
        ray_state_1 = np.append(self.state[-1][0],ray_state_1)
        ray_state_1 = np.append(ray_state_1,self.state[-1][3]*self.T())
    
        ray_state_2 = self.full_reflection()
        ray_state_2[-1] = ray_state_2[-1]*self.R()
        
        ray_state = np.append(ray_state_1,ray_state_2)
        
        
        theta = self.slope2rad(self.state[-1][2])
        y = ray_state[1]
        m = y*np.cos(theta)+np.sqrt(y*y*np.cos(theta)**2 - (y**2-r**2))
        x_0 = m*np.sin(math.pi/2 - theta)
        y_0 = m*np.sin(theta)
        
        
        print(m,x_0,y_0)
#        
#        self.state[-1][0] -= x_0
#        self.state[-1][1] -= y_0
#        self.state[-1][4] -= x_0
#        self.state[-1][5] -= y_0
        self.state.append(ray_state)
        return ray_state
       

if __name__=='__main__':
    
    def bundle(rays):
        width = 1
        distribution = np.linspace(-width,width,rays)
        bundle_of_ray = []
        for i in distribution:
            RT = RayTracing(0,i,0.05)
            RT.ray()
            RT.free_propagate(5)
            RT.curved_interface(2)
#            RT.free_propagate(5)
               
            bundle_of_ray.append(RT.state)
        return bundle_of_ray
    
    def plot(rays):
        '''
        plot_ray(self)
            Plot the ray path which is described in self.state.
        '''
        Ray_thickness = 15
        
        for j in range(len(rays)):
            for i in range(len(rays[j])-1):            
                x_1 = np.linspace(rays[j][i][0],rays[j][i+1][0])
                y_1 = np.linspace(rays[j][i][1],rays[j][i+1][1])
                plt.plot(x_1,y_1,'k-',linewidth = rays[j][i][3]*Ray_thickness)  
                
                x_2 =  np.linspace(rays[j][i][4],rays[j][i+1][4])
                y_2 = np.linspace(rays[j][i][5],rays[j][i+1][5])
                plt.plot(x_2,y_2,'k--',linewidth = rays[j][i][-1]*Ray_thickness)  
            x_main = np.linspace(0,max(x_1))
            y_main = [0]*len(x_main)
            plt.xlim(0,max(x_1))
            plt.plot(x_main,y_main,'--')
            
            #plot circle
            theta = np.linspace(0, 2*np.pi, 100)
            r = 2
            x1 = r*np.cos(theta)+5
            x2 = r*np.sin(theta)
            plt.plot(x1,x2)
            
        plt.show()
        
    def make_table(n):
        for i in range(n):            
            data = bundle(n)[i]
            df = pd.DataFrame(data,columns=list('ABCDEFGH'))
            print(df)
            
    make_table(1)
    plot(bundle(1))