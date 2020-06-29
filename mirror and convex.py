# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 09:06:39 2020

This .py file simulate the behavior of the plane mirror and lens. To trace the
ray, we use the state of ray to describe. In the end we also plot the ray 
path.

@author: Benchi Zhao
"""
import numpy as np
import matplotlib.pyplot as plt
import math

class RayTracing:

    def __init__(self,x,z,theta):
        '''
        __init__ (self,x,z,theta)
            Gives the initial state of the ray.
            
        Parameters
        ------------
        self.x: folat
            Initial x-position of the ray.
        self.z: float
            Initial z-position of the ray.
        self.theta: float
            The angle between the horizontal and the ray path (in degree). 
        self.state: list
            When ray interacting with the optical equipments, the ray state will change,
            all states are recorded in self.state.
        '''
        self.x = x
        self.z = z
        self.theta = theta
        self.state = []
        
    def ray(self):
        '''
        ray(self)
            Append the initial ray state into the total ray state.
        '''
        ray_state = np.array([self.x,self.z,self.theta])
        self.state.append(ray_state)
    
    def lens(self,focal_length,position):
        '''
        lens(self,focal_length,position)
            Simulate the behavior of lens for both convex and concave lens.
            Append the ray state into self.state after passing the lens .
            
        Parameters
        ------------
        focal_length: float
            The distance from focal point to central lens. 
            If >0 it is convex lens, if <0 it is concave lens.
        position: float
            The x position of central lens.
        '''
        slope = (self.state[-1][2]/180)*math.pi
        u = position
        f = focal_length
        hit_point = np.array([u,slope*u + self.state[-1][1]])
        
        if u == f:
            output_slope = -self.state[-1][1]/f
            output_angle = np.arctan(output_slope)/math.pi * 180
            ray_state = np.array([[hit_point[0],hit_point[1],output_angle]])         
        else: 
            v = u*f/(u-f)
            x = v + u
            y = -self.state[-1][1]/u * x + self.state[-1][1]
            image_position = np.array([x,y])
            output_slope = (image_position[1]-hit_point[1])/(image_position[0]-hit_point[0])
            output_angle = np.arctan(output_slope)/math.pi * 180
            ray_state = np.array([hit_point[0],hit_point[1],output_angle])
        self.state.append(ray_state)

    def plane_mirror(self,position):
        '''
        plane_mirror(self,position)
            Simulate the behavior of mirror.
            Append the ray state into self.state after reflected by the mirror.
        Parameters
        ------------
        position: float
            The x position of the mirror.
        '''
        slope = (self.theta/180)*math.pi
        vertical_dis = math.tan(slope)*(position-self.x)                
        ray_state = np.array([position+self.x,vertical_dis+self.z,180-self.theta])       

        self.state.append(ray_state)
        
    def plot_ray(self):
        '''
        plot_ray(self)
            Plot the ray path which is described in self.state.
        '''
        for i in range(len(self.state)):
            slope = self.state[i][2]/180 * math.pi
            if i < len(self.state)-1:             
                x = np.linspace(self.state[i][0],self.state[i+1][0])
                y = math.tan(slope)*x+self.state[i][1]
                x_main = np.linspace(0,max(x)*2)
                y_main = [0]*len(x_main)
                plt.plot(x_main,y_main,'k--')
                plt.plot(x,y)
            else:
                if self.state[i][2]<90:
                    x = np.linspace(self.state[i][0],self.state[i][0]+self.state[1][0])
                    y = math.tan(slope)*x+(self.state[i][1]-math.tan(slope)*self.state[i][0])
                else:
                    x = np.linspace(self.state[i][0]-5,self.state[i][0])
                    y = math.tan(slope)*x-2*math.tan(slope)*self.state[i][0]+self.state[0][1]
                plt.plot(x,y)
                
        for i in range(len(self.state)):
            if i == 0:
                pass
            else:
                x = self.state[i][0]
                y_axis = range(-7,7)
                x_axis = [x]*len(y_axis)
                plt.plot(x_axis,y_axis)
            
        plt.show()
        
    
if __name__=='__main__':
#    def test_mirror():
#        RT1 = RayTracing(0,2,5)
#        # Three parameters are x, z, angle
#        RT1.ray()
#        RT1.plane_mirror(15)
#        # Position of mirror
#        RT1.plot_ray()
#        print(RT1.state)
#    test_mirror()
    def test_lens():
        RT2 = RayTracing(0,2,0)
        # Three parameters are x, z, angle
        RT2.ray()
        RT2.lens(10,20) 
        # Two parameters are focal_length, position
        #There is a bug when the local_length equals to the position.   
        RT2.plot_ray()
        print(RT2.state)     
    test_lens()
        
            
        