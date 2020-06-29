# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 09:23:07 2020

This file needs three packages:numpy, matplotlib.pyplt, math.

This .py file simulate the behavior of the prism. The geometry of the prism is equilateral triangle.
To trace the ray, we use the state of ray to describe. To trace the ray while refracted by the
prism, the snell's law is used. 

In the end we also plot the ray. This file did not take the full reflection into consideration,
if the incident angle is smaller  than -2 degree, the outcoming ray will disappear.

@author: Benchi Zhao
"""

import numpy as np
import matplotlib.pyplot as plt
import math

class PrismTracing:
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
            To aviod the bug, make sure the input value is greater than -2.
        self.n_air： float
            The refractive index of air.
        self.n_glass: float
            The refractive index of glass.
        self.state: list
            When ray interacting with the optical equipments, the ray state will change,
            all states are recorded in self.state.
        self.central_point: float
            Position of central point of the prism.
        self.side_length: float
            Length of each side of the prism.
        '''
        self.x = x
        self.z = z
        self.theta = theta
        self.n_air = 1.0
        self.n_glass = 1.5
        self.state = []
        self.central_point = 0
        self.side_length = 0
        
        
    def ray(self):
        '''
        ray(self)
            Append the initial ray state into the total ray state.
        '''
        ray_state = np.array([self.x,self.z,self.theta])
        self.state.append(ray_state)
        
        
    def prism(self,side_length,central_point):
        '''
        prism(self,side_length,central_point)
            Simulate the behavior of prism.
            Append the ray state into self.state after passing the prism.
            
        Parameters
        ------------
        side_length: float
            Length of each side of the prism.
        self.central_point: float
            Position of central point of the prism.
        '''
        self.central_point = central_point
        self.side_length = side_length
        # The ray incident into the prism
        incident_slope_1 = np.tan(np.deg2rad(self.state[-1][2]))
        L = np.array([[-incident_slope_1,1],[-math.sqrt(3),1]])
        R = np.array([self.state[-1][1]-incident_slope_1*self.state[-1][0],-math.sqrt(3)*self.central_point+self.side_length/math.sqrt(3)])
        result = np.linalg.solve(L,R)
        # Calculate the position of interacting point
        incident_angle_1 = 30 + self.state[-1][2]
        out_angle_1 = np.rad2deg(np.arcsin(self.n_air/self.n_glass * np.sin(np.deg2rad(incident_angle_1))))
        ray_state = np.array([result[0],result[1],out_angle_1-30])
        self.state.append(ray_state)
        # The ray come out from the prism
        incident_slope_2 = np.tan(np.deg2rad(self.state[-1][2]))
        L = np.array([[-incident_slope_2,1],[math.sqrt(3),1]])
        R = np.array([self.state[-1][1]-incident_slope_2*self.state[-1][0],math.sqrt(3)*self.central_point+self.side_length/math.sqrt(3)])
        result = np.linalg.solve(L,R)
        # Calculate the position of interacting point
        incident_angle_2 = 60- out_angle_1
        out_angle_2 = np.rad2deg(np.arcsin(self.n_glass/self.n_air * np.sin(np.deg2rad(incident_angle_2))))
        ray_state = np.array([result[0],result[1],30-out_angle_2])
        self.state.append(ray_state)
    
    def plot_ray(self):
        '''
        plot_ray(self)
            Plot the prism and the ray path which is described in self.state.
        '''
        # plot the prism
        x1 = np.linspace(-self.side_length/2+self.central_point,0+self.central_point)
        y1 = math.sqrt(3)*(x1-self.central_point)+ 2/math.sqrt(3)*self.side_length/2
        x2 = np.linspace(0+self.central_point,self.side_length/2+self.central_point)
        y2 = -math.sqrt(3)*(x2-self.central_point) + 2/math.sqrt(3)*self.side_length/2
        x3 = np.linspace(-self.side_length/2+self.central_point,self.side_length/2+self.central_point)
        y3 = [min(y1)]*len(x3)
        plt.plot(x1,y1,'k')
        plt.plot(x2,y2,'k')
        plt.plot(x3,y3,'k')
        
        # plot ray
        for i in range(len(self.state)):
            slope = np.tan(np.deg2rad(self.state[i][2]))
            if i < len(self.state)-1:             
                x = np.linspace(self.state[i][0],self.state[i+1][0])
                y = np.linspace(self.state[i][1],self.state[i+1][1],len(x))
                plt.plot(x,y)
            else:
                x = np.linspace(self.state[i][0],self.state[i][0]+self.state[1][0])
                y = slope*x+(self.state[i][1]-slope*self.state[i][0])
                plt.plot(x,y)
        plt.show()
        
if __name__=='__main__':
    def main():
        PT = PrismTracing(0,-1,10)
        # Three parameters are x, z, angle
        PT.ray()
        PT.prism(4,6)
        # Two parameters are side_length , central position
        PT.plot_ray()
        print(PT.state)
    main()