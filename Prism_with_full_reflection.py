# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 14:05:13 2020

@author: Benchi Zhao
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 09:23:07 2020

@author: Benchi Zhao
"""

import numpy as np
import matplotlib.pyplot as plt
import math

class PrismTracing:
    '''
    Three parameters to decribe the ray path
                     1. position of x
                     2. position of y
                     3. tilt angle (in degree)
    '''
    def __init__(self,x,z,theta):
        self.x = x
        self.z = z
        self.theta = theta
        self.n_air = 1.0
        self.n_glass = 1.5
        self.state = []
        self.central_point = 0
        self.side_length = 0
        self.full_reflection = False
        self.if_in = True
        
        
    def ray(self):
        ray_state = np.array([self.x,self.z,self.theta])
        self.state.append(ray_state)
        
        
    def prism(self,side_length,central_point):
        self.central_point = central_point
        self.side_length = side_length
        
        incident_slope_1 = np.tan(np.deg2rad(self.state[-1][2]))
        L = np.array([[-incident_slope_1,1],[-math.sqrt(3),1]])
        R = np.array([self.state[-1][1]-incident_slope_1*self.state[-1][0],-math.sqrt(3)*self.central_point+self.side_length/math.sqrt(3)])
        result = np.linalg.solve(L,R)
        incident_angle_1 = 30 + self.state[-1][2]
        out_angle_1 = np.rad2deg(np.arcsin(self.n_air/self.n_glass * np.sin(np.deg2rad(incident_angle_1))))
        ray_state = np.array([result[0],result[1],out_angle_1-30])
        self.state.append(ray_state)
        
        while self.full_reflection == True:
            incident_slope_2 = np.tan(np.deg2rad(self.state[-1][2]))
            L = np.array([[-incident_slope_2,1],[math.sqrt(3),1]])
            R = np.array([self.state[-1][1]-incident_slope_2*self.state[-1][0],math.sqrt(3)*self.central_point+self.side_length/math.sqrt(3)])
            result = np.linalg.solve(L,R)
            incident_angle_2 = 60- out_angle_1
            print(incident_angle_2)
            if incident_angle_2 < np.rad2deg(np.arcsin(self.n_air/self.n_glass)):    
                out_angle_2 = np.rad2deg(np.arcsin(self.n_glass/self.n_air * np.sin(np.deg2rad(incident_angle_2))))
                
            else:
                out_angle_2 = incident_angle_2+210
                self.full_reflection = True
            
            ray_state = np.array([result[0],result[1],30-out_angle_2])
            self.state.append(ray_state)
    
    def plot_ray(self):
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
        if self.full_reflection == False:
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
        elif self.full_reflection == True:
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
            pass
        plt.show()
        
if __name__=='__main__':
    def main():
        PT = PrismTracing(0,2,-15)
        PT.ray()
        PT.prism(5,5)
        PT.plot_ray()
        print(PT.state)
    main()