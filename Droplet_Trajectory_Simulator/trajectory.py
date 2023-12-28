import numpy as np
import math


# This class instantiates an object capable of calculating the trajectory of a droplet
class trajectory :
    
    def __init__(self,binary_map,velocity_map_x,velocity_map_y,visct,step):
        self.velocity_map_x=velocity_map_x
        self.velocity_map_y=velocity_map_y
        self.binary_map=binary_map
        self.visct=visct
        self.step=step
        self.time_step=self.step/max(np.max(self.velocity_map_x),np.max(self.velocity_map_y))
    def get_map(self):
        return self.binary_map
    def get_time_step(self):
        return self.time_step
        
    # Initial velocity can be the same as the fluid if it is set to 's'
    def compute(self,droplet_density,droplet_radius,coord_x,coord_y,v_x,v_y):
        coord_x=[ coord_x]
        coord_y=[ coord_y]

        if v_x=='s' :
            v_x=[self.velocity_map_x[int(coord_x[0]/self.step)][int(coord_y[0]/self.step)]]
        else:
            v_x=[v_x]
            
        if v_y=='s' :
            v_y=[self.velocity_map_y[int(coord_x[0]/self.step)][int(coord_y[0]/self.step)]]
        else:
            v_y=[v_y]
        
        m=droplet_density*4/3*math.pi*droplet_radius**3
        
        cnt = (self.visct*(np.pi)*6*droplet_radius)/m
        
        i=0
        while self.binary_map[int(coord_x[i]/self.step)][int(coord_y[i]/self.step)]==1 :


            # Compute horizontal component of the velocity 
            v_xi=v_x[i]*math.exp(-cnt*self.time_step)+self.velocity_map_x[int(coord_x[i]/self.step)][int(coord_y[i]/self.step)]*(1-math.exp(-cnt*self.time_step))
            # Compute vertical component of the velocity 
            v_yi=v_y[i]*math.exp(-cnt*self.time_step)+self.velocity_map_y[int(coord_x[i]/self.step)][int(coord_y[i]/self.step)]*(1-math.exp(-cnt*self.time_step))
            
            cd_x=self.time_step*v_xi+coord_x[i]
            cd_y=self.time_step*v_yi+coord_y[i]

            # Conditioning the stop event
            if cd_x>=len(self.binary_map) or (cd_x/self.step)<0 or (cd_y/self.step)>=len(self.binary_map[0]) or (cd_y/self.step)<0 or i>5000 : 
                break
                
            v_x.append(v_xi)
            v_y.append(v_yi)
            coord_x.append(cd_x)
            coord_y.append(cd_y)
            i+=1

        # State variable refer to the deposition of the droplet ; 0 not deposited - 1 deposited
        state=0
        if self.binary_map[int(coord_x[i]/self.step)][int(coord_y[i]/self.step)]==2 :
            state=1
        return ([coord_x,coord_y,state])  
