import numpy as np
import math

# This class generate all binary maps
class build:
    def __init__(self,acc,angle,lo1,lo2,dia1,dia2,R,r,pi,po):
        
        angle=math.pi*angle/180  
        width=2*(lo2*math.sin(angle)+math.cos(angle)*dia2/2)
        length=lo1+lo2*math.cos(angle)+math.sin(angle)*dia2/2
        self.step=width/acc
        dim=[int(length/self.step)+5,acc]
        l1=lo1-0.5*(dia2/math.sin(angle)-dia1/math.tan(angle))-R*math.tan(angle/2)
        
        dia1/=self.step
        dia2/=self.step
        r/=self.step
        R/=self.step
        l1/=self.step
        
        
        #Initiate the maps
        self.binary_map_velocity=np.zeros((dim[0],dim[1]),dtype=float)
        self.binary_map=np.zeros((dim[0],dim[1]),dtype=float)
        self.binary_map_pressure=-1*np.ones((dim[0],dim[1]),dtype=float)

        
        k1=dim[1]/2-dia1/2
        xb=R *math.sin(angle)+l1
        yb=k1-(1-math.cos(angle))*R
        yc=yb+dia2*math.cos(angle)
        xc=xb+dia2*math.sin(angle)
        xm=r/math.sin(angle)-(dim[1]/2-yc)/math.tan(angle)+xc
        ym=k1+dia1/2
        xl=yb/math.tan(angle)+xb
     
        
        for x in range(dim[0]):
            for y in range(int(dim[1]/2)):
                
                if y>k1 or (y-yb)>(-math.tan(angle)*(x-xb)) or (x>l1 and (y-yb)>(math.tan(angle)*(x-xb)) and (x-l1)**2+(y-k1+R)**2>R**2) :
                    self.binary_map_velocity[x][y]=1
                    self.binary_map_pressure[x][y]=1
                    self.binary_map[x][y]=1  
                if y-yc>-math.tan(angle)*(x-xc) or y<(math.tan(math.pi/2-angle)*(x-xl)):
                    self.binary_map_velocity[x][y]=0
                    self.binary_map_pressure[x][y]=-1
                    self.binary_map[x][y]=0
                if y-yc<-math.tan(angle)*(x-xc) and (y-yb)>(-math.tan(angle)*(x-xb)) and y<(math.tan(math.pi/2-angle)*(x-xl)) :                   
                    self.binary_map[x][y]=2
                   
                if ((y-ym)>=(math.tan(math.pi/2-angle)*(x-xm)) and y-yc>=-math.tan(angle)*(x-xc) and (x-xm)**2+(y-ym)**2>=r**2 ) :
                    self.binary_map_velocity[x][y]=1
                    self.binary_map_pressure[x][y]=1
                    self.binary_map[x][y]=1                                   
                if (y==int(math.tan(math.pi/2-angle)*(x-xl)) or x==int(xl+y/math.tan(math.pi/2-angle))) and y-yc<-math.tan(angle)*(x-xc) :
                    self.binary_map_velocity[x][y]=-1
                    self.binary_map_pressure[x][y]=po
                    self.binary_map[x][y]=2                   
        for y in range(int(k1)+1,int(k1+dia1)-1):
            self.binary_map_velocity[0][y]=-1
            self.binary_map_pressure[0][y]=pi
            self.binary_map[0][y]=0                   
        for x in range(dim[0]):
            for y in range(int(dim[1]/2)):
                self.binary_map_velocity[x][int(dim[1]/2)+y]=self.binary_map_velocity[x][int(dim[1]/2)-y-1]
                self.binary_map_pressure[x][int(dim[1]/2)+y]=self.binary_map_pressure[x][int(dim[1]/2)-y-1]
                self.binary_map[x][int(dim[1]/2)+y]=self.binary_map[x][int(dim[1]/2)-y-1]                

    def get_binary_map_velocity(self) :
        return self.binary_map_velocity
    def get_binary_map_pressure(self) :
        return self.binary_map_pressure
    def get_binary_map(self) :
        return self.binary_map
    def get_step(self) :
        return self.step
    def get_dim(self) :
        return [len(self.binary_map),len(self.binary_map[0])]
