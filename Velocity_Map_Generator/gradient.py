import numpy as np


class gradient :
    
    def __init__(self,matrix,map,step):
        self.matrix=matrix
        self.map=map
        self.step=step
        
        self.dim=[len(matrix),len(matrix[0])]
        
    def get_gradient(self):
        result_x=self.map.copy()
        result_y=self.map.copy()
        for x in range(self.dim[0]) :
            for y in range(self.dim[1]) :
                if self.map[x][y]==1 :
                    if self.matrix[x+1][y]!=-1 :
                        
                        result_x[x][y]=(self.matrix[x+1][y]-self.matrix[x][y])/self.step

                    else :
                        result_x[x][y]=0
                        
                    if self.matrix[x][y+1]!=-1 :
                        result_y[x][y]=(self.matrix[x][y+1]-self.matrix[x][y])/self.step

                    else :
                        result_y[x][y]=0
       

        return [result_x,result_y]
