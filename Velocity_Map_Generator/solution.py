import numpy as np
from time import process_time

# This class computes the velocity map from the pressure map
class solution :
    def from_number(self,n) :
        return [int(n/self.dim[1]),n%self.dim[1]]
    def to_number(self,a,b) :
        n=0
        for z in range (len(self.exist)) :
            if self.exist[z][0]==a and self.exist[z][1]==b :
                n=1
                return z
        if n==0 :
            return "0"
          
    def __init__(self,matrix,matrix2,map,step):
        self.matrix=matrix
        self.map=map
        self.matrix2=matrix2
        self.step=step
        self.dim=[len(matrix),len(matrix[0])]
        
    def get_solution(self):
        print("\tStarting laplacien solution ...")   
        
        self.exist=[[]]
        b=[]
        for i in range (self.dim[0]*self.dim[1]) :
            if self.matrix[self.from_number(i)[0]][self.from_number(i)[1]]==1 :
                self.exist.append([self.from_number(i)[0],self.from_number(i)[1]])
                b.append(self.matrix2[self.from_number(i)[0]][self.from_number(i)[1]])
        self.exist.pop(0)
        
        #construction of matrix a
        print("\t\tStart building matrix ...")
        a = np.zeros((len(self.exist),len(self.exist)),dtype=float)
        all_progress=len(self.exist)
        pour=0
        for i in range (len(self.exist)) :
            
            if self.map[self.exist[i][0]][self.exist[i][1]]==1 :
                a[i][i]=-4
                for k in range(-1,2):
                    for l in range(-1,2):
                        if (l!=0 or k!=0)and k*l==0 :
                            
                            if self.matrix[self.exist[i][0]+k][self.exist[i][1]+l]==1 :
                                a[i][self.to_number(self.exist[i][0]+k,self.exist[i][1]+l)]=1
                            elif self.matrix[self.exist[i][0]+k][self.exist[i][1]+l]==-1 :
                                a[i][i]+=1
                            else :
                                b[i]-=self.matrix[self.exist[i][0]+k][self.exist[i][1]+l]/(2*self.step**2)

            if (100*i/all_progress-pour)>10 :
                pour=100*i/all_progress
                print("\t\t\t progress : {:.2f}%".format(pour))
        print("\t\tEnd building matrix")
        
        #Compute the solution of a*u=b
        print("\t\tStart inverting matrix ...")
        u = np.linalg.solve(a, b)
        result=self.matrix.copy()
        print("\t\tEnd inverting matrix")

        print("\t\tStart arrange matrix ...")
        all_progress=len(u)+self.dim[0]*self.dim[1]
        p=0
        pour=0
        for i in range (len(u)):
            u[i]*=(2*self.step**2)
            
            p+=1
            if (100*p/all_progress-pour)>10 :
                pour=100*p/all_progress
                print("\t\t\t progress : {:.2f}%".format(pour))

        for i in range (self.dim[0]):
            for j in range (self.dim[1]):
                if self.to_number(i,j)!="0" :
                    result[i][j]=u[self.to_number(i,j)]
                    
                p+=1
                if (100*p/all_progress-pour)>10 :
                    pour=100*p/all_progress
                    print("\t\t\t progress : {:.2f}%".format(pour))

        print("\t\tEnd arrange matrix")
        print("\tEnd laplacien solution")       
        return result
