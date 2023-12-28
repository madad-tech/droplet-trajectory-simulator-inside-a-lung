import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from map_builder import *
from solution import *
from gradient import *
from numpy import savetxt,loadtxt
import argparse

#Argument parser
parser = argparse.ArgumentParser(description='Generator of Velocity Map')
parser.add_argument('-input', dest='input_file',default="input.csv", required=False, help='CSV that contains system parameters')
parser.add_argument('-display_binary_maps', dest='display_binary_maps',default=False, required=False, help='Plot the generated binary maps')
parser.add_argument('-display_pressure_map', dest='display_pressure_map',default=False, required=False, help='Plot the generated pressure map')
parser.add_argument('-display_velocity_map', dest='display_velocity_map',default=True, required=False, help='Plot the generated velocity map')
parser.add_argument('-save', dest='save',required=False, help='Save the velocity map (Please provide the name of the destination folder)')
args = parser.parse_args()

#Load input data
input=loadtxt('./'+args.input_file, delimiter=',', usecols=(1))
d1,d2,l1,l2,Rc,rc,teta,visct,pi,po,precision=input

#Generating Binary Map
print("Start building binary maps ...")
t_start = process_time()
builder=build(int(precision),teta,l1,l2,d1,d2,Rc,rc,pi,po)
t_stop = process_time()
print("End building binary maps in {} s\n".format(t_stop-t_start))

step=builder.get_step()
dim=builder.get_dim()

#Binary Map
binary_map=builder.get_binary_map()

#Velocity Binary Map
binary_map_velocity=builder.get_binary_map_velocity()

#Pressure Binary Map
binary_map_pressure=builder.get_binary_map_pressure()


#Display Binary Maps
if args.display_binary_maps :
    plt.subplot(211)
    plt.imshow(binary_map_velocity.transpose())
    plt.title('Binary Map Velocity')
    plt.subplot(223)
    plt.imshow(binary_map_pressure.transpose())
    plt.title('Binary Map Pressure')
    plt.subplot(224)
    plt.imshow(binary_map.transpose())
    plt.title('Binary Map')
    plt.suptitle('Binary Maps', fontsize=16)  # Add a big title to the entire figure
    plt.show()


#Computing Pressure Map
#Solve equation : laplacien pression = 0
print("Start computing pressure map ....")
t_start = process_time()
pressure_map=solution(binary_map_pressure,np.zeros((dim[0],dim[1]),dtype=float),binary_map,step).get_solution()
t_stop = process_time()
print("End computing pressure map in {} s\n".format(t_stop-t_start))

#Display Pressure Map
if args.display_pressure_map :
    plt.imshow(pressure_map.transpose())
    plt.colorbar(label='Pressure (Pa)')
    plt.title('Pressure Map')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')
    plt.show()


#Computing Pressure Gradient Map
print("Start computing pressure gradient map ....")
t_start = process_time()
binary_map_grad_pression=gradient(pressure_map,binary_map,step).get_gradient()
t_stop = process_time()
print("End computing pressure gradient map in {} s\n".format(t_stop-t_start))


#Computing Velocity Map
print("Start computing velocity map ....")
t_start = process_time()
for i in range(len(binary_map_grad_pression)):
    for j in range (len(binary_map_grad_pression[0])):
        for k in range(len(binary_map_grad_pression[0][0])):
            binary_map_grad_pression[i][j][k]/=visct
print("Computing X-Velocity Component")
velocity_map_x=solution(binary_map_velocity,binary_map_grad_pression[0],binary_map,step).get_solution()
print("Computing Y-Velocity Component")
velocity_map_y=solution(binary_map_velocity,binary_map_grad_pression[1],binary_map,step).get_solution()
velocity_map=[velocity_map_x,velocity_map_y]

for i in range (len(velocity_map)) :
    for j in range (dim[0]):
        for k in range (dim[1]):
            if binary_map[j][k]!=1 :
                velocity_map[i][j][k]=0
        
t_stop = process_time()
print("End computing velocity map in {} s\n".format(t_stop-t_start))


#Diplay velocity map
if args.display_velocity_map :

    x1=np.arange(0,dim[0],1)
    y1=np.arange(0,dim[1],1)
    X1,Y1=np.meshgrid(x1,y1)
    points=np.transpose([np.repeat(x1, len(y1)),np.tile(y1, len(x1))])
    x2=np.arange(0,dim[0],0.1)
    y2=np.arange(0,dim[1],0.1)
    X2,Y2=np.meshgrid(x2,y2)

    grid_x=griddata(points, velocity_map[0].flatten(),(X2,Y2), method='cubic').transpose();
    grid_y=griddata(points, velocity_map[1].flatten(),(X2,Y2), method='cubic').transpose();
   
    n=4
    scale_factor=0.001
    lkey=50
    
    plt.figure(figsize=(12, 4))
    
    plt.subplot(121)
    # Plot quiver plot
    q_l = plt.quiver(X1[::n, ::n], Y1[::n, ::n], velocity_map[0][::n, ::n].transpose()*scale_factor, velocity_map[1][::n, ::n].transpose()*scale_factor, scale=1)
    plt.title('Velocity Map with Vector representation')
    plt.quiverkey(q_l, X=0.3, Y=1.1, U=lkey*scale_factor, label='Quiver key, length = ' + str(lkey) + ' m/s', labelpos='E')
    plt.axis('equal')
    plt.xlabel('X-axis')
    plt.ylabel('Y-axis')

    plt.subplot(122)
    # Plot norm plot
    plt.imshow(np.sqrt(grid_y.T*grid_y.T + grid_x.T*grid_x.T))
    plt.title('Velocity Norm Map')
    plt.axis('equal')
    plt.show()

#Save the results
if args.save :
    savetxt(args.save+'/binary_map.csv', binary_map, delimiter=',')
    savetxt(args.save+'/datax.csv', velocity_map[0], delimiter=',')
    savetxt(args.save+'/datay.csv', velocity_map[1], delimiter=',')
    savetxt(args.save+'/metadata.csv', [d1,d2,Rc,rc,teta,visct,pi,po,step], delimiter=',')

if args.save!=None :
    print("Results are successfully saved in '{}' folder".format(args.save))