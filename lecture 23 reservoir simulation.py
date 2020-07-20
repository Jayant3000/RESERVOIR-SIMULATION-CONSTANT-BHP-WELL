# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 19:57:28 2020

@author: Dell
"""


import numpy as np
import matplotlib.pyplot as plt 
import math

print(   "\t ########################################################################### \t"  )
print(   "\t #############                Lecture 23                                      ######## \t"  )
print(   "\t ############# Reservoir Simulation implicit solution to 1D            ############### \t"  )
print(   "\t ######## dirichlet(constant pressure) boundary consition on right side   ############ \t"  )
print(   "\t ############ left side zero flow rate Neumann boundary condition      ############### \t"  )
print(   "\t ##################           injection in first block(flow rate id given)  ############### \t"  )
print(   "\t ##################       production in third block  (prerssure given)   ############### \t"  )
print(   "\t ########################################################################### \t"  )



L = 10000
print("\n\tThe lenght of the reservoir is ", str(L))
number_nodes = 4 
print("\n\tBlock node in the reservoir is ", str(number_nodes))
P0 = 1000
print("\n\tThe intial pressure of the reservoir is ", str(P0))
P_right = 2000
print("\n\tThe pressure at the right boundary of the reservoir is in psi", str(P_right))

#P_left = 1000
#print("\n\tThe pressure at the left boundary of the reservoir is in psi", str(P_left))

dx = L/number_nodes
print("\n\tThe lenght of Blocks in the reservoir is ", str(dx))

porosity = 0.2
print("\n\tthe porosity value of the reservoir is ", str(porosity))

permeability  = 50
print("\n\tthe permeability(md) value of the reservoir is ", str(permeability))

viscosity = 1
print("\n\tthe viscosity value is ", str(viscosity))

area = 200000
print("\n\tCross sectional area of the reservoir ", str(area))

compressibility = 1*10**(-6)
print("\n\tcompressibility of the reservoir is ", str(compressibility))

Bw = 1
print("\n\t water formation volume factor is " +str(Bw)+ " rb/stb" )

radius = 0.25
print("\n\t production and injection well radius " + str(radius) + " ft" )

thickness = 250
print("\n\t reservoir thickness " + str(thickness) + " ft" )

skin = 0 
print("\n\t skin " + str(skin) + " . " )

production_thickness = (area / dx)
print("\n\t production_thickness " + str(production_thickness) + " ft" )

equivalent_radius = dx*0.2
print("\n\t equivalent_radius " + str(equivalent_radius) + " ft" )

productivity_index =  (( 2*(math.pi)*permeability*production_thickness) / (viscosity*Bw*math.log(equivalent_radius/radius) ))*6.33*10**(-3)                     

print( " productivity index at both well " + str(productivity_index) + "md-ft/cp")

production_pressure_3 = 800
print( " production pressure at 3 " + str(production_pressure_3) + "psi")

###############################################################################
###############################################################################
#################### final time for simulation is  ############################

t_final = 3
print("\n\t the reservoir simulation time in days is  " +  str(t_final) + "days")

#################### time step  ###############################################

dt_increment = 1
print("\n\t the reservoir simulation incremental time step in days is "+ str(dt_increment)+ "day")

###############################################################################
###############################################################################
###############################################################################
############            pressure and  boundary condition          #############

pressure_previous = np.ones([number_nodes,1])*P0
print("\n############## pressure distribution ################\n")
print("pressure distribution at day 0 is\n", str(pressure_previous))


###############################################################################
###############################################################################


#production_rate = productivity_index(pressure_previous[2] -  )


productivity_index_matrix = np.zeros([number_nodes , number_nodes])
productivity_index_matrix[2][2] =  productivity_index
#for i in range(1,number_nodes,1):
#    #print("i value" , i)
#    productivity_index_matrix[i][i] = 2
#print("\n transmisibility_matrix is\n", str(productivity_index_matrix))

print("\n############## productivity_index_matrix ################\n")
print(" productivity_index_matrix is\n ", str(productivity_index_matrix))


###############################################################################
###############################################################################
print("\n############## transmisibility_matrix ################\n")

transmisibility_matrix = np.zeros([number_nodes , number_nodes]) 
print("\n transmisibility_matrix is\n", str(transmisibility_matrix))


T = ((permeability*area )/(viscosity*Bw*dx ))*6.33*10**(-3)

print("\n transmissibility = " + str(T) )

for i in range(1,number_nodes,1):
    #print("i value" , i)
    transmisibility_matrix[i][i] = 2*T
    transmisibility_matrix[i][i-1] = - T
    transmisibility_matrix[i-1][i] = - T
    transmisibility_matrix[0][0]=  T
    transmisibility_matrix[number_nodes-1][number_nodes-1]= 3*T
print("\n transmisibility_matrix is\n", str(transmisibility_matrix))

    
###############################################################################

print("\n############## B_matrix ################\n")
    
B_matrix = np.zeros([number_nodes , number_nodes])
print("\n B_matrix is\n", str(B_matrix))

B = porosity*compressibility*area*dx

for i in range(1,number_nodes,1):
    #print("i value" , i)
    B_matrix[i][i] = B
    B_matrix[0][0]= B
    B_matrix[number_nodes-1][number_nodes-1]= B
print("\n B_matrix is\n", str(B_matrix))
###############################################################################


productivity_index_matrix_pressure = np.zeros([number_nodes,1])
productivity_index_matrix_pressure[2][0] =  productivity_index*production_pressure_3
productivity_index_matrix_pressure[3][0] =  2*T*P_right

print("\n############## productivity_index_matrix ################\n")
print(" productivity_index_matrix_pressure is\n ", str(productivity_index_matrix_pressure))


###############################################################################
print("\n############## Q_matrix ################\n")      
Q_matrix = np.zeros([number_nodes , 1])
print("\n Q_matrix is\n", str(Q_matrix))

#*6.33*10**(-3)*transmissibility_2[0]
Q_matrix[0][0]= 1000
Q_matrix_index = np.zeros([number_nodes , 1])
for i in range(0,number_nodes,1):
    #print(i)
    print(Q_matrix[i][0])
    Q_matrix_index[i][0] = Q_matrix[i][0] + productivity_index_matrix_pressure[i][0]    

#Q_matrix[number_nodes-1][number_nodes-1]= B
print("\n Q_matrix is\n", str(Q_matrix))
print("\n Q_matrix_index is\n", str(Q_matrix_index))


###############################################################################
print("\n############## N_plus_1_matrix ################\n")

inverse_dt =  ( 1 / dt_increment )

N_plus_1_matrix = ( transmisibility_matrix +  productivity_index_matrix + inverse_dt*B_matrix )

print("\n N_plus_1_matrix is\n", str(N_plus_1_matrix))

###############################################################################
print("\n############## N__matrix ################\n")

inverse_dt =  ( 1 / dt_increment )

N__matrix = (inverse_dt*B_matrix)

print("\n N__matrix is\n", str(N__matrix))

print("\n",pressure_previous)
for k in range(0, t_final, 1):
    
    #print("\n time step value",k)
#    boundary_condition_array[0][0] = 2*neta*P_left
    pressure_previous = np.dot(np.linalg.inv(N_plus_1_matrix), (np.dot(N__matrix , pressure_previous) + Q_matrix_index))
    print("\n",pressure_previous)

print("\nfinal pressure \n", pressure_previous )

print(N__matrix.shape)
print(pressure_previous.shape)






























