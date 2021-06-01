# -*- coding: utf-8 -*-
"""
%Anderson Vinicius de Oliveira Rosa - 15159441
$spring-based systems
%UFSC CTJ
%August 2019
"""
import numpy as np

#create spring element stifness array
def Spring_element(element_stiffness,node_start,node_end,n_of_nodes):
    #creating generic system sized generic zeros array
    K = np.zeros([n_of_nodes,n_of_nodes])
    #setting array values
    K[node_start,node_start] = element_stiffness
    K[node_end,node_end] = element_stiffness
    K[node_start,node_end] = -element_stiffness
    K[node_end,node_start] = -element_stiffness
    return K

#create global stifness array
def Spring_global_stifness(springs,n_of_nodes):
    K_global = np.zeros([n_of_nodes,n_of_nodes])
    for i in range(0,len(springs)):
        #creating stifness array for element i
        K_temp = Spring_element(springs[i,0],int(springs[i,1]),int(springs[i,2]),n_of_nodes)
        #sum created stiffness array in the global array
        K_global = K_global + K_temp
    return K_global

#solve spring system
def Spring_system_displacement(K_global,Boundarie_conditions,fixed_node,n_of_nodes): 
    #creating vectors for force and displacement
    P = np.zeros(n_of_nodes)
    U = np.zeros(n_of_nodes)
    
    for i in range (0,len(Boundarie_conditions)):
        if Boundarie_conditions[i][2] == 'load':
            P[Boundarie_conditions[i][1]] = Boundarie_conditions[i][0]
        if Boundarie_conditions[i][2] == 'disp':
            U[Boundarie_conditions[i][1]] = Boundarie_conditions[i][0]
    
    #reducing global array
    P = np.delete(P, (fixed_node), axis=0)
    #U = np.delete(U, (fixed_node), axis=0)
    K_red = np.delete(K_global, (fixed_node), axis=0)
    K_red = np.delete(K_red, (fixed_node), axis=1)
    print(K_red)
    
    #calculating displacements
    U_zero = np.zeros(1)
    U = np.concatenate((U_zero,np.linalg.solve(K_red,P)),axis=0) #create displacement vector
    
    return U

#calculating internal reactions for each spring element
def Spring_internal_forces(springs,node_deflection):  
    #creating generic zeros vector
    f_int = np.zeros(len(springs))
    #calculating internal force in each spring
    for i in range(0,len(springs)):
        f_int[i] = springs[i,0] * (node_deflection[springs[i,2]] - node_deflection[springs[i,1]])
    return f_int
