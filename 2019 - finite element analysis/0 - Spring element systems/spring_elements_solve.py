# -*- coding: utf-8 -*-
"""
%Anderson Vinicius de Oliveira Rosa - 15159441
%T2 Computational solid mechanics
%UFSC CTJ
%2019-08-24
"""

#import matplotlib.pyplot as plt
import numpy as np
import spring_elements #import elements

n_of_nodes = 5  #total number of nodes
fix_node = 0    #fixed node number

#stifness values of the springs [N/mm]
spring_k = [200, 300, 250, 500, 350]

#creating a list with element stifness, start node and end node for each spring
springs = np.array([[spring_k[0],0,1],  #spring k1 between node 0 and 1
                    [spring_k[0],0,1],  #spring k1 between node 0 and 1
                    [spring_k[1],1,2],  #spring k2 between node 1 and 2
                    [spring_k[2],2,3],  #spring k3 between node 2 and 3
                    [spring_k[2],2,3],  #spring k3 between node 2 and 3
                    [spring_k[3],0,3],  #spring k4 between node 0 and 3
                    [spring_k[4],3,4]]) #spring k5 between node 3 and 4

#creating the global stifness array
K_glb = spring_elements.Spring_global_stifness(springs,n_of_nodes) #global stiffness array

#creating a list with boundarie conditions (fixed node value equal 0 automatically aplied in fix_node position)
Bnd_cond = ([10E3,4,'load'],)      #Load applied in node 4

#solving the system for displacement values (item a)
U = spring_elements.Spring_system_displacement(K_glb,Bnd_cond,fix_node,n_of_nodes)

#finding the internal forces for each spring element
f_int = spring_elements.Spring_internal_forces(springs,U)

#Bnd_cond = ([10E3,4,'load'],[40E3,2,'load'])

#solving the system for displacement values (item b)
#U_b = spring_elements.Spring_system_displacement(K_glb,Bnd_cond,fix_node,n_of_nodes)
