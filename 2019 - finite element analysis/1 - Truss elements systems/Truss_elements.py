# -*- coding: utf-8 -*-
"""
Anderson Vinicius de Oliveira Rosa - 15159441
UFSC CTJ
2019-09-10
fuction list:
0 - #creates the global sitffness array for the system
1 - #solves the system for displacement
2 - #solves the system for reactions
3 - #calculates normal stress in each element and plot displacements
"""

import numpy as np                  #fundamental package for scientific computing in python
import matplotlib.pyplot as plt     #graph package for python

##0
def global_system_stiffness(node_coordinates,element_list):
    global_system_stiffness = np.zeros((len(node_coordinates)*2,len(node_coordinates)*2))
    for element_number in range (0,len(element_list)):
        #setting start and end nodes
        NS = int(element_list[element_number,0]) #start node
        NE = int(element_list[element_number,1]) #end node
        #calculate distance in x and y, and element length
        dx = node_coordinates[NE,0] - node_coordinates[NS,0] #element displacement in x
        dy = node_coordinates[NE,1] - node_coordinates[NS,1] #element displacement in y
        L = (dx**2+dy**2)**0.5 #element length
        #using dx and dy to find alpha
        if dx==0:
            if dy>0:
                alpha = np.pi/2
            else:
                alpha = 3*np.pi/2
        else:
            alpha = np.arctan(dy/dx)
        #creating empty transfer array
        TA = np.zeros((len(node_coordinates)*2,len(node_coordinates)*2))
        #area and young modulus
        EAL = element_list[element_number,3]*element_list[element_number,2]/L
        #inverting nodes if necessary
        if NE<NS:
            a=NE
            b=NS
            NE=b
            NS=a
        #transfer values
        C2 = np.cos(alpha)**2
        S2 = np.sin(alpha)**2
        SC = np.sin(alpha)*np.cos(alpha)
        #DOFs of first to first node
        TA[2*NS,2*NS] = C2
        TA[2*NS,2*NS+1] = SC
        TA[2*NS+1,2*NS] = SC
        TA[2*NS+1,2*NS+1] = S2
        #DOFs of first to second node
        TA[2*NS,2*NE] = -C2
        TA[2*NS,2*NE+1] = -SC
        TA[2*NS+1,2*NE] = -SC
        TA[2*NS+1,2*NE+1] = -S2
        #DOFs of second to second node
        TA[2*NE,2*NE] = C2
        TA[2*NE,2*NE+1] = SC
        TA[2*NE+1,2*NE] = SC
        TA[2*NE+1,2*NE+1] = S2
        #DOFs of second to first node
        TA[2*NE,2*NS] = -C2
        TA[2*NE,2*NS+1] = -SC
        TA[2*NE+1,2*NS] = -SC
        TA[2*NE+1,2*NS+1] = -S2
        #returning global element stiffness matrix
        global_system_stiffness = global_system_stiffness + EAL*TA
    return global_system_stiffness

##1
def disp_solve(node_coordinates,global_stiffness_array,bnd_cond_load,bnd_cond_disp):
    #creating 
    load_vector = np.zeros(len(node_coordinates)*2)
    #inserting load boundarie conditions
    for i in range(0,len(bnd_cond_load)):
        load_vector[int(bnd_cond_load[i,0])*2] = bnd_cond_load[i,2]
    #reducing global stiffness array and load vector
    for i in range(len(bnd_cond_disp)-1,-1,-1):
        a = bnd_cond_disp[i,0]*2 + bnd_cond_disp[i,1]
        load_vector = np.delete(load_vector, a, axis=0)
        global_stiffness_array = np.delete(global_stiffness_array, a, axis=0)
        global_stiffness_array = np.delete(global_stiffness_array, a, axis=1)
    #solving for displacements
    disp = np.linalg.solve(global_stiffness_array,load_vector)
    #insert fixed values in the displacement vector
    for i in range(0,len(bnd_cond_disp)):
        a = bnd_cond_disp[i,0]*2 + bnd_cond_disp[i,1]
        disp = np.insert(disp,a,0)
    return disp

##2
def load_solve(global_stiffness_array,system_disp):
    load = np.dot(global_stiffness_array,system_disp)
    return load

##3
def stress(node_coordinates,element_list,disp_vector,lV):
    #graph
    #
    plt.figure(num = None)                               #creates new figure, (figsize  = (6,4), dpi = 500)
    
    x = np.zeros(len(node_coordinates[:,0]))                        #creates vector with x coordinate values
    y = np.zeros(len(node_coordinates[:,1]))                        #creates vector with y coordinate values
    plt.plot(node_coordinates[:,0], node_coordinates[:,1], 'o')     #plot nodes
    df = 2e2                                                        #scale factor for plot
    #displaced nodes
    for i in range(0,len(y)):
        x[i] = node_coordinates[i,0] + df*disp_vector[2*i]
        y[i] = node_coordinates[i,1] + df*disp_vector[2*i+1]
    plt.plot(x, y, 'd')
    #
    #endgraph
    
                #plot displaced nodes
    
    
    #empty vector for stress values
    stress = np.zeros(len(element_list))
    for element_number in range (0,len(element_list)):
        #setting start and end nodes
        NS = int(element_list[element_number,0])                    #start node
        NE = int(element_list[element_number,1])                    #end node
        #calculate distance in x and y, and element length
        dx = node_coordinates[NE,0] - node_coordinates[NS,0]        #element displacement in x
        dy = node_coordinates[NE,1] - node_coordinates[NS,1]        #element displacement in y
        L = (dx**2+dy**2)**0.5 #element length
        #using dx and dy to find alpha
        if dx==0:
            if dy>0:
                alpha = np.pi/2
            else:
                alpha = 3*np.pi/2
        else:
            alpha = np.arctan(dy/dx)
        #area and young modulus
        EAL = element_list[element_number,3]*element_list[element_number,2]/L
        #transfer values
        C = np.cos(alpha)
        S = np.sin(alpha)
        TA = np.array([C,S])
        #displacements
        dU = disp_vector[2*NE] - disp_vector[2*NS]
        dV = disp_vector[2*NE+1] - disp_vector[2*NS+1]
        delta = np.array([dU,dV])
        #calculating stres
        stress[element_number] = np.dot(TA,delta)*EAL/element_list[element_number,2]
        
        #graph
        #
        for i in range (0,len(element_list)):
            if stress[element_number] >= 0:
                color = 'r--'
            if stress[element_number] < 0:
                color = 'b--'
            XS = node_coordinates[NS,0]
            XE = node_coordinates[NE,0]
            YS = node_coordinates[NS,1]
            YE = node_coordinates[NE,1]
            XSD = x[NS]
            XED = x[NE]
            YSD = y[NS]
            YED = y[NE]
            plt.plot([XE,XS],[YE,YS],"black")
            plt.plot([XED,XSD],[YED,YSD],color)
            plt.xlabel('x coordinate', style='italic')
            plt.ylabel('y coordinate', style='italic')
            plt.legend(['original node','displaced node','element','traction','compression'], loc=0)
            plt.title(lV)
        #
        #endplot
        
    plt.show()
    return stress