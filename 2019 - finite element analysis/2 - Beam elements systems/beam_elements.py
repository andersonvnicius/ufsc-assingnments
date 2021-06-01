# -*- coding: utf-8 -*-
"""
Anderson Vinicius de Oliveira Rosa - 15159441
UFSC CTJ
oct-2019
"""

import numpy as np                  #fundamental package for scientific computing in python
import matplotlib.pyplot as plt     #graph package for python

##00
def global_stiffness_array(elements,nodes):  
    #element distances in global CS
    dx = np.array(nodes.x[elements.node_end]) - np.array(nodes.x[elements.node_start])
    dy = np.array(nodes.y[elements.node_end]) - np.array(nodes.y[elements.node_start])
    #element length
    l = np.sqrt(dx**2 + dy**2)
    #element stiffness for normal and bending loads
    a = elements.A * elements.E /l
    b = elements.E * elements.I /l**3
    #sin and cos of alpha
    S = dy/l
    C = dx/l
    #creating transfer array for each element
    T = np.zeros([len(elements),3*len(nodes),3*len(nodes)])
    for i in range (0,len(elements)):
        NS = elements.node_start[i]
        NE = elements.node_end[i]
        T[i,3*NS,3*NS] = T[i,3*NE,3*NE] = C[i]
        T[i,3*NS+1,3*NS+1] = T[i,3*NE+1,3*NE+1] = C[i]
        T[i,3*NS+2,3*NS+2] = T[i,3*NE+2,3*NE+2] = 1
        T[i,3*NS+1,3*NS] = T[i,3*NE+1,3*NE] = S[i]
        T[i,3*NS,3*NS+1] = T[i,3*NE,3*NE+1] = -1*S[i]
    #creating stiffness array for each element
    K = np.zeros([len(elements),3*len(nodes),3*len(nodes)])
    K_gcs = np.zeros([len(elements),3*len(nodes),3*len(nodes)])
    for i in range (0,len(elements)):
        NS = elements.node_start[i]
        NE = elements.node_end[i]
        #a values
        K[i,3*NS,3*NS] = K[i,3*NE,3*NE] = a[i]
        K[i,3*NS,3*NE] = K[i,3*NE,3*NS] = -1*a[i]
        #main values with b
        K[i,3*NS+1,3*NS+1] = K[i,3*NE+1,3*NE+1] = 12*b[i]
        K[i,3*NS+2,3*NS+2] = K[i,3*NE+2,3*NE+2] = 4*b[i]*l[i]**2
        K[i,3*NS+1,3*NS+2] = 6*b[i]*l[i]
        K[i,3*NE+1,3*NE+2] = -6*b[i]*l[i]
        K[i,3*NS+2,3*NS+1] = 6*b[i]*l[i]
        K[i,3*NE+2,3*NE+1] = -6*b[i]*l[i]
        #secondary values with b
        K[i,3*NS+1,3*NE+1] = -12*b[i]
        K[i,3*NE+1,3*NS+1] = -12*b[i]
        K[i,3*NS+2,3*NE+2] = 2*b[i]*l[i]**2
        K[i,3*NE+2,3*NS+2] = 2*b[i]*l[i]**2
        K[i,3*NS+1,3*NE+2] = 6*b[i]*l[i]
        K[i,3*NE+1,3*NS+2] = -6*b[i]*l[i]
        K[i,3*NS+2,3*NE+1] = -6*b[i]*l[i]
        K[i,3*NE+2,3*NS+1] = 6*b[i]*l[i]
        #creating stiffness array in global cs
        K_gcs[i] = np.dot(np.dot(T[i].transpose(),K[i]),T[i])
    #Global system stiffness array
    global_system_stiffness = np.sum(K_gcs,axis=0)
    return global_system_stiffness

##1
def disp_solve(nodes,global_stiffness_array,bnc_cond):
    #creating vector for load values
    load_vector = np.zeros(len(nodes)*3)
    #inserting load boundarie conditions
    for i in range(0,len(bnc_cond)):
        if bnc_cond.bc_type[i] == 'load':
            load_vector[3*bnc_cond.node[i] + bnc_cond.direction[i]] = bnc_cond.value[i]
    #reducing global stiffness array and load vector
    for i in range(len(bnc_cond)-1,-1,-1):
        if bnc_cond.bc_type[i] == 'disp':
            a = bnc_cond.node[i]*3 + bnc_cond.direction[i]
            load_vector = np.delete(load_vector, a, axis=0)
            global_stiffness_array = np.delete(global_stiffness_array, a, axis=0)
            global_stiffness_array = np.delete(global_stiffness_array, a, axis=1)
    #solving displacements
    disp = np.linalg.solve(global_stiffness_array,load_vector)
    #insert fixed values in the displacement vector
    for i in range(0,len(bnc_cond)):
        if bnc_cond.bc_type[i] == 'disp':
            a = bnc_cond.node[i]*3 + bnc_cond.direction[i]
            disp = np.insert(disp,a,0)
    return disp

##2
def load_solve(global_stiffness_array,system_disp):
    load = np.dot(global_stiffness_array,system_disp)
    return load

##3
def stress(nodes,el,disp,load,Title):
    #graph parameters
    plt.figure(num=None)
    
    #plot original nodes
    plt.plot(nodes.x, nodes.y, 'o')
    
    #plot diaplaced nodes
    df = 1e2 #scale factor for plot
    x=np.zeros(len(nodes));y=np.zeros(len(nodes))
    for i in range(0,len(y)):
        x[i] = nodes.x[i] + df*disp[3*i]
        y[i] = nodes.y[i] + df*disp[3*i+1]
    plt.plot(x, y, 'o')
    
    for element_number in range (0,len(el)):
        #setting start and end nodes
        NS = el.node_start[element_number] #start node
        NE = el.node_end[element_number]   #end node
        XS = nodes.x[NS]
        XE = nodes.x[NE]
        YS = nodes.y[NS]
        YE = nodes.y[NE]
        XSD = x[NS]
        XED = x[NE]
        YSD = y[NS]
        YED = y[NE]
        #plot original elements
        plt.plot([XE,XS],[YE,YS],"black")
        #plot displaced elements
        plt.plot([XED,XSD],[YED,YSD],"k--")
        plt.xlabel('x coordinate [in]', style='italic')
        plt.ylabel('y coordinate [in]', style='italic')
        plt.legend(['original node','displaced node','element','disp element'], loc=6)
        plt.title(Title)
    plt.show()
    
    #element distances in global CS
    dx = np.array(nodes.x[el.node_end]) - np.array(nodes.x[el.node_start])
    dy = np.array(nodes.y[el.node_end]) - np.array(nodes.y[el.node_start])
    #element length
    l = np.sqrt(dx**2 + dy**2)
    #element stiffness for normal and bending loads
    eal = el.A * el.E /l
    #sin and cos of alpha
    S = dy/l
    C = dx/l
    #displacements in local CS
    dU = disp[3*el.node_end] - disp[3*el.node_start]
    dV = disp[3*el.node_end+1] - disp[3*el.node_start+1]
    delta = dU*C + dV*S
    #calculate stress for each element
    stress = delta*eal
    
    #plot shear and normal diagrams
    #sorting loads
    loads = np.array([load[3*el.node_start],load[3*el.node_end],load[3*el.node_start+1],load[3*el.node_end+1]])
    #loads in global CS
    load_global = np.array([load[3*el.node_start],load[3*el.node_end],load[3*el.node_start+1],load[3*el.node_end+1]])
    #loads in local CS
    load_N_start = load_global[0,:]*C + load_global[2,:]*S
    load_N_end = load_global[1,:]*C + load_global[3,:]*S
    load_V_start = load_global[0,:]*S + load_global[2,:]*C
    load_V_end = load_global[1,:]*S + load_global[3,:]*C
    for i in range(0,len(loads[0,:])):
        #plot diagrams for each element
        plt.figure()
        plt.title("%s %s" % ('Element', i+1))
        plt.plot([0,l[i]],[load_N_start[i],load_N_end[i]],'blue')
        plt.plot([0,l[i]],[load_V_start[i],load_V_end[i]],'red')
        plt.ylabel('Load [lb]', style='italic')
        plt.xlabel('length [in]', style='italic')
        plt.legend(['Normal load','Shear load'], loc=0)
    
    return stress