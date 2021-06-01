# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 13:05:20 2019
Anderson Vinicius de Oliveira Rosa - 15159441
UFSC CTJ
"""

import pandas as pd #python package for sheet file processing
import beam_elements #functions for beam element FEM analisys

"""___input_data___""" #the names of nodes and elements are given by their line position (indexes)
nodes = pd.read_csv("nodes.csv") #csv file structured as: x_node_value, y_node_value
el = pd.read_csv("elements.csv") #csv file structured as: node_start, node_end, A, I, E
bc = pd.read_csv("bndcond.csv")  #csv file structured as: bc_type, node, direction, value

"""___item_a___"""
global_stiffness_array = beam_elements.global_stiffness_array(el,nodes)
disp = beam_elements.disp_solve(nodes,global_stiffness_array,bc)
load = beam_elements.load_solve(global_stiffness_array,disp)
stress = beam_elements.stress(nodes,el,disp,load,'Displacement result item a')

"""___item_b___"""
bc_b = bc; bc_b.value = bc_b.value*2;
disp_b = beam_elements.disp_solve(nodes,global_stiffness_array,bc_b)
load_b = beam_elements.load_solve(global_stiffness_array,disp_b)
stress = beam_elements.stress(nodes,el,disp_b,load_b,'Displacement result item b')
