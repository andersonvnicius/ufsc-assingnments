# -*- coding: utf-8 -*-
"""
Anderson Vinicius de Oliveira Rosa - 15159441
P3 introduction to FEM
UFSC CTJ
2019-09-10
"""

import numpy as np         #fundamental package for scientific computing in python
import Truss_elements      #user-created package for truss-element system solving

#input
AA = 0.004; EE = 70E9; #area[m2] and young modulus[N/m2]
node_coordinates = np.array(([0,0],[2.4,1],[2.4,0],[4.8,1]))
element_list = np.array(([0,2,AA,EE],[0,1,AA,EE],[2,1,AA,EE],[1,3,AA,EE],[2,3,AA,EE],))
bnd_cond_load = np.array(([3,0,30E3],))             #([node,direction,load[N]],) 0: X, 1: Y
bnd_cond_disp = np.array(([0,0,],[0,1],[3,1],))     #([node,direction],) 0: X, 1: Y

#item a
global_stiffness_array = Truss_elements.global_system_stiffness(node_coordinates,element_list)
displacements = Truss_elements.disp_solve(node_coordinates, global_stiffness_array, bnd_cond_load,bnd_cond_disp)
reactions = Truss_elements.load_solve(global_stiffness_array,displacements)
element_stress = Truss_elements.stress(node_coordinates,element_list,displacements,'P = 30kN, displaced nodes out of scale')
#end item a

#item b
bnd_cond_load_b = np.array(([3,0,300E3],))
displacements_b = Truss_elements.disp_solve(node_coordinates, global_stiffness_array, bnd_cond_load_b,bnd_cond_disp)
reactions_b = Truss_elements.load_solve(global_stiffness_array,displacements_b)
element_stress_b = Truss_elements.stress(node_coordinates,element_list,displacements_b,'P = 300kN, displaced nodes out of scale')
#end item b