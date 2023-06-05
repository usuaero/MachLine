import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, './dev/helper_scripts')
from wedge_wing_generator import mesh_creator


nodes = [20, 29, 42, 62, 90]
grids = ['ultra_coarse', 'coarse', 'medium', 'fine', 'ultra_fine']


for i, (node, grid) in enumerate(zip(nodes, grids)):
    filename = 'studies/supersonic_love_delta_wing/meshes/delta_wing_clustered_mesh_{0}.vtk'.format(grid)
            
    mesh_creator(filename, R_T=0, x_cr=0.18, x_ct=0.18, b_2c=1.00652, node_ratio=0.18, cw_nodes=node, sw_nodes=node, cluster_cw=True, LE_sweep=44.85, mirror_xy=True)
