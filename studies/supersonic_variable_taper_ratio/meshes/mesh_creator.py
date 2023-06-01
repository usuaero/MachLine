import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, './dev/helper_scripts')
from wedge_wing_generator import mesh_creator


nodes = [20,40,60,80]
grids = ['ultra_coarse', 'coarse', 'medium', 'fine']
R_ts = [0,0.25,0.5,0.75,1]
clustering = [True,False]

for i, (node, grid) in enumerate(zip(nodes, grids)):
    for j, R_t in enumerate(R_ts):
        LE_sweep = 180 * np.arctan2(1-R_t,1)/np.pi
        for k, clustered in enumerate(clustering):
            if clustered :
                filename = 'studies/supersonic_variable_taper_ratio/meshes/wing_R_T_{0:.0f}_{1}_clustered.vtk'.format(R_t*100, grid)
            else:
                filename = 'studies/supersonic_variable_taper_ratio/meshes/wing_R_T_{0:.0f}_{1}_non-clustered.vtk'.format(R_t*100, grid)
            

            mesh_creator(filename, R_T=R_t, cw_nodes=node, sw_nodes=node, cluster_cw=clustered, LE_sweep=LE_sweep, mirror_xy=False)
