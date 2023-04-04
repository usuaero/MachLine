import numpy as np
import machupX as mx
import matplotlib.pyplot as plt
from dev.helper_scripts.stl_to_vtk import *


def gen_knight_wing_mesh(N_spanwise, N_chordwise, N_round, filename):
    """Generates an stl of Knight's straight wing with the given parameters."""

    # Declare MachUpX input
    V = 59.5
    alpha = -8.0
    wing_input = {
        "weight" : 1.0,
        "airfoils" : {
            "G387FB" : {
                "type" : "linear",
                "geometry" : {
                    "outline_points" : "studies/incompressible_knight_wing/data/G387FB.csv"
                }
            }
        },
        "wings" : {
            "main_wing" : {
                "ID" : 1,
                "side" : "right",
                "is_main" : True,
                "semispan" : [15.0, "in"],
                "chord" : [5.0, "in"],
                "airfoil" : "G387FB",
                "grid" : {
                    "N" : N_spanwise
                },
                "CAD_options" : {
                    "round_wing_tip" : True,
                    "n_rounding_sections" : N_round
                }
            }
        }
    }

    # Export mesh
    scene = mx.Scene()
    scene.add_aircraft("wing", wing_input, state={"velocity" : V})
    scene.export_stl(filename=filename, section_resolution=N_chordwise)


if __name__ == "__main__":

    densities = ['coarse', 'medium', 'fine']
    N_spanwises = [15, 30, 60]
    N_chordwises = [21, 41, 81]
    N_rounds = [4, 8, 16]

    for i, density in enumerate(densities):

        stl_file = "studies/incompressible_knight_wing/meshes/knight_wing_{0}.stl".format(density)
        gen_knight_wing_mesh(N_spanwises[i], N_chordwises[i], N_rounds[i], stl_file)
        convert_stl_to_vtk(stl_file)