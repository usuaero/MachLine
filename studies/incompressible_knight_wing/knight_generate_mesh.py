import numpy as np
import machupX as mx
import matplotlib.pyplot as plt
import pypan as pp

import matplotlib
import json
import copy


if __name__ == "__main__":

    # Declare MachUpX input
    V = 59.5
    alpha = -8.0
    wing_input = {
        "weight" : 1.0,
        "airfoils" : {
            "G387FB" : {
                "type" : "linear",
                "geometry" : {
                    "outline_points" : "G387FB.csv"
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
                    "N" : 30
                },
                "CAD_options" : {
                    "round_wing_tip" : True,
                    "n_rounding_sections" : 6
                }
            }
        }
    }

    # Export mesh
    scene = mx.Scene()
    scene.add_aircraft("wing", wing_input, state={"velocity" : V})
    scene.export_stl(filename="knight_wing.stl", section_resolution=51)

    my_mesh = pp.Mesh(name="", mesh_file="knight_wing.stl", mesh_file_type="STL")
    my_mesh.export_vtk("knight_wing.vtk")