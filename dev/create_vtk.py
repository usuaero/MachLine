import json
import machupX as MX
import pypan as pp
import numpy as np


if __name__=="__main__":

    # Specify airplane
    airplane_dict = {
        "weight" : 50.0,
        "units" : "English",
        "airfoils" : {
            "NACA_0010" : {
                "geometry" : {
                    "NACA" : "0005",
                    "NACA_closed_te" : True
                }
            },
            "NACA_2410" : {
                "geometry" : {
                    "NACA" : "2410",
                    "NACA_closed_te" : True
                }
            }
        },
        "wings" : {
            "main_wing" : {
                "ID" : 1,
                "side" : "both",
                "is_main" : True,
                "semispan" : 2.0,
                "airfoil" : "NACA_0010",
                "chord" : [[0.0, 1.0],
                           [0.2, 1.0],
                           [1.0, 1.0]],
                "sweep" : [[0.0, 00.0],
                           [0.2, 00.0],
                           [0.2, 00.0],
                           [1.0, 00.0]],
                "grid" : {
                    "N" : 10,
                },
                "CAD_options" :{
                    "round_wing_tip" : True,
                    "round_wing_root" : False,
                    "n_rounding_sections" : 6
                }
            #},
            #"h_stab" : {
            #    "ID" : 2,
            #    "side" : "both",
            #    "is_main" : False,
            #    "connect_to" : {
            #        "ID" : 1,
            #        "location" : "root",
            #        "dx" : -4.0,
            #    },
            #    "semispan" : 1.0,
            #    "chord" : [[0.0, 0.8],
            #               [1.0, 0.5]],
            #    "sweep" : [[0.0, 25.0],
            #               [1.0, 25.0]],
            #    "grid" : {
            #        "N" : 20,
            #    },
            #    "CAD_options" :{
            #        "round_wing_tip" : True,
            #        "round_wing_root" : False,
            #        "n_rounding_sections" : 20
            #    }
            }
        }
    }

    # Load scene
    state = {
        "velocity" : 100.0
    }
    scene = MX.Scene()
    scene.add_aircraft("plane", airplane_dict, state=state)

    # Export stl
    stl_file = "dev/meshes/straight_wing_med_res.stl"
    scene.export_stl(filename=stl_file, section_resolution=41)

    # Load into PyPan and export vtk
    my_mesh = pp.Mesh(name="", mesh_file=stl_file, mesh_file_type="STL")
    my_mesh.export_vtk(stl_file.replace(".stl", ".vtk"))