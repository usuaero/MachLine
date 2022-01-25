import json
import subprocess as sp
import flow54 as fl
import numpy as np


if __name__=="__main__":

    # Parameters
    M = 2.0
    c_inf = 500.0
    gamma = 1.4
    alpha = 2.0

    # Declare MachLine input
    input_dict = {
        "flow": {
            "freestream_velocity": [M*c_inf*np.cos(np.radians(alpha)), 0.0, M*c_inf*np.cos(np.radians(alpha))],
            "gamma" : gamma,
            "freestream_mach_number" : M
        },
        "geometry": {
            "file": "dev/meshes/diamond_5_deg_coarse.vtk",
            "mirror_about": "xz",
            "spanwise_axis" : "+y",
            "wake_model": {
                "append_wake" : False,
            },
            "reference": {
                "area": 1.0
            }
        },
        "solver": {
            "formulation": "morino",
            "control_point_offset": 1.1e-05,
            "influence_calculations" : "epton-magnus"
        },
        "output": {
            "body_file": "dev/results/diamond_5_deg_coarse.vtk"
        }
    }

    # Dump
    input_file = "dev/diamond_input.json"
    with open(input_file, 'w') as input_handle:
        json.dump(input_dict, input_handle)

    # Run
    sp.run(["./machline.exe", input_file])
    
    # Run shock-expansion comparison
    rho = 1.225
    p_inf = 1.0e5
    airfoil = fl.DiamondAirfoil(5.0, 1.0)
    airfoil.set_state(M, alpha, gamma, p_inf, 300.0, c_inf, rho, 1.0e-5, 800.0, 0.7)
    p2, p3, p4, p5 = airfoil.get_pressures()
    x = 0.5*gamma*M**2
    Cp2 = (p2/p_inf-1.0)/x
    Cp3 = (p3/p_inf-1.0)/x
    Cp4 = (p4/p_inf-1.0)/x
    Cp5 = (p5/p_inf-1.0)/x

    # Print table
    print("{0:<20}{1:<20}".format("Surface", "Cp"))
    print("".join(["-"]*40))
    print("{0:<20}{1:<20}".format(2, Cp2))
    print("{0:<20}{1:<20}".format(3, Cp3))
    print("{0:<20}{1:<20}".format(4, Cp4))
    print("{0:<20}{1:<20}".format(5, Cp5))