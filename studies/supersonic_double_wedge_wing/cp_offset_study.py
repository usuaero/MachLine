import subprocess as sp
import numpy as np
import matplotlib.pyplot as plt
import json

if __name__=="__main__":

    # Sphere
    sphere_input = {
        "flow" : {
            "freestream_velocity" : [0.0, 0.0, 10.0]
        },
        "geometry" : {
            "file" : "sphere.vtk",
            "wake_model" : {
                "wake_present" : False
            }
        },
        "solver" : {
            "formulation" : "source-free",
            "control_point_offset" : 1.0e-5
        },
        "post_processing" : {
        },
        "output" : {
            "report_file" : "report.json"
        }
    }

    # Double-wedge
    wing_input = {
        "flow": {
            "freestream_velocity": [ 694.4475790151479, 0.0, 0.0 ],
            "gamma": 1.4,
            "freestream_mach_number": 2.0
        },
        "geometry": {
            "file": "diamond_5_deg_full_medium.stl",
            "spanwise_axis": "+y",
            "wake_model": {
                "append_wake": False
            },
            "reference": {
                "area": 4.0
            }
        },
        "solver": {
            "formulation": "source-free",
            "control_point_offset": 1.1e-08
        },
        "post_processing": {
            "pressure_rules": {
                "isentropic": True
            }
        },
        "output": {
            "report_file" : "report.json"
        }
    }

    # Set up offsets
    offsets = np.logspace(-11.999, 0.001, 13)
    Cx = np.zeros(13)
    Cy = np.zeros(13)
    Cz = np.zeros(13)

    # Choose which file to use
    input_dict = wing_input

    for i, offset in enumerate(offsets):

        # Set offset
        input_dict["solver"]["control_point_offset"] = offset
        with open("input.json", 'w') as input_file_handle:
            json.dump(input_dict, input_file_handle, indent=4)

        # Run
        result = sp.run(["./machline.exe", "input.json"], capture_output=True, text=True)
        
        # Load results
        if "successfully" in result.stdout:
            with open("report.json", 'r') as results_file:
                results = json.load(results_file)
            Cx[i] = results["total_forces"]["Cx"]
            Cy[i] = results["total_forces"]["Cy"]
            Cz[i] = results["total_forces"]["Cz"]
        else:
            Cx[i] = np.nan
            Cy[i] = np.nan
            Cz[i] = np.nan

    # Plot
    plt.figure()
    #plt.plot(offsets, np.abs(Cx), 'k-', label="$C_x$")
    plt.plot(offsets, np.abs(Cy), 'k--', label="$C_z$")
    plt.plot(offsets, np.abs(Cz), 'k-.', label="$C_y$")
    plt.xlabel('Control Point Offset')
    plt.ylabel('Force Coefficient')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend()
    plt.show()

    # Plot
    plt.figure()
    plt.plot(offsets, np.abs(Cx), 'k-', label="$C_x$")
    #plt.plot(offsets, np.abs(Cy), 'k--', label="$C_z$")
    #plt.plot(offsets, np.abs(Cz), 'k-.', label="$C_y$")
    plt.xlabel('Control Point Offset')
    plt.ylabel('Force Coefficient')
    plt.xscale('log')
    #plt.yscale('log')
    plt.legend()
    plt.show()