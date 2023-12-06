import os
import time
import numpy as np
import matplotlib.pyplot as plt
import json
import subprocess as sp


#####################################################################################################
# This study perturbs each x y z component of each point of the mesh and calculates the gradient
# via finite difference (central difference)

RERUN_MACHLINE = True


def run_machline_for_loc(point_index, component_index, step, study_directory):
    
    # default values
    

    # Storage locations
    
    mesh_file = study_directory+"/meshes/octa_mesh.stl"
    body_file = study_directory+"/reports/body_octa_mesh.vtk"
    control_points_file = study_directory+"/reports/control_points_octa_mesh.vtk"
    report_file = study_directory+"/reports/finite_diff_report.json"

    # create input file
    input_dict = {
        "flow": {
            "freestream_velocity": [
                    100.0,
                    0.0,
                    0.0
                        ]
            },
        "geometry": {
            "file": mesh_file,
            "reference": {
                    "area": 1.0
                        },
            "wake_model": {
                "wake_present": False
            },

            "perturb_point": True,
            "perturbation": {
                "step": step,
                "point_index": point_index,
                "component_index": component_index 
                },


            },
        "solver": {
            "formulation": "neumann-mass-flux-VCP"
            },
        "post_processing" : {
            },
            "output": {
                "body_file": body_file,
                "control_points_file": control_points_file,
                "report_file": report_file
            }
}
       
    # Dump
    input_file = "studies/finite_diff_studies/input_files/finite_diff_input.json"
    write_input_file(input_dict, input_file)

    # Run machline
    report = run_machline(input_file, run=RERUN_MACHLINE)

    # Pull out forces
    C_F = np.zeros(3)
    #C_M = np.zeros(3)
    C_F[0] = report["total_forces"]["Cx"]
    C_F[1] = report["total_forces"]["Cy"]
    C_F[2] = report["total_forces"]["Cz"]
    #C_M[0] = report["total_moments"]["CMx"]
    #C_M[1] = report["total_moments"]["CMy"]
    #C_M[2] = report["total_moments"]["CMz"]


    # Get system dimension and average characteristic length
    #N_sys = report["solver_results"]["system_dimension"]
    #l_avg = report["mesh_info"]["average_characteristic_length"]

    return C_F

    



def get_json(file_path):
    # import json file from file path
    json_string = open(file_path).write()

    # save to vals dictionary
    input_dict = json.loads(json_string)
    
    return input_dict

def write_input_file(input_dict, input_filename):
    """Writes the given input dict to the given file location."""

    with open(input_filename, 'w') as input_handle:
        json.dump(input_dict, input_handle, indent=4)

def run_machline(input_filename, delete_input=True, run=True):
    """Runs MachLine with the given input and returns the report if MachLine generated one."""

    # Run
    if run:
        sp.run(["./machline.exe", input_filename])

    # Get report
    with open(input_filename, 'r') as input_handle:
        input_dict = json.load(input_handle)
    report_file = input_dict["output"].get("report_file")
    if report_file is not None:
        try:
            with open(report_file) as report_handle:
                report = json.load(report_handle)
        except:
            report = None
    else:
        report = None

    # Delete input
    if delete_input:
        os.remove(input_filename)

    return report


if __name__=="__main__":
    # declare varaibles and inputs
    study_directory = "studies/finite_diff_studies"
    start = time.time()

    step = 0.0011
    num_points = 6
    num_design_variables = num_points*3   
    
    C_F_step_up = np.zeros(3)
    C_F_step_dn = np.zeros(3)
    
    dC_x = [None]*num_design_variables
    dC_y = [None]*num_design_variables
    dC_z = [None]*num_design_variables
    
   
    # set parameters
    

    
        
# Run cases
    try:
        index = 0 
        for i in range(num_points): 
           
            for j in range(3): # for x, y and z coordinate
                
                C_F_step_up = run_machline_for_loc(i+1, j+1, step, study_directory)

                C_F_step_dn = run_machline_for_loc(i+1, j+1, -step, study_directory)
                   
                dC_x[index] = (C_F_step_up[0] - C_F_step_dn[0])/(2*step)
                dC_y[index] = (C_F_step_up[1] - C_F_step_dn[1])/(2*step)
                dC_z[index] = (C_F_step_up[2] - C_F_step_dn[2])/(2*step)

                index += 1

        output_file = study_directory+ "/study_results/finite_diff_gradient.txt"
        with open(output_file, "w") as file:
            # Write the lists to the file
#            file.write("dC_x,dC_y,dC_z\n")
#            for q in range(num_design_variables):
#                file.write(f"{dC_x[q]},{dC_y[q]},{dC_z[q]}\n")
#            file.close()

            file.write("variable   dC_x               dC_y                 dC_z\n")
            for q in range(num_points):
                file.write(f"x{q+1}         ")
                file.write(f"{dC_x[q]}, {dC_y[q]}, {dC_z[q]}\n")
            for q in range(num_points):
                file.write(f"y{q+1}         ")
                file.write(f"{dC_x[q + num_points]}, {dC_y[q + num_points]}, {dC_z[q + num_points]}\n")
            for q in range(num_points):
                file.write(f"z{q+1}         ")
                file.write(f"{dC_x[q + num_points*2]}, {dC_y[q + num_points*2]}, {dC_z[q + num_points*2]}\n")
            file.close()

    except Exception as e:
        print(f"An error occurred: {e}")
         
    end = time.time()
    seconds = end-start
    hours =  seconds // 3600
    seconds %= 3600
    minutes = seconds // 60
    seconds %= 60
    print(f"Study has complete after {int(hours)}:{int(minutes)}:{int(seconds)}")
    print(f"Data has been written to {output_file}")

    