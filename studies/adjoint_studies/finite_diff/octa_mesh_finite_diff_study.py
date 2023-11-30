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


def run_machline_for_loc(point_index, component_index,study_directory):
    
    # default values
    

    # Storage locations
    
    mesh_file = study_directory+"/meshes/octa_mesh.stl"
    report_file = study_directory+"/reports/octa_mesh_results.json"

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
            "perturb_point": True,
            "perturbation": {
                "step": 0.00001,
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
                "report_file": report_file
            }
}
       
    # Dump
    input_file = "studies/adjoint_studies/finite_diff/input_files/octa_mesh_finite_diff_input.json"
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
    study_directory = "studies/adjoint_studies/finite_diff"
    start = time.time()
    #num_cases = 30
    #N_sys = list(range(num_cases))
    #l_avg = list(range(num_cases))
    C_F = list(range(num_cases))
    C_x = list(range(num_cases))
    C_y = list(range(num_cases))
    C_z = list(range(num_cases))
    #C_M = list(range(num_cases))
    #C_mx = list(range(num_cases))
    #C_my = list(range(num_cases))
    #C_mz = list(range(num_cases))

    # set parameters
    

    
        
# Run cases
    try:
        for i in range(8): # 8 points in the octa mesh
            i += 1  # change to fortran 1 based indexing
            for j in range(3):
                j += 1  # change to fortran 1 based indexing 
                
                C_F[] = run_machline_for_loc(i,j,study_directory)
            
            
            for p in range(num_cases):
                C_x[p] = C_F[p][0]
                C_y[p] = C_F[p][1]
                C_z[p] = C_F[p][2]
                #C_mx[p] = C_M[p][0]
                #C_my[p] = C_M[p][1]
                #C_mz[p] = C_M[p][2]
            output_file = study_directory
            output_file += "/imp_test_" + str(wake_type) + "_" + str(freestream_mach) + ".txt"
            with open(output_file, "w") as file:
                # Write the lists to the file
                file.write("z_dist,C_x,C_y,C_z,C_mx,C_my,C_mz\n")
                for q in range(num_cases):
                    file.write(f"{z_dists[q]},{C_x[q]},{C_y[q]},{C_z[q]},{C_mx[q]},{C_my[q]},{C_mz[q]}\n")
                file.close()

    except KeyboardInterrupt:
        end = time.time()
        seconds = end-start
        hours =  seconds // 3600
        seconds %= 3600
        minutes = seconds // 60
        seconds %= 60
        error_filename = study_directory+"/Keyboard_exception_output.txt"
        print(f"Code has stopped due to keyboard interrupt at {hours}:{minutes}:{seconds}, data has been written to {error_filename}.")
        for p in range(num_cases):
            try:
                C_x[p] = C_F[p][0]
                C_y[p] = C_F[p][1]
                C_z[p] = C_F[p][2]
                C_mx[p] = C_M[p][0]
                C_my[p] = C_M[p][1]
                C_mz[p] = C_M[p][2]
            except:
                break
        with open(error_filename, "w") as file:
            # Write the lists to the file
            file.write("z_dist,C_x,C_y,C_z,C_mx,C_my,C_mz\n")
            for q in range(len(C_x)):
                file.write(f"{z_dists[q]},{C_x[q]},{C_y[q]},{C_z[q]},{C_mx[q]},{C_my[q]},{C_mz[q]}\n")
            file.close()
    except Exception as e:
        end = time.time()
        seconds = end-start
        hours =  seconds // 3600
        seconds %= 3600
        minutes = seconds // 60
        seconds %= 60
        error_filename = study_directory + "/exception_output.txt"
        print(f"Code has stopped due to exception: {e} at {hours}:{minutes}:{seconds}")
        for p in range(len(C_F)):
                    C_x[p] = C_F[p][0]
                    C_y[p] = C_F[p][1]
                    C_z[p] = C_F[p][2]
                    C_mx[p] = C_M[p][0]
                    C_my[p] = C_M[p][1]
                    C_mz[p] = C_M[p][2]
        with open(error_filename, "w") as file:
            # Write the lists to the file
            file.write("z_dist,C_x,C_y,C_z,C_mx,C_my,C_mz\n")
            for q in range(len(C_x)):
                file.write(f"{z_dists[q]},{C_x[q]},{C_y[q]},{C_z[q]},{C_mx[q]},{C_my[q]},{C_mz[q]}\n")
            file.close()
        print(f"Data has been written to {error_filename}.")
         
    end = time.time()
    seconds = end-start
    hours =  seconds // 3600
    seconds %= 3600
    minutes = seconds // 60
    seconds %= 60
    print(f"Study has complete after {int(hours)}:{int(minutes)}:{int(seconds)}")
    print(f"Data has been written to {output_file}")