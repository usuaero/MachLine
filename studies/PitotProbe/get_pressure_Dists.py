import os
import time
import numpy as np
import matplotlib.pyplot as plt
import json
import subprocess as sp
import generate_meshes


RERUN_MACHLINE = True


def run_machline_for_loc(x_dist,y_dist,z_dist,index,wake_type,freestream_mach,study_directory):
    
    # default values
    alpha = 5
    #generate geometry
    area = generate_meshes.gen_multi_wing_geom(x_dist,y_dist,z_dist,study_directory,index)
    # Storage locations
    
    mesh_file = study_directory+"/meshes/multi_lifting_surface_"+str(index)+".stl"
    results_file = study_directory+"/results/"+str(index)+".vtk"
    wake_file = study_directory+"/results/"+str(index)+"_wake.vtk"
    report_file = study_directory+"/reports/"+str(index)+".json"

    # create input file
    input_dict = {
        "flow": {
            "freestream_velocity": [np.cos(np.radians(alpha)), 0.0, np.sin(np.radians(alpha))],
            "freestream_mach_number" : freestream_mach
            
        },
        "geometry": {
            "file": mesh_file,
            "spanwise_axis" : "+y",
            "reference" : {
                "area" : area
            },
            "wake_model" : {
                "wake_present" : True,
                "append_wake" : True,
                "trefftz_distance": 60,
                "wake_type": wake_type
            }
        },
    
        "solver": {
        "formulation" : "neumann-mass-flux",
        "ormulation" : "dirichlet-morino",
        "matrix_solver" : "GMRES",
        "write_A_and_b" : False
        },
        "post_processing" : {
            "pressure_rules" : {
                "ncompressible" : True
            }
        },
        "output" : {
            "body_file" : results_file,
            "wake_file" : wake_file,
            "report_file" : report_file
        }
    }
    # Dump
    input_file = "studies/filament_studies/aft_surfaces/input.json"
    write_input_file(input_dict, input_file)

    # Run machline
    report = run_machline(input_file, run=RERUN_MACHLINE)

    # Pull out forces
    C_F = np.zeros(3)
    C_M = np.zeros(3)
    C_F[0] = report["total_forces"]["Cx"]
    C_F[1] = report["total_forces"]["Cy"]
    C_F[2] = report["total_forces"]["Cz"]
    C_M[0] = report["total_moments"]["CMx"]
    C_M[1] = report["total_moments"]["CMy"]
    C_M[2] = report["total_moments"]["CMz"]


    # Get system dimension and average characteristic length
    N_sys = report["solver_results"]["system_dimension"]
    l_avg = report["mesh_info"]["average_characteristic_length"]

    return N_sys, l_avg, C_F, C_M

    



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
    study_directory = "studies/filament_studies/aft_surfaces"
    start = time.time()
    num_cases = 30
    N_sys = list(range(num_cases))
    l_avg =list(range(num_cases))
    C_F = list(range(num_cases))
    C_x = list(range(num_cases))
    C_y = list(range(num_cases))
    C_z = list(range(num_cases))
    C_M = list(range(num_cases))
    C_mx = list(range(num_cases))
    C_my = list(range(num_cases))
    C_mz = list(range(num_cases))

    # set parameters
    x_dist = 1.5
    y_dist = 0
    z_dists = np.linspace(0,30,num_cases)
    wake_types = ["panels", "filaments"]
    freestream_machs = [0.25, 1.7]
    # Run cases

    
        
# Run cases
    try:
        for i in range(len(freestream_machs)):
            for j in range(len(wake_types)):
                for k in range(num_cases):
                    z_dist = z_dists[k]
                    wake_type = wake_types[j]
                    freestream_mach = freestream_machs[i]
                    index = i*len(wake_type)*num_cases + j*num_cases + k
                    N_sys[k], l_avg[k], C_F[k], C_M[k] = run_machline_for_loc(x_dist,y_dist,z_dist,index,wake_type,freestream_mach,study_directory)
                for p in range(num_cases):
                    C_x[p] = C_F[p][0]
                    C_y[p] = C_F[p][1]
                    C_z[p] = C_F[p][2]
                    C_mx[p] = C_M[p][0]
                    C_my[p] = C_M[p][1]
                    C_mz[p] = C_M[p][2]
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