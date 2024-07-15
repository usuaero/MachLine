from tkinter import *
from tkinter import ttk
from tkinter import filedialog
import numpy as np
import json
def check_inputs():
    # declare global variables
    global mach, alpha, beta, trefftz_distance, N_panels, cp_offset, reference_area,errorList
    if not machS.get():
        errorList.append("Mach number must be specified")
    if not alphaS.get():
        errorList.append("Alpha must be specified")
    if not betaS.get():
        errorList.append("Beta must be specified")
    if not trefftz_distanceS.get():
        errorList.append("Trefftz distance must be specified")
    if not N_panelsS.get():
        errorList.append("Number of panels must be specified")
    if not mesh_filename.get():
        errorList.append("Mesh filename must be specified")
    if not formulation.get():
        errorList.append("Formulation must be specified")
    if not cp_offsetS.get():
        errorList.append("Control point offset must be specified")
    if not matrix_solver.get():
        errorList.append("Matrix solver must be specified")
    if not output_directory.get():
        errorList.append("Output directory must be specified")
    if not spanwise_axis.get():
        errorList.append("Spanwise axis must be specified")
    if not reference_areaS.get():
        errorList.append("Reference area must be specified")
    if not wake_type.get():
        errorList.append("Wake type must be specified")
    if not wake_present.get():
        errorList.append("Wake present must be specified")
    if not run_checks.get():
        errorList.append("Run checks must be specified")
    if not pressure_rules.get():
        errorList.append("Pressure rules must be specified")
    # if not study_directory.get():
    #     errorList.append("Study directory must be specified")
    
    # check for valid inputs
    try:
        mach = float(machS.get())
    except ValueError:
        errorList.append("Mach number must be a number")
    try:
        alpha = float(alphaS.get())
    except ValueError:
        errorList.append("Alpha must be a number")
    try:
        beta = float(betaS.get())
    except ValueError:
        errorList.append("Beta must be a number")
    try:
        trefftz_distance = float(trefftz_distanceS.get())
    except ValueError:
        errorList.append("Trefftz distance must be a number")
    try:
        N_panels = int(N_panelsS.get())
    except ValueError:
        errorList.append("Number of panels must be an integer")
    try:
        cp_offset = float(cp_offsetS.get())
    except ValueError:
        errorList.append("Control point offset must be a number")
    try:
        reference_area = float(reference_areaS.get())
    except ValueError:
        errorList.append("Reference area must be a number")
    try:
        wake_present_bool = bool(wake_present.get())
    except ValueError:
        errorList.append("Wake present must be a boolean")
    try:
        run_checks_bool = bool(run_checks.get())
    except ValueError:
        errorList.append("Run checks must be a boolean")
    try:
        write_A_and_b_bool = bool(write_A_and_b.get())
    except ValueError:
        errorList.append("Write A and b must be a boolean")
    try:
        pressure_rules_list = pressure_rules.get()
    except ValueError:
        errorList.append("Pressure rules must be a list")
    # try:
    #     study_directory_str = study_directory.get()
    # except ValueError:
    #     errorList.append("Study directory must be a string")
    try:
        formulation_str = formulation.get()
    except ValueError:
        errorList.append("Formulation must be a string")
    try:
        matrix_solver_str = matrix_solver.get()
    except ValueError:
        errorList.append("Matrix solver must be a string")
    try:
        spanwise_axis_str = spanwise_axis.get()
    except ValueError:
        errorList.append("Spanwise axis must be a string")
    try:
        chordwise_axis_str = chordwise_axis.get()
    except ValueError:
        errorList.append("Chordwise axis must be a string")
    try:
        wake_type_str = wake_type.get()
    except ValueError:
        errorList.append("Wake type must be a string")
    try:
        mesh_filename_str = mesh_filename.get()
    except ValueError:
        errorList.append("Mesh filename must be a string")
    try:
        output_directory_str = output_directory.get()
    except ValueError:
        errorList.append("Output directory must be a string")

    errorText.set(errorList)
    
def get_freestream_velocity(alpha_L, beta_L, spanwise_axis, chordwise_axis):  
    # calculate freestream velocity with spanwise=+y and chordwise=+x
    freestream_velocity = np.zeros(3)
    freestream_velocity[0] = np.cos(np.deg2rad(alpha_L))*np.cos(np.deg2rad(beta_L))
    freestream_velocity[1] = np.sin(np.deg2rad(beta_L))
    freestream_velocity[2] = np.sin(np.deg2rad(alpha_L)*np.cos(np.deg2rad(beta_L)))
    forward = [0.,0.,0.]
    
    # a = about z axis, b = about y axis, c = about x axis
    if chordwise_axis == "+x":
        forward = [1,0,0]
    elif chordwise_axis == "-x":
        forward = [-1,0,0]
    elif chordwise_axis == "+y":
        forward = [1,0,0]
    elif chordwise_axis == "-y":
        forward = [1,0,0]
    elif chordwise_axis == "+z":
        forward = [0,0,1]
    elif chordwise_axis == "-z":
        forward = [0,0,-1]
    else:
        errorList.append("Invalid chordwise axis")
    if spanwise_axis == "+y":
        span = [0,1,0]
    elif spanwise_axis == "-y":
        span = [0,-1,0]
    elif spanwise_axis == "+x":
        span = [1,0,0]
    elif spanwise_axis == "-x":
        span = [-1,0,0]
    elif spanwise_axis == "+z":
        span = [0,0,1]
    elif spanwise_axis == "-z":
        span = [0,0,-1]
    else:
        errorList.append("Invalid spanwise axis")

    span = np.array(span)
    forward = np.array(forward)
    # check to see if span or forward are paralell
    if np.dot(forward, span) == 1:
        errorList.append("Forward and span vectors are parallel")
    down = np.cross(forward,span)
    forward = forward * freestream_velocity[0]
    span = span * freestream_velocity[1]
    down = down * freestream_velocity[2]
    print(forward)
    print(span)
    print(down)

    freestream_velocity = forward + span + down
    return freestream_velocity



def generate_machline_input(*args):
    global errorList
    errorList = []
    # check inputs
    check_inputs()
    # Storage locations
    rootDir = output_directory.get()
    mesh_file = mesh_filename.get()
    results_file = rootDir+"/results.vtk"
    wake_file = rootDir+"/wake_results.vtk"
    report_file = rootDir+"/report.json"
    control_point_file = rootDir+"/control_points.vtk"

    # get velocity
    velocity = get_freestream_velocity(alpha,beta,spanwise_axis.get(),chordwise_axis.get())

    # create input file
    input_dict = {
        "flow": {
            "freestream_velocity": velocity.tolist(),
            "freestream_mach_number": mach
        },
        "geometry": {
            "file": mesh_file,
            "spanwise_axis" : spanwise_axis.get(),
            "reference" : {
                "area" : reference_area
            },
            "wake_model" : {
                "wake_present" : wake_present.get(),
                "append_wake" : True,
                "trefftz_distance" : trefftz_distance,
                "wake_type": wake_type.get(),
                "N_panels" : N_panels
            }
        },
    
        "solver": {
            "formulation" : formulation.get(),
            "control_point_offset": cp_offset,
            "matrix_solver" : matrix_solver.get(),
            "write_A_and_b" : write_A_and_b.get(),
            "run_checks": run_checks.get(),
        },
        "post_processing" : {
            "pressure_rules" : {
                "incompressible" : incompressible,
                "isentropic" : isentropic,
                "second-order": second_order,
                "slender-body":slender_body,
                "linear":linear
            }
        },
        "output" : {
            "body_file" : results_file,
            "wake_file" : wake_file,
            "report_file" : report_file,
            "control_point_file" : control_point_file
        }
    }
    # Dump
    input_file = rootDir+"/input.json"
    write_input_file(input_dict, input_file)


def browse_study_directory():
    
    study_directory = filedialog.askdirectory()
    study_directory_entry.delete(0, END)
    study_directory_entry.insert(0, study_directory)

def browse_input_filename():
    input_filename = filedialog.askopenfilename(filetypes=[("JSON files", "*.json")])
    input_filename_entry.delete(0, END)
    input_filename_entry.insert(0, input_filename)

def browse_mesh_filename():
    mesh_filename = filedialog.askopenfilename(filetypes=[(".vtk", "*.vtk"),(".stl", "*.stl"),(".tri", "*.tri"),("All files", "*.*")])
    mesh_filename_entry.delete(0, END)
    mesh_filename_entry.insert(0, mesh_filename)

def browse_output_directory():
    output_directory = filedialog.askdirectory()
    output_directory_entry.delete(0, END)
    output_directory_entry.insert(0, output_directory)

def write_input_file(input_dict, input_filename):
    """Writes the given input dict to the given file location."""
    global errorList
    try:
        with open(input_filename, 'w') as input_handle:
            json.dump(input_dict, input_handle, indent=4,)
    except:
        errorList.append("Invalid Directory")
        errorText.set(errorList)
        
root = Tk()
root.title("MachLine")
small = PhotoImage(file="32Logo.png")
big = PhotoImage(file="64Logo.png")
smaller = PhotoImage(file="16Logo.png")
root.iconphoto(False, smaller, smaller)
# root.iconbitmap("AeroLabLogo.ico")
# root.iconphoto(False, PhotoImage(file="AeroLabLogo.png"))

mainframe = ttk.Frame(root, padding="3 3 12 12")
mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
mainframe.columnconfigure(0, weight=1)
mainframe.rowconfigure(0, weight=1)

### inputs
# study directory
# study_directory = StringVar()
# study_directory.set("studies/default")
# ttk.Label(mainframe, text="Study Directory").grid(column=1, row=1, sticky=E)
# study_directory_entry = ttk.Entry(mainframe, width=40,textvariable=study_directory)
# study_directory_entry.grid(column=2, row=1, sticky=(W, E))
# ttk.Button(mainframe, text="Browse", command=browse_study_directory).grid(column=3, row=1, sticky=W)
# input file
input_filename = StringVar()
input_filename.set("default.json")
# ttk.Label(mainframe, text="Input filename").grid(column=1, row=1, sticky=E)
# input_filename_entry = ttk.Entry(mainframe, width=40,textvariable=input_filename)
# input_filename_entry.grid(column=2, row=1, sticky=(W, E))
# ttk.Button(mainframe, text="Browse", command=browse_input_filename).grid(column=3, row=1, sticky=W)

# flow section
sectionRow = 0
flowFrame=ttk.Labelframe(mainframe, text="Flow Properties")
flowFrame.grid(column=0, row=sectionRow, sticky=(W,E))
flowFrame.columnconfigure(0, weight=1)
flowFrame.rowconfigure(0, weight=1)
row = 0

# mach
row += 1
ttk.Label(flowFrame, text="Mach").grid(column=0, row=row, sticky=E)
machS = StringVar(value=2.0)
mach_entry = ttk.Entry(flowFrame, width=7, textvariable=machS).grid(column=1, row=row, sticky=(W, E))
# alpha
row += 1
ttk.Label(flowFrame, text="Alpha").grid(column=0, row=row, sticky=E)
alphaS = StringVar(value=0.0)
alpha_entry = ttk.Entry(flowFrame, width=7, textvariable=alphaS)
alpha_entry.grid(column=1, row=row, sticky=(W, E))
# beta
row += 1
ttk.Label(flowFrame, text="Beta").grid(column=0, row=row, sticky=E)
betaS = StringVar(value=0.0)
beta_entry = ttk.Entry(flowFrame, width=7, textvariable=betaS)
beta_entry.grid(column=1, row=row, sticky=(W, E))

## wake model
sectionRow += 1
wakeFrame = ttk.Labelframe(mainframe, text="Wake Model")
wakeFrame.grid(column=0, row=sectionRow, sticky=(W,E))
wakeFrame.columnconfigure(0, weight=1)
wakeFrame.rowconfigure(0, weight=1)
row =0
# wake present
row += 1
ttk.Label(wakeFrame, text="Wake Present").grid(column=0, row=row, sticky=E)
wake_present = BooleanVar(value=True)
wake_present_entry = ttk.Checkbutton(wakeFrame, variable=wake_present)
wake_present_entry.grid(column=1, row=row, sticky=(W, E))
# treffitz distance
row += 1
ttk.Label(wakeFrame, text="Trefftz Distance").grid(column=0, row=row, sticky=E)
trefftz_distanceS = StringVar(value="20")
trefftz_distance_entry = ttk.Entry(wakeFrame, width=7, textvariable=trefftz_distanceS)
trefftz_distance_entry.grid(column=1, row=row, sticky=(W, E))
# wake type
row += 1
ttk.Label(wakeFrame, text="Wake Type").grid(column=0, row=row, sticky=E)
wake_type = StringVar(value="panels")
wake_type_entry = ttk.Combobox(wakeFrame, textvariable=wake_type)
wake_type_entry['values'] = ("panels", "filaments")
wake_type_entry.grid(column=1, row=row, sticky=(W, E))
# N panels
row += 1
n_panels_label = StringVar(value="Number of Wake Elements")
ttk.Label(wakeFrame, textvariable=n_panels_label).grid(column=0, row=row, sticky=E)
N_panelsS = StringVar(value="1")
N_panels_entry = ttk.Entry(wakeFrame, width=7, textvariable=N_panelsS)
N_panels_entry.grid(column=1, row=row, sticky=(W, E))


# geometry section
sectionRow += 1
geometryFrame = ttk.Labelframe(mainframe, text="Geometry")
geometryFrame.grid(column=0, row=sectionRow, sticky=(W,E))
geometryFrame.columnconfigure(0, weight=1)
geometryFrame.rowconfigure(0, weight=1)
# mesh file
row = 0
ttk.Label(geometryFrame, text="Mesh filename").grid(column=0, row=row, sticky=E)
mesh_filename = StringVar(value="default.vtk")
mesh_filename_entry = ttk.Entry(geometryFrame, width=40, textvariable=mesh_filename)
mesh_filename_entry.grid(column=1, row=row, sticky=(W, E))
ttk.Button(geometryFrame, text="Browse", command=browse_mesh_filename).grid(column=3, row=row, sticky=W)

# spanwise axis
row += 1
ttk.Label(geometryFrame, text="Spanwise Axis").grid(column=0, row=row, sticky=E)
spanwise_axis = StringVar(value="+y")
spanwise_axis_entry = ttk.Combobox(geometryFrame, textvariable=spanwise_axis)
spanwise_axis_entry['values'] = ("+x", "-x", "+y", "-y", "+z", "-z")
spanwise_axis_entry.grid(column=1, row=row, sticky=(W, E))

# chordwise axis
row += 1
ttk.Label(geometryFrame, text="Chordwise Axis").grid(column=0, row=row, sticky=E)
chordwise_axis = StringVar(value="+x")
chordwise_axis_entry = ttk.Combobox(geometryFrame, textvariable=chordwise_axis)
chordwise_axis_entry['values'] = ("+x", "-x", "+y", "-y", "+z", "-z")
chordwise_axis_entry.grid(column=1, row=row, sticky=(W, E))

# reference area
row += 1
ttk.Label(geometryFrame, text="Reference Area").grid(column=0, row=row, sticky=E)
reference_areaS = StringVar(value="1.0")
reference_area_entry = ttk.Entry(geometryFrame, width=7, textvariable=reference_areaS)
reference_area_entry.grid(column=1, row=row, sticky=(W, E))



## solver
sectionRow += 1
solverFrame = ttk.Labelframe(mainframe, text="Solver")
solverFrame.grid(column=0, row=sectionRow, sticky=(W,E))
solverFrame.columnconfigure(0, weight=1)
solverFrame.rowconfigure(0, weight=1)


# formulation
row = 0
ttk.Label(solverFrame, text="Formulation").grid(column=0, row=row, sticky=E)
formulation = StringVar(value="dirichlet-morino")
formulation_entry = ttk.Combobox(solverFrame, textvariable=formulation)
formulation_entry['values'] = ("dirichlet-morino", "neumann-mass-flux", "neumann-mass-flux-VCP")
formulation_entry.grid(column=1, row=row, sticky=(W, E))

# control point offset
row += 1
ttk.Label(solverFrame, text="Control Point Offset").grid(column=0, row=row, sticky=E)
cp_offsetS = StringVar(value="1e-7")
cp_offset_entry = ttk.Entry(solverFrame, width=7, textvariable=cp_offsetS)
cp_offset_entry.grid(column=1, row=row, sticky=(W, E))

# matrix solver
row += 1
ttk.Label(solverFrame, text="Matrix Solver").grid(column=0, row=row, sticky=E)
matrix_solver = StringVar(value="HHLS")
matrix_solver_entry = ttk.Combobox(solverFrame, textvariable=matrix_solver)
matrix_solver_entry['values'] = ("HHLS", "GMRES")
matrix_solver_entry.grid(column=1, row=row, sticky=(W, E))

# checks
row += 1
ttk.Label(solverFrame, text="Run Checks").grid(column=0, row=row, sticky=E)
run_checks = BooleanVar(value=True)
run_checks_entry = ttk.Checkbutton(solverFrame, variable=run_checks)
run_checks_entry.grid(column=1, row=row, sticky=(W, E))

# write A and b
row += 1
ttk.Label(solverFrame, text="Write A and b").grid(column=0, row=row, sticky=E)
write_A_and_b = BooleanVar(value=False)
write_A_and_b_entry = ttk.Checkbutton(solverFrame, variable=write_A_and_b)
write_A_and_b_entry.grid(column=1, row=row, sticky=(W, E))

# post processing
sectionRow += 1
postProcessingFrame = ttk.Labelframe(mainframe, text="Post Processing")
postProcessingFrame.grid(column=0, row=sectionRow, sticky=(W,E))
postProcessingFrame.columnconfigure(0, weight=1)
postProcessingFrame.rowconfigure(0, weight=1)

#pressure rules
row = 0
ttk.Label(postProcessingFrame, text="Pressure Rules").grid(column=0, row=row, sticky=(N,E))
pressure_rules = StringVar(value=["incompressible","isentropic","second-order","slender-body","linear"])
pressure_rules_entry = Listbox(postProcessingFrame, listvariable=pressure_rules,height=5,selectmode=EXTENDED)
pressure_rules_entry.grid(column=1, row=row, sticky=(W, E))
if pressure_rules_entry.selection_includes(0):
    incompressible = True
elif pressure_rules_entry.selection_includes(1):
    isentropic = True
elif pressure_rules_entry.selection_includes(2):
    second_order = True
elif pressure_rules_entry.selection_includes(3):
    slender_body = True
elif pressure_rules_entry.selection_includes(4):
    linear
else:
    incompressible = False
    isentropic = True
    second_order = False
    slender_body = False
    linear = False

# output
sectionRow += 1
outputFrame = ttk.Labelframe(mainframe, text="Output")
outputFrame.grid(column=0, row=sectionRow, sticky=(W,E))
outputFrame.columnconfigure(0, weight=1)
outputFrame.rowconfigure(0, weight=1)

# directory
row = 0
ttk.Label(outputFrame, text="Output Directory").grid(column=0, row=row, sticky=E)
output_directory = StringVar()
output_directory.set("studies/default/results")
output_directory_entry = ttk.Entry(outputFrame, width=40, textvariable=output_directory)
output_directory_entry.grid(column=1, row=row, sticky=(W, E))
ttk.Button(outputFrame, text="Browse", command=browse_output_directory).grid(column=3, row=row, sticky=W)














sectionRow += 1
ttk.Button(mainframe, text="Generate Input", command=generate_machline_input).grid(column=0, row=sectionRow, sticky=W)


errorText = StringVar()
ttk.Label(mainframe,textvariable=errorText)



for child in mainframe.winfo_children(): 
    child.grid_configure(padx=5, pady=5)

# root.bind('<Return>', generate_machline_input)

root.mainloop()
