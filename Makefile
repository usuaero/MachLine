# make for MachLine

# Deactivate implicit rules
.SUFFIXES:

# Directories
SRC_DIR = src
COM_DIR = common
BIN_DIR = bin

# List common files (ordered based on dependency)
COMMON_FILES = helpers.f90 linked_list.f90 math.f90 linalg.f90 json.f90 json_xtnsn.f90 sort.f90
COMMON_PATHS = $(addprefix $(COM_DIR)/, $(COMMON_FILES))

# List source files (ordered based on dependency)
SRC_FILES = flow.f90 base_geom.f90 panel.f90 mesh.f90 filament_segment.f90 stl.f90 vtk.f90 tri.f90 wake_strip.f90 filament.f90 wake_mesh.f90 filament_mesh.f90 surface_mesh.f90 panel_solver.f90
SRC_PATHS = $(addprefix $(SRC_DIR)/, $(SRC_FILES))

# Main
MAIN_PATH = src/main.f90

# Compiler
COMPILER = gfortran

# Flags
FLAGS = -O2 -fdefault-real-8
OMP_FLAG = -fopenmp
DEBUG_FLAGS = -fbounds-check -fbacktrace -g

# Program name
PROGRAM = machline.exe

# Default make
default:
	$(COMPILER) $(FLAGS) $(OMP_FLAG) -o $(PROGRAM) $(COMMON_PATHS) $(SRC_PATHS) $(MAIN_PATH)

# Debug option
debug:
	$(COMPILER) $(FLAGS) $(OMP_FLAG) $(DEBUG_FLAGS) -o $(PROGRAM) $(COMMON_PATHS) $(SRC_PATHS) $(MAIN_PATH)

# Debug with all warnings
wall:
	$(COMPILER) $(FLAGS) $(OMP_FLAG) $(DEBUG_FLAGS) -Wall -o $(PROGRAM) $(COMMON_PATHS) $(SRC_PATHS) $(MAIN_PATH)

# Debug with all warnings
debug-serial:
	$(COMPILER) $(FLAGS) $(DEBUG_FLAGS) -o $(PROGRAM) $(COMMON_PATHS) $(SRC_PATHS) $(MAIN_PATH)

# Serial compilation (without OpenMP)
serial:
	$(COMPILER) $(FLAGS) -o $(PROGRAM) $(COMMON_PATHS) $(SRC_PATHS) $(MAIN_PATH)

# Cleanup
clean:
	rm -rf *.mod *.exe $(SRC_DIR)/*.mod $(COM_DIR)/*.mod