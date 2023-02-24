# make for MachLine

# Deactivate implicit rules
.SUFFIXES:

# Directories
SRC_DIR = src
COM_DIR = common
BIN_DIR = bin

# List common files (ordered based on dependency)
COMMON_FILES = helpers.f95 linked_list.f95 math.f95 linalg.f95 json.f95 json_xtnsn.f95 sort.f95
COMMON_PATHS = $(addprefix $(COM_DIR)/, $(COMMON_FILES))

# List source files (ordered based on dependency)
SRC_FILES = flow.f95 base_geom.f95 panel.f95 mesh.f95 stl.f95 vtk.f95 tri.f95 wake_strip.f95 wake_mesh.f95 surface_mesh.f95 panel_solver.f95
SRC_PATHS = $(addprefix $(SRC_DIR)/, $(SRC_FILES))

# Main
MAIN_PATH = src/main.f95

# Compiler
COMPILER = gfortran

# Flags
FLAGS = -O2 -fdefault-real-8
OMP_FLAG = -fopenmp
DEBUG_FLAGS = -fbounds-check -fbacktrace

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

# Serial compilation (without OpenMP)
serial:
	$(COMPILER) $(FLAGS) -o $(PROGRAM) $(COMMON_PATHS) $(SRC_PATHS) $(MAIN_PATH)

# Cleanup
clean:
	rm -rf *.mod *.exe $(SRC_DIR)/*.mod $(COM_DIR)/*.mod