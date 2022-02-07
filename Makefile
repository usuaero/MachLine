# make for MachLine

# Directories
SRC_DIR = ./src
COM_DIR = ./common
BIN_DIR = ./bin

# Find source files
SRCS = $(wildcard $(COM_DIR)/*.f95 $(SRC_DIR)/*.f95)

# Compiler
COMPILER = gfortran

# Flags
FLAGS = -O2 -fdefault-real-8 -fbounds-check #-ffpe-trap=invalid,zero

# Program name
PROGRAM = machline.exe

default:
	$(COMPILER) $(FLAGS) -o $(PROGRAM) \
	common/helpers.f95 \
	common/linked_list.f95 \
	common/math.f95 \
	common/linalg.f95 \
	common/json.f95 \
	common/json_xtnsn.f95 \
	src/flow.f95 \
	src/vertex.f95 \
	src/edge.f95 \
	src/panel.f95 \
	src/vtk.f95 \
	src/stl.f95 \
	src/wake_mesh.f95 \
	src/surface_mesh.f95 \
	src/panel_solver.f95 \
	src/main.f95

# Debug option
debug:
	@echo "SRCS=$(SRCS)"

# Cleanup
clean:
	rm -f *.mod *.exe $(SRC_DIR)/*.mod $(COM_DIR)/*.mod