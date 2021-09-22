all:
	gfortran -O2 -fdefault-real-8 -fbounds-check -o mftran.exe \
	src/linked_list.f95 \
	src/math.f95 \
	src/json.f95 \
	src/json_xtnsn.f95 \
	src/flow.f95 \
	src/vertex.f95 \
	src/kutta_edge.f95 \
	src/panel.f95 \
	src/adt.f95 \
	src/vtk.f95 \
	src/surface_mesh.f95 \
	src/panel_solver.f95 \
	src/main.f95