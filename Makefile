all:
	gfortran -O2 -fdefault-real-8 -o mftran.exe src/math.f95 src/json.f95 src/json_xtnsn.f95 src/geometry.f95 src/vtk.f95 src/mesh.f95 src/main.f95