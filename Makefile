all:
	gfortran -O2 -o mftran.exe src/math.f95 src/json.f95 src/json_xtnsn.f95 src/main.f95