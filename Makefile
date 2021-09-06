all:
	gfortran -O2 -o mftran.exe src/main.f95 src/math.f95 src/json.f95 src/myjson.f95