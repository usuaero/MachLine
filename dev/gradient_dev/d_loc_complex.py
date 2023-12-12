import subprocess as sp
# this script tests the partial derivatives of the mesh vertices x,y,z coordinates with respect 
# to all design variables X(beta)

if __name__=="__main__":

    # Compile
    sp.run(["gfortran", "-fdefault-real-8", "-fbounds-check", "common/json.f90", "common/json_xtnsn.f90", "common/linalg.f90", "common/helpers.f90", "common/math.f90", "common/linked_list.f90", "common/complexify.f90", "src/base_geom.f90", "dev/gradient_dev/d_loc_complex.f90"])

    # Run
    sp.run(["a.exe"])