import subprocess as sp
# this script tests the partial derivatives of the mesh vertices x,y,z coordinates with respect 
# to all design variables X(beta)

if __name__=="__main__":

    # Compile
    sp.run(["gfortran", "-fdefault-real-8", "-fbounds-check", "common/json.f95", "common/json_xtnsn.f95", "common/linalg.f95", "common/helpers.f95", "common/math.f95", "common/linked_list.f95", "src/flow.f95", "src/base_geom.f95", "src/panel.f95", "dev/AIC_sample.f95"])

    # Run
    sp.run(["./a.out"])