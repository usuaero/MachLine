import subprocess as sp
# this script tests various panel sensitivities (dependent on flow) in the adjoint gradient calculation, compares them to FD
# to all design variables X(beta)

if __name__=="__main__":

     # Compile
    sp.run(["gfortran", "-fdefault-real-8", "-fbounds-check", "common/json.f90", "common/json_xtnsn.f90", "common/linalg.f90", "common/helpers.f90", "common/math.f90", "common/linked_list.f90", "src/adjoint.f90", "src/base_geom.f90", "src/panel.f90", "src/flow.f90", "dev/unit_tests/gradient_tests/panel_flow_test/panel_flow_test.f90"])

    # Run
    sp.run(["a.exe"])