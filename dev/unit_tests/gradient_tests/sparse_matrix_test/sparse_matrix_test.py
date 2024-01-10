import subprocess as sp
# this script tests sparse matrix operations in adjoint.f90

if __name__=="__main__":

    # Compile
    sp.run(["gfortran", "-fdefault-real-8", "-fbounds-check", "common/json.f90", "common/json_xtnsn.f90", "common/linalg.f90", "common/helpers.f90", "common/math.f90","src/adjoint.f90", "dev/unit_tests/gradient_tests/sparse_matrix_test/sparse_matrix_test.f90"])

    # Run
    sp.run(["a.exe"])