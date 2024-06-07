import subprocess as sp
# this script tests calc basic F integral sensitivities in the adjoint gradient calculation, compares them to FD
# to all design variables X(beta)

if __name__=="__main__":

     # Compile
    sp.run(["gfortran", "-fdefault-real-8", "-fbounds-check", "common/json.f90", 
            "common/json_xtnsn.f90", "common/linalg.f90", "common/helpers.f90", 
            "common/math.f90", "common/sort.f90","common/linked_list.f90", 
            "src/adjoint.f90", "src/flow.f90", "src/base_geom.f90", "src/panel.f90", 
            "src/mesh.f90", "src/filament_segment.f90", "src/filament.f90", 
            "src/stl.f90", "src/vtk.f90", "src/filament_mesh.f90", "src/tri.f90", 
            "src/wake_strip.f90", "src/wake_mesh.f90", "src/surface_mesh.f90", 
            "src/panel_solver.f90", 
            
            "dev/unit_tests/dirichlet_gradient_tests/dirichlet_test11_inf_adjoint2/dirichlet_test11_inf_adjoint2.f90"])

    # Run
    sp.run(["a.exe"])