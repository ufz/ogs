if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(
        PROJECTFILE Parabolic/T/PicardDamping/undamped.prj
        PROPERTIES WILL_FAIL true
    )
    OgsTest(
        PROJECTFILE Parabolic/T/PicardDamping/damped_039.xml
        RUNTIME 1
    )
    OgsTest(PROJECTFILE Parabolic/T/1D_freezing_column_Stefan/Stefan_problem.prj RUNTIME 2)
    OgsTest(PROJECTFILE Parabolic/T/1D_freezing_column_Stefan/Stefan_problem_homogen.prj RUNTIME 1)
    OgsTest(PROJECTFILE Parabolic/T/2D_freezing_disk/circle_disk.prj RUNTIME 2)
    OgsTest(PROJECTFILE Parabolic/T/2D_Robin/square_1e4_robin.prj RUNTIME 1)
    OgsTest(PROJECTFILE Parabolic/T/2D_Robin/square_1e4_robin_newton.xml RUNTIME 1)
    OgsTest(PROJECTFILE Parabolic/T/2D_Ice_melting-forming_manuf_solution/ManSol3_IceWaterMix_Scaled.prj RUNTIME 1)
    OgsTest(PROJECTFILE Parabolic/T/1D_Two-phase_Stefan_problem_for_ice_melting/Two-phase_Stefan_problem.prj RUNTIME 2)
    OgsTest(PROJECTFILE Parabolic/T/2D_Soil_freezing_round_BHE/m16m15projectB.prj RUNTIME 10)
    OgsTest(PROJECTFILE Parabolic/T/TimeDecayBC/time_decay_bc.prj RUNTIME 1)
    OgsTest(PROJECTFILE Parabolic/T/1D_dirichlet/line_60_heat.prj)
    # variation of the test above with assembly optimization
    OgsTest(PROJECTFILE Parabolic/T/1D_dirichlet/linear/line_60_heat.xml)
    # another variation with more optimizations
    OgsTest(PROJECTFILE Parabolic/T/1D_dirichlet/linear_compute_only_on_dt_change/line_60_heat.xml)
    # a variation with a varying time step size
    OgsTest(PROJECTFILE Parabolic/T/1D_dirichlet/varying_dt/line_60_heat.xml)
    # a variation with a varying time step size and <linear> optimization
    OgsTest(PROJECTFILE Parabolic/T/1D_dirichlet/varying_dt_linear/line_60_heat.xml)
    # a variation with a varying time step size and more optimizations
    OgsTest(PROJECTFILE Parabolic/T/1D_dirichlet/varying_dt_linear_compute_only_on_dt_change/line_60_heat.xml)
    # SQUARE 1x1 HEAT CONDUCTION TEST -- AXIALLY SYMMETRIC
    # The results were compared to an analytical solution (method of manufactured
    # solutions). The vtkdiff comparison is against the numerical solution.
    OgsTest(PROJECTFILE Parabolic/T/2D_axially_symmetric/square_1e2_axi.prj)
    # WEDGE 1x1 HEAT CONDUCTION TEST -- same setup as above test but in cartesian coordinates
    OgsTest(PROJECTFILE Parabolic/T/2D_axially_symmetric/wedge_1e2_axi_ang_0.02.prj)
    # The 25 BHE array benchmark
    # test results are compared to 2D simulation result
    OgsTest(PROJECTFILE Parabolic/T/2D_BHE_array/bhe2d.prj RUNTIME 47)
    OgsTest(PROJECTFILE Parabolic/T/t1_1Dsource/t1_1Dsource.prj)
    OgsTest(PROJECTFILE Parabolic/T/t1_1Dsteady/t1_1Dsteady.prj)
    OgsTest(PROJECTFILE Parabolic/T/t2_1D1bt/t2_1D1bt.prj)
    OgsTest(PROJECTFILE Parabolic/T/t2_1D2bt/t2_1D2bt.prj)
    OgsTest(PROJECTFILE Parabolic/T/1D_line_source_term_tests/line_source_term.prj)
    OgsTest(PROJECTFILE Parabolic/T/1D_line_source_term_tests/moving_source_term.prj)
    # tests for line source term implementation
    OgsTest(PROJECTFILE Parabolic/T/2D_source_term_tests/line_source_term_left/source_term_left.prj)
    # tests for line source term implementation with inclined elements
    OgsTest(PROJECTFILE Parabolic/T/2D_source_term_tests/line_source_term_left/source_term_left_r.prj)
    # For the special setup with a 'dirac' line source term at x=0.5 the
    # analytical solution in 2 dimensions is valid:
    # u(x,y) = -ln(sqrt((x-0.5)^2))/(2 * Pi)
    OgsTest(PROJECTFILE Parabolic/T/2D_source_term_tests/line_source_term_x_0.5/line_source_term_x_0.5.prj)
    OgsTest(PROJECTFILE Parabolic/T/2D_source_term_tests/line_source_term_x_0.5_restricted_to_middle/line_source_term_x_0.5.prj)
    # tests for line source term implementation on a cubic domain
    OgsTest(PROJECTFILE Parabolic/T/3D_line_source_term_tests/3D_line_source_term_middle/line_source_term_x_0.5_y_0.5.prj)
    # tests for line source term implementation on a cubic domain
    OgsTest(PROJECTFILE Parabolic/T/3D_line_source_term_tests/3D_line_source_term_middle_restricted/line_source_term_x_0.5_y_0.5_restricted.prj)
    # tests for line source term implementation on a cylindrical domain
    # For the special setup with a line source term at position (xi, eta) the
    # analytical solution in 2 dimensions is valid:
    # u(x,y) = -ln(sqrt((x-xi)^2+(y-eta)^2))/(2 * Pi)
    OgsTest(PROJECTFILE Parabolic/T/3D_line_source_term_tests/3D_line_source_term_in_cylinder/49k_prisms/line_source_term_in_cylinder.prj)
    OgsTest(
        PROJECTFILE Parabolic/T/3D_line_source_term_tests/3D_line_source_term_in_cylinder/286k_prisms/line_source_term_in_cylinder.prj
        RUNTIME 10
    )
    # Both tests above are executed in this notebook (without diff check). Maybe remove regular tests later?
    NotebookTest(NOTEBOOKFILE Parabolic/T/3D_line_source_term_tests/3D_line_source_term_in_cylinder/heatconduction-line_source_term.py RUNTIME 15)
endif()
if(OGS_USE_MPI)
    OgsTest(PROJECTFILE Parabolic/T/1D_neumann/petsc_newtonls.prj)
    # test the source term on a subdomain - parallel version
    OgsTest(
        PROJECTFILE
            Parabolic/T/1D_dirichlet_source-term/2/line_1_line_1e2_source_term.prj
        WRAPPER mpirun -np 2
    )
endif()
if(NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE Parabolic/T/1D_neumann/picard.prj)
    OgsTest(PROJECTFILE Parabolic/T/1D_neumann/picard_masslumping.prj)
endif()

if(NOT (OGS_USE_MPI OR OGS_USE_LIS) AND NOT ${CMAKE_HOST_SYSTEM_PROCESSOR} STREQUAL "arm64")
    OgsTest(PROJECTFILE Parabolic/T/1D_neumann/newton.prj)
    OgsTest(PROJECTFILE Parabolic/T/1D_neumann/newton_masslumping.prj)
endif()

# test the source term on a subdomain
if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE Parabolic/T/1D_dirichlet_source-term/line_1_line_1e2_source_term.prj)
endif()

# test the source term on a subdomain with the PETSc embedded executable file
if(OGS_USE_MPI)
    OgsTest(
        PROJECTFILE Parabolic/T/1D_dirichlet_source-term/line_1_line_1e2_source_term_petsc.xml
        WRAPPER mpirun -np 1
    )
endif()

# failing test - mesh not found
OgsTest(
    PROJECTFILE Parabolic/T/1D_dirichlet_source-term/line_1_line_1e2_source_term_fail_mesh_not_found.xml
    NO_TEST_DEFINITION
    PROPERTIES PASS_REGULAR_EXPRESSION "Could not read mesh from '.*' file[.] No mesh added[.]"
)
if(NOT OGS_USE_PETSC)
    NotebookTest(NOTEBOOKFILE ../../web/content/docs/tutorials/moving_bc/notebook-moving_bc.py RUNTIME 13)
endif()
