if (NOT (OGS_USE_MPI OR OGS_USE_LIS))
    if (OGS_USE_MFRONT)
        OgsTest(PROJECTFILE TH2M/M/MultiMaterialEhlers/square_1e1_2_matIDs.prj RUNTIME 1)
        OgsTest(PROJECTFILE TH2M/M/MultiMaterialEhlers/square_1e1_2_matIDs_restart.prj RUNTIME 1)
    endif()
    OgsTest(PROJECTFILE TH2M/M/M_2d_neumann/M_2d_neumann_fd_jac.xml RUNTIME 0.2)
    OgsTest(PROJECTFILE TH2M/HM/flow_fully_saturated.prj RUNTIME 1)
    OgsTest(PROJECTFILE TH2M/HM/flow_fully_saturated_newton.xml RUNTIME 1)
    OgsTest(PROJECTFILE TH2M/HM/flow_fully_saturated_gas.prj RUNTIME 1)
    OgsTest(PROJECTFILE TH2M/HM/flow_fully_saturated_gas_newton.xml RUNTIME 1)
    OgsTest(PROJECTFILE TH2M/HM/Confined_Compression/HM_confined_compression_gas.prj RUNTIME 9)
    OgsTest(PROJECTFILE TH2M/HM/Confined_Compression/HM_confined_compression_liquid.prj RUNTIME 9)
    OgsTest(PROJECTFILE TH2M/HM/Porosity/deformation_dependent_porosity.prj RUNTIME 18)
    OgsTest(PROJECTFILE TH2M/HM/Porosity/deformation_dependent_porosity_swelling.prj RUNTIME 25)
    OgsTest(PROJECTFILE TH2M/HM/Porosity/deformation_temperature_dependent_porosity_swelling.prj RUNTIME 61)
    OgsTest(PROJECTFILE TH2M/THM/Confined_Compression/THM_confined_compression_gas.prj RUNTIME 10)
    OgsTest(PROJECTFILE TH2M/THM/Confined_Compression/THM_confined_compression_liquid.prj RUNTIME 9)
    OgsTest(PROJECTFILE TH2M/THM/sphere/point_heatsource_fd_jac.xml RUNTIME 4)
    OgsTest(PROJECTFILE TH2M/TH/idealGasLaw/compression_gas.prj RUNTIME 1)
    OgsTest(PROJECTFILE TH2M/TH/temperature_PengRobinson/temperature_preos.prj RUNTIME 1)
    OgsTest(PROJECTFILE TH2M/H/pressure_PengRobinson/pressure_preos.prj RUNTIME 1)
    OgsTest(PROJECTFILE TH2M/HM/compression_PengRobinson/compression_preos.prj RUNTIME 1)
    OgsTest(PROJECTFILE TH2M/H2M/EmbeddedFracturePermeability/IfG.prj RUNTIME 9)
    NotebookTest(NOTEBOOKFILE TH2M/H2M/Liakopoulos/ogs-jupyter-lab-h2m-2d-liakopoulos.py RUNTIME 11)
    if(NOT ENABLE_ASAN)
        OgsTest(PROJECTFILE TH2M/H2M/Liakopoulos/liakopoulos_newton.xml RUNTIME 2)
    endif()
    OgsTest(PROJECTFILE TH2M/H2M/OrthotropicSwelling/square.prj RUNTIME 1)
    OgsTest(PROJECTFILE TH2M/H2/mcWhorter/mcWhorter_h2.prj RUNTIME 12)
    OgsTest(PROJECTFILE TH2M/H2/mcWhorter/mcWhorter_h2_fd_jac.xml RUNTIME 7)
    OgsTest(PROJECTFILE TH2M/H2/mcWhorter/mcWhorter_h2_newton.xml RUNTIME 2)
    OgsTest(PROJECTFILE TH2M/TH2/unicube/unicube.prj RUNTIME 8)
    OgsTest(PROJECTFILE TH2M/TH2/heatpipe/heat_pipe_rough.prj RUNTIME 90)
    OgsTest(PROJECTFILE TH2M/TH2/heatpipe/heat_pipe_strict.prj RUNTIME 11)
    OgsTest(PROJECTFILE TH2M/H2/dissolution_diffusion/continuous_injection.prj RUNTIME 13)
    OgsTest(PROJECTFILE TH2M/H2/dissolution_diffusion/bourgeat.prj RUNTIME 8)
    NotebookTest(NOTEBOOKFILE TH2M/H2/dissolution_diffusion/phase_appearance.py RUNTIME 48)
    NotebookTest(NOTEBOOKFILE TH2M/H2/mcWhorter/mcWhorter.py RUNTIME 16)
    OgsTest(PROJECTFILE TH2M/H/diffusion/diffusion.prj RUNTIME 2)
    NotebookTest(NOTEBOOKFILE TH2M/H/diffusion/diffusion.py RUNTIME 18)
    OgsTest(PROJECTFILE TH2M/TH/Ogata-Banks/ogata-banks.prj RUNTIME 2)
    NotebookTest(NOTEBOOKFILE TH2M/TH/Ogata-Banks/Ogata-Banks.py RUNTIME 33)
    NotebookTest(NOTEBOOKFILE TH2M/TH/idealGasLaw/confined_gas_compression.py RUNTIME 6)
    # submesh residuum output
    OgsTest(PROJECTFILE TH2M/submesh_residuum_assembly/T.xml RUNTIME 1)
    OgsTest(PROJECTFILE TH2M/submesh_residuum_assembly/p_G.xml RUNTIME 1)
    OgsTest(PROJECTFILE TH2M/submesh_residuum_assembly/p_cap.xml RUNTIME 1)
    OgsTest(PROJECTFILE TH2M/submesh_residuum_assembly/u.xml RUNTIME 1)
    OgsTest(PROJECTFILE TH2M/TotalInitialStress/total_initial_stress_H2M.prj RUNTIME 1)
endif()

if(NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE TH2M/H2M/Liakopoulos/liakopoulos_TH2M.prj RUNTIME 4)
    OgsTest(PROJECTFILE TH2M/H2M/Liakopoulos/liakopoulos_fd_jac.xml RUNTIME 4)
endif()

# TH2M 1d heat diffusion w/ Dirichlet-BC
if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE TH2M/T/T_1d_dirichlet/T_1d_dirichlet.prj RUNTIME 3)
endif()
AddTest(
    NAME TH2M_T_1d_dirichlet_newton
    PATH TH2M/T/T_1d_dirichlet
    EXECUTABLE ogs
    EXECUTABLE_ARGS T_1d_dirichlet_newton.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_MPI OR OGS_USE_LIS)
    DIFF_DATA
    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu T_1d_dirichlet_newton_ts_34_t_4000.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 3e-6 1e-8
    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu T_1d_dirichlet_newton_ts_34_t_4000.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8
    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu T_1d_dirichlet_newton_ts_34_t_4000.000000.vtu temperature_interpolated temperature_interpolated 1e-8 1e-8
    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu T_1d_dirichlet_newton_ts_34_t_4000.000000.vtu displacement displacement 1e-8 1e-8
    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu T_1d_dirichlet_newton_ts_34_t_4000.000000.vtu liquid_pressure_interpolated liquid_pressure_interpolated 3e-6 1e-8
    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu T_1d_dirichlet_newton_ts_34_t_4000.000000.vtu velocity_gas velocity_gas 1e-8 1e-8
    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu T_1d_dirichlet_newton_ts_34_t_4000.000000.vtu velocity_liquid velocity_liquid 1e-8 1e-8
    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu T_1d_dirichlet_newton_ts_34_t_4000.000000.vtu sigma sigma 4e-6 1e-8
    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu T_1d_dirichlet_newton_ts_34_t_4000.000000.vtu epsilon epsilon 1e-8 1e-8
    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu T_1d_dirichlet_newton_ts_34_t_4000.000000.vtu liquid_density liquid_density 1e-8 1e-8
    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu T_1d_dirichlet_newton_ts_34_t_4000.000000.vtu gas_density gas_density 1e-8 1e-8
    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu T_1d_dirichlet_newton_ts_34_t_4000.000000.vtu porosity porosity 1e-8 1e-8
    result_TH2M_T_dirichlet_ts_34_t_4000.000000.vtu T_1d_dirichlet_newton_ts_34_t_4000.000000.vtu saturation saturation 1e-8 1e-8
)

# TH2M 2d linear elastic mechanics w/ neumann BC
if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE TH2M/M/M_2d_neumann/M_2d_neumann.prj)
endif()

AddTest(
    NAME TH2M_M_2d_neumann_newton
    PATH TH2M/M/M_2d_neumann
    EXECUTABLE ogs
    EXECUTABLE_ARGS M_2d_neumann_newton.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_MPI OR OGS_USE_LIS)
    DIFF_DATA
    result_TH2M_M_ts_2_t_2.000000.vtu M_2d_neumann_newton_ts_2_t_2.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 1e-8 1e-8
    result_TH2M_M_ts_2_t_2.000000.vtu M_2d_neumann_newton_ts_2_t_2.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8
    result_TH2M_M_ts_2_t_2.000000.vtu M_2d_neumann_newton_ts_2_t_2.000000.vtu temperature_interpolated temperature_interpolated 1e-8 1e-8
    result_TH2M_M_ts_2_t_2.000000.vtu M_2d_neumann_newton_ts_2_t_2.000000.vtu displacement displacement 1e-8 1e-8
    result_TH2M_M_ts_2_t_2.000000.vtu M_2d_neumann_newton_ts_2_t_2.000000.vtu liquid_pressure_interpolated liquid_pressure_interpolated 1e-8 1e-8
    result_TH2M_M_ts_2_t_2.000000.vtu M_2d_neumann_newton_ts_2_t_2.000000.vtu velocity_gas velocity_gas 1e-8 1e-8
    result_TH2M_M_ts_2_t_2.000000.vtu M_2d_neumann_newton_ts_2_t_2.000000.vtu velocity_liquid velocity_liquid 1e-8 1e-8
    result_TH2M_M_ts_2_t_2.000000.vtu M_2d_neumann_newton_ts_2_t_2.000000.vtu sigma sigma 8e-7 1e-8
    result_TH2M_M_ts_2_t_2.000000.vtu M_2d_neumann_newton_ts_2_t_2.000000.vtu epsilon epsilon 1e-8 1e-8
    result_TH2M_M_ts_2_t_2.000000.vtu M_2d_neumann_newton_ts_2_t_2.000000.vtu liquid_density liquid_density 1e-8 1e-8
    result_TH2M_M_ts_2_t_2.000000.vtu M_2d_neumann_newton_ts_2_t_2.000000.vtu gas_density gas_density 1e-8 1e-8
    result_TH2M_M_ts_2_t_2.000000.vtu M_2d_neumann_newton_ts_2_t_2.000000.vtu porosity porosity 1e-8 1e-8
    result_TH2M_M_ts_2_t_2.000000.vtu M_2d_neumann_newton_ts_2_t_2.000000.vtu saturation saturation 1e-8 1e-8
)

# TH2M THM point_heatsource benchmark
if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE TH2M/THM/sphere/point_heatsource.prj RUNTIME 7)
endif()

AddTest(
    NAME TH2M_THM_point_heatsource_newton
    PATH TH2M/THM/sphere
    RUNTIME 4
    EXECUTABLE ogs
    EXECUTABLE_ARGS point_heatsource_newton.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_MPI OR OGS_USE_LIS)
    DIFF_DATA
    result_point_heatsource_ts_5_t_100000.000000.vtu point_heatsource_newton_ts_5_t_100000.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 8e-7 9e-8
    result_point_heatsource_ts_5_t_100000.000000.vtu point_heatsource_newton_ts_5_t_100000.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8
    result_point_heatsource_ts_5_t_100000.000000.vtu point_heatsource_newton_ts_5_t_100000.000000.vtu temperature_interpolated temperature_interpolated 1e-8 1e-8
    result_point_heatsource_ts_5_t_100000.000000.vtu point_heatsource_newton_ts_5_t_100000.000000.vtu displacement displacement 1e-8 1e-8
    result_point_heatsource_ts_5_t_100000.000000.vtu point_heatsource_newton_ts_5_t_100000.000000.vtu liquid_pressure_interpolated liquid_pressure_interpolated 8e-7 9e-8
    result_point_heatsource_ts_5_t_100000.000000.vtu point_heatsource_newton_ts_5_t_100000.000000.vtu velocity_gas velocity_gas 1e-8 1e-8
    result_point_heatsource_ts_5_t_100000.000000.vtu point_heatsource_newton_ts_5_t_100000.000000.vtu velocity_liquid velocity_liquid 1e-8 1e-8
    result_point_heatsource_ts_5_t_100000.000000.vtu point_heatsource_newton_ts_5_t_100000.000000.vtu sigma sigma 4e-7 1e-8
    result_point_heatsource_ts_5_t_100000.000000.vtu point_heatsource_newton_ts_5_t_100000.000000.vtu epsilon epsilon 1e-8 1e-8
    result_point_heatsource_ts_5_t_100000.000000.vtu point_heatsource_newton_ts_5_t_100000.000000.vtu liquid_density liquid_density 1e-8 1e-8
    result_point_heatsource_ts_5_t_100000.000000.vtu point_heatsource_newton_ts_5_t_100000.000000.vtu gas_density gas_density 1e-8 1e-8
    result_point_heatsource_ts_5_t_100000.000000.vtu point_heatsource_newton_ts_5_t_100000.000000.vtu porosity porosity 1e-8 1e-8
    result_point_heatsource_ts_5_t_100000.000000.vtu point_heatsource_newton_ts_5_t_100000.000000.vtu saturation saturation 1e-8 1e-8
)

# TH2M Thermohydromechanics in a slab
if(NOT (OGS_USE_MPI OR OGS_USE_LIS))
    OgsTest(PROJECTFILE TH2M/THM/slab/THM_1d_dirichlet.prj RUNTIME 2)
endif()
AddTest(
    NAME TH2M_THM_THM_1d_dirichlet_newton
    PATH TH2M/THM/slab
    RUNTIME 1
    EXECUTABLE ogs
    EXECUTABLE_ARGS THM_1d_dirichlet_newton.xml
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT (OGS_USE_MPI OR OGS_USE_LIS)
    DIFF_DATA
    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu THM_1d_dirichlet_newton_ts_5_t_100000.000000.vtu gas_pressure_interpolated gas_pressure_interpolated 3e-6 2e-8
    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu THM_1d_dirichlet_newton_ts_5_t_100000.000000.vtu capillary_pressure_interpolated capillary_pressure_interpolated 1e-8 1e-8
    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu THM_1d_dirichlet_newton_ts_5_t_100000.000000.vtu temperature_interpolated temperature_interpolated 1e-8 1e-8
    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu THM_1d_dirichlet_newton_ts_5_t_100000.000000.vtu displacement displacement 1e-8 1e-8
    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu THM_1d_dirichlet_newton_ts_5_t_100000.000000.vtu liquid_pressure_interpolated liquid_pressure_interpolated 3e-6 2e-8
    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu THM_1d_dirichlet_newton_ts_5_t_100000.000000.vtu velocity_gas velocity_gas 1e-8 1e-8
    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu THM_1d_dirichlet_newton_ts_5_t_100000.000000.vtu velocity_liquid velocity_liquid 1e-8 1e-8
    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu THM_1d_dirichlet_newton_ts_5_t_100000.000000.vtu sigma sigma 4e-6 1e-8
    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu THM_1d_dirichlet_newton_ts_5_t_100000.000000.vtu epsilon epsilon 1e-8 1e-8
    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu THM_1d_dirichlet_newton_ts_5_t_100000.000000.vtu liquid_density liquid_density 1e-8 1e-8
    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu THM_1d_dirichlet_newton_ts_5_t_100000.000000.vtu gas_density gas_density 1e-8 1e-8
    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu THM_1d_dirichlet_newton_ts_5_t_100000.000000.vtu porosity porosity 1e-8 1e-8
    result_1d_dirichlet_slab_ts_5_t_100000.000000.vtu THM_1d_dirichlet_newton_ts_5_t_100000.000000.vtu saturation saturation 1e-8 1e-8
)

if(NOT OGS_USE_MPI)
    OgsTest(
        PROJECTFILE
            TH2M/H2M/StrainDependentPermeability/Strain_Dependent_Permeability_Test.prj
        RUNTIME 19
    )
endif()

# PETSc
if(OGS_USE_MPI)
    OgsTest(PROJECTFILE TH2M/THM/sphere/point_heatsource_mpi.xml WRAPPER mpirun
                                                                         -np 3
            RUNTIME 58
    )
endif()

if(OGS_USE_MPI)
    OgsTest(PROJECTFILE TH2M/TH2/heatpipe/PETSc/heat_pipe_strict.prj
            WRAPPER mpirun -np 3 RUNTIME 15
    )
endif()

if(NOT (OGS_USE_PETSC OR OGS_USE_LIS))
    NotebookTest(NOTEBOOKFILE TH2M/TH2/heatpipe/heatpipe.py RUNTIME 194)
    NotebookTest(NOTEBOOKFILE TH2M/H2/mcWhorter_interactive/mcWhorter_interactive.py RUNTIME 18)
    if (OGS_USE_MFRONT)
        NotebookTest(NOTEBOOKFILE TH2M/ExcavationTH2M/excavation_th2m.py RUNTIME 117)
    endif()
endif()
