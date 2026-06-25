# LIE; HydroMechanics
if (NOT (OGS_USE_LIS OR OGS_USE_MPI))
    OgsTest(PROJECTFILE LIE/HydroMechanics/single_fracture_constK.prj RUNTIME 7)
    OgsTest(PROJECTFILE LIE/HydroMechanics/single_fracture_CZ_kf_const.prj RUNTIME 191)
    OgsTest(PROJECTFILE LIE/HydroMechanics/single_fracture_CZ_kf_cubic.prj RUNTIME 86)
    OgsTest(PROJECTFILE LIE/HydroMechanics/single_fracture_LE_kf_const.prj RUNTIME 67)
    OgsTest(PROJECTFILE LIE/HydroMechanics/single_fracture_LE_kf_cubic.prj RUNTIME 43)
endif()

if(NOT (OGS_USE_LIS OR OGS_USE_MPI))
    OgsTest(PROJECTFILE LIE/HydroMechanics/single_fracture.prj RUNTIME 6)

    OgsTest(PROJECTFILE LIE/HydroMechanics/single_fracture_3D.prj RUNTIME 75)

    OgsTest(PROJECTFILE LIE/HydroMechanics/TaskB.prj RUNTIME 19)

    OgsTest(
        PROJECTFILE LIE/HydroMechanics/single_fracture_3compartments_flow.prj
    )

    OgsTest(
        PROJECTFILE
            LIE/HydroMechanics/single_fracture_3compartments_flow_linear_aperture0.prj
    )
endif()

# Same as the LIE_HM_single_fracture_3compartments_flow_linear_aperture0 but with
# aperture0 defined on the elements (and discontinuous on the nodes).
if(NOT (OGS_USE_LIS OR OGS_USE_MPI))
    OgsTest(
        PROJECTFILE
            LIE/HydroMechanics/single_fracture_3compartments_flow_linear_aperture0_e.prj
    )

    OgsTest(
        PROJECTFILE
            LIE/HydroMechanics/single_fracture_3compartments_flow_CHZ.prj
        RUNTIME 10
    )

    OgsTest(
        PROJECTFILE
            LIE/HydroMechanics/single_fracture_3compartments_flow_CHZ_sigma0.prj
    )

    OgsTest(
        PROJECTFILE
            LIE/HydroMechanics/GreatCellWithTrianglularMesh/great_cell_2d_HM_LIE_embedded_fracture.prj
        RUNTIME 2
    )
endif()

# Note: OGS_USE_LIS takes longer runtime.
AddTest(
    NAME LIE_HydroMechanics_GreatCellWithBBar/HM2b_LIE_F
    PATH LIE/HydroMechanics/GreatCellWithBBar
    EXECUTABLE ogs
    EXECUTABLE_ARGS HM2b_LIE_F.prj -m mesh_GreatCell_embeddedFracture
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 30
    DIFF_DATA
    HM2b_LIE_F_Greywacke_ts_69_t_3500.000000.vtu HM2b_LIE_F_Greywacke_ts_69_t_3500.000000.vtu pressure pressure 2e-8 5.2e-14
    HM2b_LIE_F_Greywacke_ts_69_t_3500.000000.vtu HM2b_LIE_F_Greywacke_ts_69_t_3500.000000.vtu displacement displacement 1e-15 1e-15
    HM2b_LIE_F_Greywacke_ts_69_t_3500.000000.vtu HM2b_LIE_F_Greywacke_ts_69_t_3500.000000.vtu displacement_jump1 displacement_jump1 1e-15 1e-15
    HM2b_LIE_F_Greywacke_ts_69_t_3500.000000.vtu HM2b_LIE_F_Greywacke_ts_69_t_3500.000000.vtu fracture_stress fracture_stress 3.5e-9 0
    HM2b_LIE_F_Greywacke_ts_69_t_3500.000000.vtu HM2b_LIE_F_Greywacke_ts_69_t_3500.000000.vtu fracture_permeability fracture_permeability 1e-15 1e-15
    HM2b_LIE_F_Greywacke_ts_69_t_3500.000000.vtu HM2b_LIE_F_Greywacke_ts_69_t_3500.000000.vtu fracture_aperture fracture_aperture 1e-15 1e-15
    HM2b_LIE_F_Greywacke_ts_69_t_3500.000000.vtu HM2b_LIE_F_Greywacke_ts_69_t_3500.000000.vtu sigma sigma 5.1e-7 0
    HM2b_LIE_F_Greywacke_ts_69_t_3500.000000.vtu HM2b_LIE_F_Greywacke_ts_69_t_3500.000000.vtu epsilon epsilon 1e-15 1e-15
    HM2b_LIE_F_Greywacke_ts_69_t_3500.000000.vtu HM2b_LIE_F_Greywacke_ts_69_t_3500.000000.vtu velocity velocity 1e-15 1e-15
)

if(NOT (OGS_USE_PETSC OR OGS_USE_LIS))
    NotebookTest(NOTEBOOKFILE LIE/HydroMechanics/GreatCellWithBBar/great_cell_LIE.py RUNTIME 60)
    NotebookTest(NOTEBOOKFILE LIE/Mechanics/GreatCelljupyterNotebook/GreatCellM.py RUNTIME 280)
    NotebookTest(NOTEBOOKFILE LIE/HydroMechanics/GreatCelljupyterNotebook/GreatCellHM.py RUNTIME 851)
endif()
