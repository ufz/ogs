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
if(NOT OGS_USE_MPI)
    OgsTest(
        PROJECTFILE LIE/HydroMechanics/GreatCellWithBBar/HM2b_LIE_F.prj
        RUNTIME 30
        EXECUTABLE_ARGS
            -m ${Data_SOURCE_DIR}/LIE/HydroMechanics/GreatCellWithBBar/mesh_GreatCell_embeddedFracture
    )
endif()

if(NOT (OGS_USE_PETSC OR OGS_USE_LIS))
    NotebookTest(NOTEBOOKFILE LIE/HydroMechanics/GreatCellWithBBar/great_cell_LIE.py RUNTIME 60)
    NotebookTest(NOTEBOOKFILE LIE/Mechanics/GreatCelljupyterNotebook/GreatCellM.py RUNTIME 280)
    NotebookTest(NOTEBOOKFILE LIE/HydroMechanics/GreatCelljupyterNotebook/GreatCellHM.py RUNTIME 851)
endif()
