# HMPhaseField; Phase-field fracture, Small deformation, Hydromechanics,
# Hydraulic fracture

# Beam fracture for M and PhaseField processes
if(OGS_USE_MPI)
    OgsTest(PROJECTFILE HMPhaseField/Beam/beam.prj WRAPPER mpirun -np 1 RUNTIME 18)
    OgsTest(PROJECTFILE HMPhaseField/Beam/AT2_iso_tensile.xml WRAPPER mpirun -np 1 RUNTIME 120)
    OgsTest(PROJECTFILE HMPhaseField/Beam/AT1_vd_tensile.xml WRAPPER mpirun -np 1 RUNTIME 18)
    OgsTest(PROJECTFILE HMPhaseField/Beam/bar_COHESIVE_linear.xml RUNTIME 100
            WRAPPER mpirun -np 1
    )
    OgsTest(PROJECTFILE HMPhaseField/Beam/bar_COHESIVE_exponential.xml WRAPPER mpirun -np 1 RUNTIME 80)
endif()

# Small deformation for M process
OgsTest(PROJECTFILE HMPhaseField/Mechanics/bar.prj)
if(OGS_USE_MPI)
    OgsTest(PROJECTFILE HMPhaseField/Mechanics/bar_3D.xml)
    OgsTest(PROJECTFILE HMPhaseField/Mechanics/disc_with_hole.xml)
    # Fluid flow in a fracture for H process
    OgsTest(PROJECTFILE HMPhaseField/LiquidFlowInFracture/fractureflow.prj
            RUNTIME 10 WRAPPER mpirun -np 1
    )
    # Hydromechanics (1D Consolidation) for H and M processes
    OgsTest(PROJECTFILE HMPhaseField/Consolidation/consolidation.prj RUNTIME 20)
    # KGD (toughness-dominated hydraulic fracture) for H, M and PhaseField processes
    OgsTest(PROJECTFILE HMPhaseField/KGD/KGD.prj WRAPPER mpirun -np 4 RUNTIME 300)
    OgsTest(PROJECTFILE HMPhaseField/KGD/KGD_cohesive_linear_softening.xml WRAPPER mpirun -np 4 RUNTIME 300)
endif()

if(OGS_USE_PETSC)
    NotebookTest(NOTEBOOKFILE HMPhaseField/GreatCell/GreatCellHM_VPF.py
        RUNTIME 1300
        PROPERTIES PROCESSORS 2
    )
    NotebookTest(NOTEBOOKFILE HMPhaseField/GreatCell/GreatCellHM_VPF_propagating.py
        RUNTIME 3800
        PROPERTIES PROCESSORS 4
    )
endif()
