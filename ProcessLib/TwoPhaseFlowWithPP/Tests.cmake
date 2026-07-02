# Run time takes longer if OGS_USE_LIS is set
if(NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE Parabolic/TwoPhaseFlowPP/Liakopoulos/TwoPhase_Lia_quad_short.prj RUNTIME 2)
    OgsTest(PROJECTFILE Parabolic/TwoPhaseFlowPP/Liakopoulos/TwoPhase_Lia_quad_large.prj RUNTIME 6)
    OgsTest(PROJECTFILE Parabolic/TwoPhaseFlowPP/McWhorter/TwoPhase_mcwt_line.prj RUNTIME 5)
endif()
