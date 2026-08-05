# Disabled until issue https://gitlab.opengeosys.org/ogs/ogs/-/issues/3486 is
# resolved.
set(ISSUE_3486_RESOLVED FALSE)
if((NOT OGS_USE_MPI) AND ISSUE_3486_RESOLVED)
    OgsTest(PROJECTFILE Parabolic/ThermalTwoPhaseFlowPP/HeatPipe/Twophase_HeatPipe_quad_curve_small.prj RUNTIME 10)
endif()

if((NOT OGS_USE_MPI) AND ISSUE_3486_RESOLVED)
    OgsTest(PROJECTFILE Parabolic/ThermalTwoPhaseFlowPP/HeatPipe/Twophase_HeatPipe_quad_curve_large.prj RUNTIME 170)
endif()

if(NOT OGS_USE_MPI)
    OgsTest(PROJECTFILE Parabolic/ThermalTwoPhaseFlowPP/TCEDiffusion/Twophase_TCE_diffusion_1D_small.prj RUNTIME 2)
endif()
