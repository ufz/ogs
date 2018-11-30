# Mechanics; Small deformations Non-Local, Ehlers-damage Bar
AddTest(
    NAME Mechanics_PlasticDamageNonlocalModel_NNLED_Ehlers_Damage_BarCoarse
    PATH Mechanics/EhlersDamageNonLocal/bar/coarse
    EXECUTABLE ogs
    EXECUTABLE_ARGS bar.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu displacement displacement 2e-11 1e-13
    out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu damage damage 2e-11 0
    out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu sigma sigma 5e-8 0
    out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu epsilon epsilon 2e-11 0

    out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu displacement displacement 2e-11 1e-13
    out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu damage damage 2e-11 0
    out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu sigma sigma 1.5e-7 0
    out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu epsilon epsilon 2e-11 0

    out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu displacement displacement 2e-11 1e-13
    out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu damage damage 2e-11 0
    out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu sigma sigma 3e-7 0
    out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu epsilon epsilon 2e-11 0
)

# same as bar/coarse but with tangent type being Elastic.
AddTest(
    NAME
    Mechanics_PlasticDamageNonlocalModel_NNLED_Ehlers_Damage_BarCoarseElasticTangent
    PATH Mechanics/EhlersDamageNonLocal/bar/coarse_elastic_tangent
    EXECUTABLE ogs
    EXECUTABLE_ARGS bar.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu displacement displacement 2e-11 1e-13
    out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu damage damage 3e-10 0
    out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu sigma sigma 5e-3 0
    out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu epsilon epsilon 2e-11 0

    out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu displacement displacement 2e-11 1e-13
    out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu damage damage 2e-10 0
    out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu sigma sigma 5e-3 0
    out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu epsilon epsilon 2e-11 0

    out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu displacement displacement 2e-11 1e-13
    out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu damage damage 2e-10 0
    out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu sigma sigma 2e-3 0
    out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu epsilon epsilon 2e-11 0
)

AddTest(
    NAME Mechanics_PlasticDamageNonlocalModel_NNLED_Ehlers_Damage_BarMedium
    PATH Mechanics/EhlersDamageNonLocal/bar/medium
    EXECUTABLE ogs
    EXECUTABLE_ARGS bar.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    out_bar_medium_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_medium_ED_pcs_0_ts_30_t_0.300000.vtu displacement displacement 2e-11 1e-13
    out_bar_medium_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_medium_ED_pcs_0_ts_30_t_0.300000.vtu damage damage 2e-11 0
    out_bar_medium_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_medium_ED_pcs_0_ts_30_t_0.300000.vtu sigma sigma 1e-7 0
    out_bar_medium_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_medium_ED_pcs_0_ts_30_t_0.300000.vtu epsilon epsilon 2e-11 0

    out_bar_medium_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_medium_ED_pcs_0_ts_60_t_0.600000.vtu displacement displacement 2e-11 1e-13
    out_bar_medium_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_medium_ED_pcs_0_ts_60_t_0.600000.vtu damage damage 2e-11 0
    out_bar_medium_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_medium_ED_pcs_0_ts_60_t_0.600000.vtu sigma sigma 2e-7 0
    out_bar_medium_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_medium_ED_pcs_0_ts_60_t_0.600000.vtu epsilon epsilon 2e-11 0

    out_bar_medium_ED_pcs_0_ts_90_t_0.900000.vtu out_bar_medium_ED_pcs_0_ts_90_t_0.900000.vtu displacement displacement 2e-11 1e-13
    out_bar_medium_ED_pcs_0_ts_90_t_0.900000.vtu out_bar_medium_ED_pcs_0_ts_90_t_0.900000.vtu damage damage 2e-11 0
    out_bar_medium_ED_pcs_0_ts_90_t_0.900000.vtu out_bar_medium_ED_pcs_0_ts_90_t_0.900000.vtu sigma sigma 3e-7 0
    out_bar_medium_ED_pcs_0_ts_90_t_0.900000.vtu out_bar_medium_ED_pcs_0_ts_90_t_0.900000.vtu epsilon epsilon 2e-11 0
)

AddTest(
    NAME LARGE_Mechanics_PlasticDamageNonlocalModel_NNLED_Ehlers_Damage_BarFine
    PATH Mechanics/EhlersDamageNonLocal/bar/fine
    EXECUTABLE ogs
    EXECUTABLE_ARGS bar.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    out_bar_fine_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_fine_ED_pcs_0_ts_30_t_0.300000.vtu displacement displacement 2e-11 1e-13
    out_bar_fine_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_fine_ED_pcs_0_ts_30_t_0.300000.vtu damage damage 2e-11 0
    out_bar_fine_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_fine_ED_pcs_0_ts_30_t_0.300000.vtu sigma sigma 3e-7 0
    out_bar_fine_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_fine_ED_pcs_0_ts_30_t_0.300000.vtu epsilon epsilon 2e-11 0

    out_bar_fine_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_fine_ED_pcs_0_ts_60_t_0.600000.vtu displacement displacement 2e-11 1e-13
    out_bar_fine_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_fine_ED_pcs_0_ts_60_t_0.600000.vtu damage damage 2e-11 0
    out_bar_fine_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_fine_ED_pcs_0_ts_60_t_0.600000.vtu sigma sigma 5e-7 0
    out_bar_fine_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_fine_ED_pcs_0_ts_60_t_0.600000.vtu epsilon epsilon 2e-11 0

    out_bar_fine_ED_pcs_0_ts_90_t_0.900000.vtu out_bar_fine_ED_pcs_0_ts_90_t_0.900000.vtu displacement displacement 2e-11 1e-13
    out_bar_fine_ED_pcs_0_ts_90_t_0.900000.vtu out_bar_fine_ED_pcs_0_ts_90_t_0.900000.vtu damage damage 2e-11 0
    out_bar_fine_ED_pcs_0_ts_90_t_0.900000.vtu out_bar_fine_ED_pcs_0_ts_90_t_0.900000.vtu sigma sigma 6e-7 0
    out_bar_fine_ED_pcs_0_ts_90_t_0.900000.vtu out_bar_fine_ED_pcs_0_ts_90_t_0.900000.vtu epsilon epsilon 2e-11 0
)

AddTest(
    NAME LARGE_Mechanics_PlasticDamageNonlocalModel_NNLED_Ehlers_Damage_BarVeryFine
    PATH Mechanics/EhlersDamageNonLocal/bar/veryfine
    EXECUTABLE ogs
    EXECUTABLE_ARGS bar.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    out_bar_veryfine_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_veryfine_ED_pcs_0_ts_30_t_0.300000.vtu displacement displacement 2e-11 1e-13
    out_bar_veryfine_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_veryfine_ED_pcs_0_ts_30_t_0.300000.vtu damage damage 2e-11 0
    out_bar_veryfine_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_veryfine_ED_pcs_0_ts_30_t_0.300000.vtu sigma sigma 5e-7 0
    out_bar_veryfine_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_veryfine_ED_pcs_0_ts_30_t_0.300000.vtu epsilon epsilon 2e-11 0

    out_bar_veryfine_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_veryfine_ED_pcs_0_ts_60_t_0.600000.vtu displacement displacement 2e-11 1e-13
    out_bar_veryfine_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_veryfine_ED_pcs_0_ts_60_t_0.600000.vtu damage damage 2e-11 0
    out_bar_veryfine_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_veryfine_ED_pcs_0_ts_60_t_0.600000.vtu sigma sigma 1e-6 0
    out_bar_veryfine_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_veryfine_ED_pcs_0_ts_60_t_0.600000.vtu epsilon epsilon 2e-11 0

    out_bar_veryfine_ED_pcs_0_ts_70_t_0.700000.vtu out_bar_veryfine_ED_pcs_0_ts_70_t_0.700000.vtu displacement displacement 2e-11 1e-13
    out_bar_veryfine_ED_pcs_0_ts_70_t_0.700000.vtu out_bar_veryfine_ED_pcs_0_ts_70_t_0.700000.vtu damage damage 2e-11 0
    out_bar_veryfine_ED_pcs_0_ts_70_t_0.700000.vtu out_bar_veryfine_ED_pcs_0_ts_70_t_0.700000.vtu sigma sigma 1e-6 0
    out_bar_veryfine_ED_pcs_0_ts_70_t_0.700000.vtu out_bar_veryfine_ED_pcs_0_ts_70_t_0.700000.vtu epsilon epsilon 2e-11 0
)

AddTest(
    NAME Mechanics_PlasticDamageNonlocalModel_NNLED_Ehlers_Damage_BarCoarse_p1
    PATH Mechanics/EhlersDamageNonLocal/bar/p1
    EXECUTABLE ogs
    EXECUTABLE_ARGS bar.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu displacement displacement 2e-11 1e-13
    out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu damage damage 2e-11 0
    out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu sigma sigma 5e-8 0
    out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu epsilon epsilon 2e-11 0

    out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu displacement displacement 2e-11 1e-13
    out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu damage damage 2e-11 0
    out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu sigma sigma 2e-7 0
    out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu epsilon epsilon 2e-11 0

    out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu displacement displacement 2e-11 1e-13
    out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu damage damage 2e-11 0
    out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu sigma sigma 2e-7 0
    out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu epsilon epsilon 2e-11 0
)

AddTest(
    NAME Mechanics_PlasticDamageNonlocalModel_NNLED_Ehlers_Damage_BarCoarse_p2
    PATH Mechanics/EhlersDamageNonLocal/bar/p2
    EXECUTABLE ogs
    EXECUTABLE_ARGS bar.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu displacement displacement 2e-11 1e-13
    out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu damage damage 2e-11 0
    out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu sigma sigma 7e-8 0
    out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu epsilon epsilon 2e-11 0

    out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu displacement displacement 2e-11 1e-13
    out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu damage damage 2e-11 0
    out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu sigma sigma 2e-7 0
    out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu epsilon epsilon 2e-11 0

    out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu displacement displacement 2e-11 1e-13
    out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu damage damage 2e-11 0
    out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu sigma sigma 2e-7 0
    out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu epsilon epsilon 2e-11 0
)

AddTest(
    NAME Mechanics_PlasticDamageNonlocalModel_NNLED_Ehlers_Damage_BarCoarse_p3
    PATH Mechanics/EhlersDamageNonLocal/bar/p3
    EXECUTABLE ogs
    EXECUTABLE_ARGS bar.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu displacement displacement 2e-11 1e-13
    out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu damage damage 2e-11 0
    out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu sigma sigma 7e-8 0
    out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu out_bar_coarse_ED_pcs_0_ts_30_t_0.300000.vtu epsilon epsilon 2e-11 0

    out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu displacement displacement 2e-11 1e-13
    out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu damage damage 2e-11 0
    out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu sigma sigma 1e-7 0
    out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_bar_coarse_ED_pcs_0_ts_60_t_0.600000.vtu epsilon epsilon 2e-11 0

    out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu displacement displacement 2e-11 1e-13
    out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu damage damage 2e-11 0
    out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu sigma sigma 2e-7 0
    out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu out_bar_coarse_ED_pcs_0_ts_90_t_0.900000.vtu epsilon epsilon 2e-11 0
)

AddTest(
    NAME Mechanics_PlasticDamageNonlocalModel_NNLED_Ehlers_Damage_BarCoarse_t1
    PATH Mechanics/EhlersDamageNonLocal/bar/t1
    EXECUTABLE ogs
    EXECUTABLE_ARGS bar.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    out_bar_coarse_ED_pcs_0_ts_300_t_0.300000.vtu out_bar_coarse_ED_pcs_0_ts_300_t_0.300000.vtu displacement displacement 2e-11 1e-13
    out_bar_coarse_ED_pcs_0_ts_300_t_0.300000.vtu out_bar_coarse_ED_pcs_0_ts_300_t_0.300000.vtu damage damage 2e-11 0
    out_bar_coarse_ED_pcs_0_ts_300_t_0.300000.vtu out_bar_coarse_ED_pcs_0_ts_300_t_0.300000.vtu sigma sigma 5e-8 0
    out_bar_coarse_ED_pcs_0_ts_300_t_0.300000.vtu out_bar_coarse_ED_pcs_0_ts_300_t_0.300000.vtu epsilon epsilon 2e-11 0

    out_bar_coarse_ED_pcs_0_ts_600_t_0.600000.vtu out_bar_coarse_ED_pcs_0_ts_600_t_0.600000.vtu displacement displacement 2e-11 1e-13
    out_bar_coarse_ED_pcs_0_ts_600_t_0.600000.vtu out_bar_coarse_ED_pcs_0_ts_600_t_0.600000.vtu damage damage 2e-11 0
    out_bar_coarse_ED_pcs_0_ts_600_t_0.600000.vtu out_bar_coarse_ED_pcs_0_ts_600_t_0.600000.vtu sigma sigma 2e-7 0
    out_bar_coarse_ED_pcs_0_ts_600_t_0.600000.vtu out_bar_coarse_ED_pcs_0_ts_600_t_0.600000.vtu epsilon epsilon 2e-11 0

    out_bar_coarse_ED_pcs_0_ts_900_t_0.900000.vtu out_bar_coarse_ED_pcs_0_ts_900_t_0.900000.vtu displacement displacement 2e-11 1e-13
    out_bar_coarse_ED_pcs_0_ts_900_t_0.900000.vtu out_bar_coarse_ED_pcs_0_ts_900_t_0.900000.vtu damage damage 2e-11 0
    out_bar_coarse_ED_pcs_0_ts_900_t_0.900000.vtu out_bar_coarse_ED_pcs_0_ts_900_t_0.900000.vtu sigma sigma 2e-7 0
    out_bar_coarse_ED_pcs_0_ts_900_t_0.900000.vtu out_bar_coarse_ED_pcs_0_ts_900_t_0.900000.vtu epsilon epsilon 2e-11 0
)

AddTest(
    NAME Mechanics_PlasticDamageNonlocalModel_NNLED_Ehlers_Damage_BarCoarse_t2
    PATH Mechanics/EhlersDamageNonLocal/bar/t2
    EXECUTABLE ogs
    EXECUTABLE_ARGS bar.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    out_bar_coarse_ED_pcs_0_ts_600_t_0.300000.vtu out_bar_coarse_ED_pcs_0_ts_600_t_0.300000.vtu displacement displacement 2e-11 1e-13
    out_bar_coarse_ED_pcs_0_ts_600_t_0.300000.vtu out_bar_coarse_ED_pcs_0_ts_600_t_0.300000.vtu damage damage 2e-11 0
    out_bar_coarse_ED_pcs_0_ts_600_t_0.300000.vtu out_bar_coarse_ED_pcs_0_ts_600_t_0.300000.vtu sigma sigma 5e-8 0
    out_bar_coarse_ED_pcs_0_ts_600_t_0.300000.vtu out_bar_coarse_ED_pcs_0_ts_600_t_0.300000.vtu epsilon epsilon 2e-11 0

    out_bar_coarse_ED_pcs_0_ts_1200_t_0.600000.vtu out_bar_coarse_ED_pcs_0_ts_1200_t_0.600000.vtu displacement displacement 2e-11 1e-13
    out_bar_coarse_ED_pcs_0_ts_1200_t_0.600000.vtu out_bar_coarse_ED_pcs_0_ts_1200_t_0.600000.vtu damage damage 2e-11 0
    out_bar_coarse_ED_pcs_0_ts_1200_t_0.600000.vtu out_bar_coarse_ED_pcs_0_ts_1200_t_0.600000.vtu sigma sigma 2e-7 0
    out_bar_coarse_ED_pcs_0_ts_1200_t_0.600000.vtu out_bar_coarse_ED_pcs_0_ts_1200_t_0.600000.vtu epsilon epsilon 2e-11 0

    out_bar_coarse_ED_pcs_0_ts_1800_t_0.900000.vtu out_bar_coarse_ED_pcs_0_ts_1800_t_0.900000.vtu displacement displacement 2e-11 1e-13
    out_bar_coarse_ED_pcs_0_ts_1800_t_0.900000.vtu out_bar_coarse_ED_pcs_0_ts_1800_t_0.900000.vtu damage damage 2e-11 0
    out_bar_coarse_ED_pcs_0_ts_1800_t_0.900000.vtu out_bar_coarse_ED_pcs_0_ts_1800_t_0.900000.vtu sigma sigma 2e-7 0
    out_bar_coarse_ED_pcs_0_ts_1800_t_0.900000.vtu out_bar_coarse_ED_pcs_0_ts_1800_t_0.900000.vtu epsilon epsilon 2e-11 0
)

AddTest(
    NAME Mechanics_PlasticDamageNonlocalModel_NNLED_Ehlers_Damage_BeamCoarse
    PATH Mechanics/EhlersDamageNonLocal/beam/coarse
    EXECUTABLE ogs
    EXECUTABLE_ARGS beam.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu displacement displacement 2e-11 1e-13
    out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu damage damage 2e-11 0
    out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu sigma sigma 2e-7 0
    out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu epsilon epsilon 2e-11 0

    out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu displacement displacement 2e-11 1e-13
    out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu damage damage 2e-11 0
    out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu sigma sigma 3e-7 0
    out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu epsilon epsilon 2e-11 0

    out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu displacement displacement 2e-11 1e-13
    out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu damage damage 2e-11 0
    out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu sigma sigma 4e-7 0
    out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu epsilon epsilon 2e-11 0

    out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu displacement displacement 2e-11 1e-13
    out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu damage damage 2e-11 0
    out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu sigma sigma 4e-7 0
    out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu epsilon epsilon 2e-11 0
)

AddTest(
    NAME Mechanics_PlasticDamageNonlocalModel_NNLED_Ehlers_Damage_BeamMedium
    PATH Mechanics/EhlersDamageNonLocal/beam/medium
    EXECUTABLE ogs
    EXECUTABLE_ARGS beam.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    out_beam_medium_ED_pcs_0_ts_40_t_0.400000.vtu out_beam_medium_ED_pcs_0_ts_40_t_0.400000.vtu displacement displacement 2e-11 1e-13
    out_beam_medium_ED_pcs_0_ts_40_t_0.400000.vtu out_beam_medium_ED_pcs_0_ts_40_t_0.400000.vtu damage damage 2e-11 0
    out_beam_medium_ED_pcs_0_ts_40_t_0.400000.vtu out_beam_medium_ED_pcs_0_ts_40_t_0.400000.vtu sigma sigma 4e-7 0
    out_beam_medium_ED_pcs_0_ts_40_t_0.400000.vtu out_beam_medium_ED_pcs_0_ts_40_t_0.400000.vtu epsilon epsilon 2e-11 0

    out_beam_medium_ED_pcs_0_ts_60_t_0.600000.vtu out_beam_medium_ED_pcs_0_ts_60_t_0.600000.vtu displacement displacement 2e-11 1e-13
    out_beam_medium_ED_pcs_0_ts_60_t_0.600000.vtu out_beam_medium_ED_pcs_0_ts_60_t_0.600000.vtu damage damage 2e-11 0
    out_beam_medium_ED_pcs_0_ts_60_t_0.600000.vtu out_beam_medium_ED_pcs_0_ts_60_t_0.600000.vtu sigma sigma 7e-7 0
    out_beam_medium_ED_pcs_0_ts_60_t_0.600000.vtu out_beam_medium_ED_pcs_0_ts_60_t_0.600000.vtu epsilon epsilon 2e-11 0

    out_beam_medium_ED_pcs_0_ts_80_t_0.800000.vtu out_beam_medium_ED_pcs_0_ts_80_t_0.800000.vtu displacement displacement 2e-11 1e-13
    out_beam_medium_ED_pcs_0_ts_80_t_0.800000.vtu out_beam_medium_ED_pcs_0_ts_80_t_0.800000.vtu damage damage 2e-11 0
    out_beam_medium_ED_pcs_0_ts_80_t_0.800000.vtu out_beam_medium_ED_pcs_0_ts_80_t_0.800000.vtu sigma sigma 6e-7 0
    out_beam_medium_ED_pcs_0_ts_80_t_0.800000.vtu out_beam_medium_ED_pcs_0_ts_80_t_0.800000.vtu epsilon epsilon 2e-11 0

    out_beam_medium_ED_pcs_0_ts_100_t_1.000000.vtu out_beam_medium_ED_pcs_0_ts_100_t_1.000000.vtu displacement displacement 2e-11 1e-13
    out_beam_medium_ED_pcs_0_ts_100_t_1.000000.vtu out_beam_medium_ED_pcs_0_ts_100_t_1.000000.vtu damage damage 2e-11 0
    out_beam_medium_ED_pcs_0_ts_100_t_1.000000.vtu out_beam_medium_ED_pcs_0_ts_100_t_1.000000.vtu sigma sigma 9e-7 0
    out_beam_medium_ED_pcs_0_ts_100_t_1.000000.vtu out_beam_medium_ED_pcs_0_ts_100_t_1.000000.vtu epsilon epsilon 2e-11 0
)

AddTest(
    NAME LARGE_Mechanics_PlasticDamageNonlocalModel_NNLED_Ehlers_Damage_BeamFine
    PATH Mechanics/EhlersDamageNonLocal/beam/fine
    EXECUTABLE ogs
    EXECUTABLE_ARGS beam.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    out_beam_fine_ED_pcs_0_ts_40_t_0.400000.vtu out_beam_fine_ED_pcs_0_ts_40_t_0.400000.vtu displacement displacement 2e-11 1e-13
    out_beam_fine_ED_pcs_0_ts_40_t_0.400000.vtu out_beam_fine_ED_pcs_0_ts_40_t_0.400000.vtu damage damage 2e-11 0
    out_beam_fine_ED_pcs_0_ts_40_t_0.400000.vtu out_beam_fine_ED_pcs_0_ts_40_t_0.400000.vtu sigma sigma 7e-7 0
    out_beam_fine_ED_pcs_0_ts_40_t_0.400000.vtu out_beam_fine_ED_pcs_0_ts_40_t_0.400000.vtu epsilon epsilon 2e-11 0

    out_beam_fine_ED_pcs_0_ts_60_t_0.600000.vtu out_beam_fine_ED_pcs_0_ts_60_t_0.600000.vtu displacement displacement 2e-11 1e-13
    out_beam_fine_ED_pcs_0_ts_60_t_0.600000.vtu out_beam_fine_ED_pcs_0_ts_60_t_0.600000.vtu damage damage 2e-11 0
    out_beam_fine_ED_pcs_0_ts_60_t_0.600000.vtu out_beam_fine_ED_pcs_0_ts_60_t_0.600000.vtu sigma sigma 2e-6 0
    out_beam_fine_ED_pcs_0_ts_60_t_0.600000.vtu out_beam_fine_ED_pcs_0_ts_60_t_0.600000.vtu epsilon epsilon 2e-11 0

    out_beam_fine_ED_pcs_0_ts_80_t_0.800000.vtu out_beam_fine_ED_pcs_0_ts_80_t_0.800000.vtu displacement displacement 2e-11 1e-13
    out_beam_fine_ED_pcs_0_ts_80_t_0.800000.vtu out_beam_fine_ED_pcs_0_ts_80_t_0.800000.vtu damage damage 2e-11 0
    out_beam_fine_ED_pcs_0_ts_80_t_0.800000.vtu out_beam_fine_ED_pcs_0_ts_80_t_0.800000.vtu sigma sigma 2e-6 0
    out_beam_fine_ED_pcs_0_ts_80_t_0.800000.vtu out_beam_fine_ED_pcs_0_ts_80_t_0.800000.vtu epsilon epsilon 2e-11 0

    out_beam_fine_ED_pcs_0_ts_100_t_1.000000.vtu out_beam_fine_ED_pcs_0_ts_100_t_1.000000.vtu displacement displacement 2e-11 1e-13
    out_beam_fine_ED_pcs_0_ts_100_t_1.000000.vtu out_beam_fine_ED_pcs_0_ts_100_t_1.000000.vtu damage damage 2e-11 0
    out_beam_fine_ED_pcs_0_ts_100_t_1.000000.vtu out_beam_fine_ED_pcs_0_ts_100_t_1.000000.vtu sigma sigma 2e-6 0
    out_beam_fine_ED_pcs_0_ts_100_t_1.000000.vtu out_beam_fine_ED_pcs_0_ts_100_t_1.000000.vtu epsilon epsilon 2e-11 0
)

AddTest(
    NAME Mechanics_PlasticDamageNonlocalModel_NNLED_Ehlers_Damage_BeamCoarseQuad
    PATH Mechanics/EhlersDamageNonLocal/beam/coarse_quad
    EXECUTABLE ogs
    EXECUTABLE_ARGS beam.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu displacement displacement 2e-11 1e-13
    out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu damage damage 2e-11 0
    out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu sigma sigma 2e-6 0
    out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu epsilon epsilon 2e-11 0

    out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu displacement displacement 2e-11 1e-13
    out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu damage damage 2e-11 0
    out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu sigma sigma 2e-6 0
    out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu epsilon epsilon 2e-11 0

    out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu displacement displacement 2e-11 1e-13
    out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu damage damage 2e-11 0
    out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu sigma sigma 3e-6 0
    out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu epsilon epsilon 2e-11 0

    out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu displacement displacement 2e-11 1e-13
    out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu damage damage 2e-11 0
    out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu sigma sigma 4e-6 0
    out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu epsilon epsilon 2e-11 0
)

AddTest(
    NAME Mechanics_PlasticDamageNonlocalModel_NNLED_Ehlers_Damage_BeamCoarseTria
    PATH Mechanics/EhlersDamageNonLocal/beam/coarse_tria
    EXECUTABLE ogs
    EXECUTABLE_ARGS beam.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu displacement displacement 2e-11 1e-13
    out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu damage damage 2e-11 0
    out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu sigma sigma 2e-7 0
    out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu epsilon epsilon 2e-11 0

    out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu displacement displacement 2e-11 1e-13
    out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu damage damage 2e-11 0
    out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu sigma sigma 5e-7 0
    out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu epsilon epsilon 2e-11 0

    out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu displacement displacement 2e-11 1e-13
    out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu damage damage 2e-11 0
    out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu sigma sigma 4e-7 0
    out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu epsilon epsilon 2e-11 0

    out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu displacement displacement 2e-11 1e-13
    out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu damage damage 2e-11 0
    out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu sigma sigma 6e-7 0
    out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu epsilon epsilon 2e-11 0
)

AddTest(
    NAME Mechanics_PlasticDamageNonlocalModel_NNLED_Ehlers_Damage_BeamCoarse_3DHexa
    PATH Mechanics/EhlersDamageNonLocal/beam/3D_hexa
    EXECUTABLE ogs
    EXECUTABLE_ARGS beam.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu displacement displacement 2e-11 1e-13
    out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu damage damage 2e-11 0
    out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu sigma sigma 3e-7 0
    out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu epsilon epsilon 2e-11 0

    out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu displacement displacement 2e-11 1e-13
    out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu damage damage 2e-11 0
    out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu sigma sigma 4e-7 0
    out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu out_beam_coarse_ED_pcs_0_ts_60_t_0.600000.vtu epsilon epsilon 2e-11 0

    out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu displacement displacement 2e-11 1e-13
    out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu damage damage 2e-11 0
    out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu sigma sigma 5e-7 0
    out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu out_beam_coarse_ED_pcs_0_ts_80_t_0.800000.vtu epsilon epsilon 2e-11 0

    out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu displacement displacement 2e-11 1e-13
    out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu damage damage 2e-11 0
    out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu sigma sigma 5e-7 0
    out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu out_beam_coarse_ED_pcs_0_ts_100_t_1.000000.vtu epsilon epsilon 2e-11 0
)

AddTest(
    NAME LARGE_Mechanics_PlasticDamageNonlocalModel_NNLED_Ehlers_Damage_BeamCoarse_3DTetra
    PATH Mechanics/EhlersDamageNonLocal/beam/3D_tetr
    EXECUTABLE ogs
    EXECUTABLE_ARGS beam.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    out_beam_coarse_ED_pcs_0_ts_10_t_0.100000.vtu out_beam_coarse_ED_pcs_0_ts_10_t_0.100000.vtu displacement displacement 2e-11 1e-13
    out_beam_coarse_ED_pcs_0_ts_10_t_0.100000.vtu out_beam_coarse_ED_pcs_0_ts_10_t_0.100000.vtu damage damage 2e-11 0
    out_beam_coarse_ED_pcs_0_ts_10_t_0.100000.vtu out_beam_coarse_ED_pcs_0_ts_10_t_0.100000.vtu sigma sigma 5e-8 0
    out_beam_coarse_ED_pcs_0_ts_10_t_0.100000.vtu out_beam_coarse_ED_pcs_0_ts_10_t_0.100000.vtu epsilon epsilon 2e-11 0

    out_beam_coarse_ED_pcs_0_ts_20_t_0.200000.vtu out_beam_coarse_ED_pcs_0_ts_20_t_0.200000.vtu displacement displacement 2e-11 1e-13
    out_beam_coarse_ED_pcs_0_ts_20_t_0.200000.vtu out_beam_coarse_ED_pcs_0_ts_20_t_0.200000.vtu damage damage 2e-11 0
    out_beam_coarse_ED_pcs_0_ts_20_t_0.200000.vtu out_beam_coarse_ED_pcs_0_ts_20_t_0.200000.vtu sigma sigma 7e-8 0
    out_beam_coarse_ED_pcs_0_ts_20_t_0.200000.vtu out_beam_coarse_ED_pcs_0_ts_20_t_0.200000.vtu epsilon epsilon 2e-11 0

    out_beam_coarse_ED_pcs_0_ts_30_t_0.300000.vtu out_beam_coarse_ED_pcs_0_ts_30_t_0.300000.vtu displacement displacement 2e-11 1e-13
    out_beam_coarse_ED_pcs_0_ts_30_t_0.300000.vtu out_beam_coarse_ED_pcs_0_ts_30_t_0.300000.vtu damage damage 2e-11 0
    out_beam_coarse_ED_pcs_0_ts_30_t_0.300000.vtu out_beam_coarse_ED_pcs_0_ts_30_t_0.300000.vtu sigma sigma 2e-7 0
    out_beam_coarse_ED_pcs_0_ts_30_t_0.300000.vtu out_beam_coarse_ED_pcs_0_ts_30_t_0.300000.vtu epsilon epsilon 2e-11 0

    out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu displacement displacement 2e-11 1e-13
    out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu damage damage 2e-11 0
    out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu sigma sigma 2e-7 0
    out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu out_beam_coarse_ED_pcs_0_ts_40_t_0.400000.vtu epsilon epsilon 2e-11 0
)

AddTest(
    NAME Mechanics_PlasticDamageNonlocalModel_NNLED_Ehlers_Damage_HoledBeam
    PATH Mechanics/EhlersDamageNonLocal/holed_beam
    EXECUTABLE ogs
    EXECUTABLE_ARGS holed_beam.prj
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    out_beam_hole_ED_pcs_0_ts_100_t_0.250000.vtu out_beam_hole_ED_pcs_0_ts_100_t_0.250000.vtu displacement displacement 2e-11 1e-13
    out_beam_hole_ED_pcs_0_ts_100_t_0.250000.vtu out_beam_hole_ED_pcs_0_ts_100_t_0.250000.vtu damage damage 2e-11 0
    out_beam_hole_ED_pcs_0_ts_100_t_0.250000.vtu out_beam_hole_ED_pcs_0_ts_100_t_0.250000.vtu sigma sigma 3e-7 0
    out_beam_hole_ED_pcs_0_ts_100_t_0.250000.vtu out_beam_hole_ED_pcs_0_ts_100_t_0.250000.vtu epsilon epsilon 2e-11 0

    out_beam_hole_ED_pcs_0_ts_200_t_0.500000.vtu out_beam_hole_ED_pcs_0_ts_200_t_0.500000.vtu displacement displacement 2e-11 1e-13
    out_beam_hole_ED_pcs_0_ts_200_t_0.500000.vtu out_beam_hole_ED_pcs_0_ts_200_t_0.500000.vtu damage damage 2e-11 0
    out_beam_hole_ED_pcs_0_ts_200_t_0.500000.vtu out_beam_hole_ED_pcs_0_ts_200_t_0.500000.vtu sigma sigma 5e-7 0
    out_beam_hole_ED_pcs_0_ts_200_t_0.500000.vtu out_beam_hole_ED_pcs_0_ts_200_t_0.500000.vtu epsilon epsilon 2e-11 0

    out_beam_hole_ED_pcs_0_ts_300_t_0.750000.vtu out_beam_hole_ED_pcs_0_ts_300_t_0.750000.vtu displacement displacement 2e-11 1e-13
    out_beam_hole_ED_pcs_0_ts_300_t_0.750000.vtu out_beam_hole_ED_pcs_0_ts_300_t_0.750000.vtu damage damage 2e-11 0
    out_beam_hole_ED_pcs_0_ts_300_t_0.750000.vtu out_beam_hole_ED_pcs_0_ts_300_t_0.750000.vtu sigma sigma 8e-7 0
    out_beam_hole_ED_pcs_0_ts_300_t_0.750000.vtu out_beam_hole_ED_pcs_0_ts_300_t_0.750000.vtu epsilon epsilon 2e-11 0

    out_beam_hole_ED_pcs_0_ts_400_t_1.000000.vtu out_beam_hole_ED_pcs_0_ts_400_t_1.000000.vtu displacement displacement 2e-11 1e-13
    out_beam_hole_ED_pcs_0_ts_400_t_1.000000.vtu out_beam_hole_ED_pcs_0_ts_400_t_1.000000.vtu damage damage 2e-11 0
    out_beam_hole_ED_pcs_0_ts_400_t_1.000000.vtu out_beam_hole_ED_pcs_0_ts_400_t_1.000000.vtu sigma sigma 2e-6 0
    out_beam_hole_ED_pcs_0_ts_400_t_1.000000.vtu out_beam_hole_ED_pcs_0_ts_400_t_1.000000.vtu epsilon epsilon 2e-11 0
)
