
AddTest(
    NAME MapGeometryToMeshSurface_Ammer
    PATH MeshGeoToolsLib/Ammer/
    EXECUTABLE MapGeometryToMeshSurface
    EXECUTABLE_ARGS -m Ammer-Homogen100m-Final-TopSurface.vtu -i Ammer-Rivers.gml -o ${CMAKE_BINARY_DIR}/Tests/Data/MeshGeoToolsLib/Ammer/Ammer-Rivers-Mapped.gml
    TESTER diff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA Ammer-Rivers-Mapped.gml
)

AddTest(
    NAME LARGE_MapGeometryToMeshSurface_Bode
    PATH MeshGeoToolsLib/Bode/
    EXECUTABLE MapGeometryToMeshSurface
    EXECUTABLE_ARGS -m BodeComplex.msh -i BodeEZG_Fliessgewaesser.gml -o ${CMAKE_BINARY_DIR}/Tests/Data/MeshGeoToolsLib/Bode/BodeEZG_Fliessgewaesser-Mapped.gml
    REQUIREMENTS NOT OGS_USE_MPI
    TESTER diff
    DIFF_DATA BodeEZG_Fliessgewaesser-Mapped.gml
)

AddTest(
    NAME LARGE_MapGeometryToMeshSurface_Naegelstedt
    PATH MeshGeoToolsLib/Naegelstedt
    EXECUTABLE MapGeometryToMeshSurface
    EXECUTABLE_ARGS -m SmallTest.vtu -i RiverNetwork.gml -o ${CMAKE_BINARY_DIR}/Tests/Data/MeshGeoToolsLib/Naegelstedt/RiverNetwork-Mapped.gml
    REQUIREMENTS NOT OGS_USE_MPI
    TESTER diff
    DIFF_DATA RiverNetwork-Mapped.gml
)

AddTest(
    NAME postLIE
    PATH LIE/PostProcessing
    EXECUTABLE postLIE
    EXECUTABLE_ARGS -i single_joint_pcs_0.pvd -o ${CMAKE_BINARY_DIR}/Tests/Data/LIE/PostProcessing/post_single_joint_pcs_0.pvd
    REQUIREMENTS NOT OGS_USE_MPI
    ABSTOL 1e-14 RELTOL 1e-14
    TESTER vtkdiff
    DIFF_DATA
    expected_post_single_joint_pcs_0_ts_1_t_1.000000.vtu post_single_joint_pcs_0_ts_1_t_1.000000.vtu u u
)
