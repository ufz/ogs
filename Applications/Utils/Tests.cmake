
AddTest(
    NAME MapGeometryToMeshSurface_Ammer
    PATH MeshGeoToolsLib/Ammer/
    EXECUTABLE MapGeometryToMeshSurface
    EXECUTABLE_ARGS -m Ammer-Homogen100m-Final-TopSurface.vtu -i Ammer-Rivers.gml -o ${CMAKE_BINARY_DIR}/Tests/Data/MeshGeoToolsLib/Ammer/Ammer-Rivers-Mapped.gml
    WRAPPER time
    TESTER diff
    DIFF_DATA Ammer-Rivers-Mapped.gml
)

AddTest(
    NAME LARGE_MapGeometryToMeshSurface_Bode
    PATH MeshGeoToolsLib/Bode/
    EXECUTABLE MapGeometryToMeshSurface
    EXECUTABLE_ARGS -m BodeComplex.msh -i BodeEZG_Fliessgewaesser.gml -o ${CMAKE_BINARY_DIR}/Tests/Data/MeshGeoToolsLib/Bode/BodeEZG_Fliessgewaesser-Mapped.gml
    WRAPPER time
    TESTER diff
    DIFF_DATA BodeEZG_Fliessgewaesser-Mapped.gml
)

AddTest(
    NAME LARGE_MapGeometryToMeshSurface_Naegelstedt
    PATH MeshGeoToolsLib/Naegelstedt
    EXECUTABLE MapGeometryToMeshSurface
    EXECUTABLE_ARGS -m SmallTest.vtu -i RiverNetwork.gml -o ${CMAKE_BINARY_DIR}/Tests/Data/MeshGeoToolsLib/Naegelstedt/RiverNetwork-Mapped.gml
    WRAPPER time
    TESTER diff
    DIFF_DATA RiverNetwork-Mapped.gml
)
