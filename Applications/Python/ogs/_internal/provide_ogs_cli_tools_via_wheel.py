import os
import platform
import subprocess
import sys

binaries_list = [
    "addDataToRaster",
    "AddElementQuality",
    "AddFaultToVoxelGrid",
    "AddLayer",
    "appendLinesAlongPolyline",
    "AssignRasterDataToMesh",
    "checkMesh",
    "ComputeNodeAreasFromSurfaceMesh",
    "computeSurfaceNodeIDsInPolygonalRegion",
    "constructMeshesFromGeometry",
    "convertGEO",
    "convertToLinearMesh",
    "convertVtkDataArrayToVtkDataArray",
    "CreateBoundaryConditionsAlongPolylines",
    "createIntermediateRasters",
    "createLayeredMeshFromRasters",
    "createMeshElemPropertiesFromASCRaster",
    "createNeumannBc",
    "createQuadraticMesh",
    "createRaster",
    "editMaterialID",
    "ExtractBoundary",
    "ExtractMaterials",
    "ExtractSurface",
    "generateGeometry",
    "generateMatPropsFromMatID",
    "generateStructuredMesh",
    "geometryToGmshGeo",
    "GMSH2OGS",
    "GocadSGridReader",
    "GocadTSurfaceReader",
    "identifySubdomains",
    "IntegrateBoreholesIntoMesh",
    "Layers2Grid",
    "MapGeometryToMeshSurface",
    "Mesh2Raster",
    "MoveGeometry",
    "MoveMesh",
    "moveMeshNodes",
    "mpmetis",
    "NodeReordering",
    "ogs",
    "OGS2VTK",
    "partmesh",
    "PVD2XDMF",
    "queryMesh",
    "Raster2Mesh",
    "RemoveGhostData",
    "removeMeshElements",
    "ResetPropertiesInPolygonalRegion",
    "reviseMesh",
    "scaleProperty",
    "swapNodeCoordinateAxes",
    "TecPlotTools",
    "tetgen",
    "TIN2VTK",
    "VTK2OGS",
    "VTK2TIN",
    "vtkdiff",
    "Vtu2Grid",
]

if "PEP517_BUILD_BACKEND" not in os.environ:
    if platform.system() == "Windows":
        os.add_dll_directory(os.path.join(os.path.dirname(__file__), "bin"))

    OGS_BIN_DIR = os.path.join(os.path.join(os.path.dirname(__file__), "bin"))

    def _program(name, args):
        return subprocess.run([os.path.join(OGS_BIN_DIR, name)] + args).returncode

    FUNC_TEMPLATE = """def {0}(): raise SystemExit(_program("{0}", sys.argv[1:]))"""
    for f in binaries_list:
        exec(FUNC_TEMPLATE.format(f))
