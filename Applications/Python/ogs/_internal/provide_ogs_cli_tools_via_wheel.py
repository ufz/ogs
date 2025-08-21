import os
import platform
import subprocess
import sys

from . import OGS_USE_PATH
from .get_bin_dir import get_bin_dir

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
    "CreateAnchors",
    "CreateBoundaryConditionsAlongPolylines",
    "createIntermediateRasters",
    "createLayeredMeshFromRasters",
    "createMeshElemPropertiesFromASCRaster",
    "createNeumannBc",
    "createQuadraticMesh",
    "createRaster",
    "createTetgenSmeshFromRasters",
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
    "ipDataToPointCloud",
    "Layers2Grid",
    "MapGeometryToMeshSurface",
    "mergeMeshToBulkMesh",
    "Mesh2Raster",
    "MeshMapping",
    "MoveGeometry",
    "MoveMesh",
    "moveMeshNodes",
    "mpmetis",
    "NodeReordering",
    "ogs",
    "OGS2VTK",
    "partmesh",
    "postLIE",
    "PVD2XDMF",
    "pvtu2vtu",
    "queryMesh",
    "Raster2ASC",
    "Raster2Mesh",
    "RemoveGhostData",
    "removeMeshElements",
    "RemoveUnusedPoints",
    "ReorderMesh",
    "ResetPropertiesInPolygonalRegion",
    "reviseMesh",
    "scaleProperty",
    "swapNodeCoordinateAxes",
    "TecPlotTools",
    "TIN2VTK",
    "VTK2OGS",
    "VTK2TIN",
    "vtkdiff",
    "Vtu2Grid",
    "xdmfdiff",
]


def pyproject_get_binaries():
    return {
        binary: f"ogs._internal.provide_ogs_cli_tools_via_wheel:{binary}"
        for binary in binaries_list
    }


# Not used when OGS_USE_PATH is true!
def ogs():
    raise SystemExit(ogs_with_args(sys.argv))


def ogs_with_args(argv):
    import ogs.simulator as sim  # noqa: PLC0415

    return_code = sim.initialize(argv)

    # map mangled TCLAP status to usual exit status
    if return_code == 3:  # EXIT_ARGPARSE_FAILURE
        sim.finalize()
        return 1  # EXIT_FAILURE
    if return_code == 2:  # EXIT_ARGPARSE_EXIT_OK
        sim.finalize()
        return 0  # EXIT_SUCCESS

    if return_code != 0:
        sim.finalize()
        return return_code

    return_code = sim.executeSimulation()
    sim.finalize()
    return return_code


if "PEP517_BUILD_BACKEND" not in os.environ:
    OGS_BIN_DIR = get_bin_dir()

    if platform.system() == "Windows":
        os.add_dll_directory(OGS_BIN_DIR)

    def _program(name, args):
        exe = OGS_BIN_DIR / name
        env = None  # by default use unmodified environment
        if OGS_USE_PATH:
            exe = name
            env = os.environ.copy()
            # prevent infinite recursion if OGS in PATH happens to be this very
            # script
            env["OGS_USE_PATH"] = "0"
            print(f"OGS_USE_PATH is true: {name} from $PATH is used!")
        return subprocess.run([exe] + args, env=env).returncode  # noqa: PLW1510

    FUNC_TEMPLATE = """def {0}(): raise SystemExit(_program("{0}", sys.argv[1:]))"""
    for f in binaries_list:
        if f == "ogs" and not OGS_USE_PATH:
            continue  # provided by separate function
        # When OGS_USE_PATH is true then ogs()-function above is not used!
        exec(FUNC_TEMPLATE.format(f))
