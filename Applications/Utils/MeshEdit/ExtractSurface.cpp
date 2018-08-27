/**
 * @brief Extracts the surface from the given mesh.
 *
 * @copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "BaseLib/BuildInfo.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"

#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"

#include "MathLib/Vector3.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSurfaceExtraction.h"

int main (int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd(
        "Tool extracts the surface of the given mesh. The documentation is "
        "available at "
        "https://docs.opengeosys.org/docs/tools/meshing/extract-surface.\n\n"
        "OpenGeoSys-6 software, version " +
            BaseLib::BuildInfo::git_describe +
            ".\n"
            "Copyright (c) 2012-2018, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', BaseLib::BuildInfo::git_describe);
    TCLAP::ValueArg<std::string> mesh_in(
        "i", "mesh-input-file",
        "the name of the file containing the input mesh", true, "",
        "file name of input mesh");
    cmd.add(mesh_in);
    TCLAP::ValueArg<std::string> mesh_out(
        "o", "mesh-output-file",
        "the name of the file the surface mesh should be written to", false, "",
        "file name of output mesh");
    cmd.add(mesh_out);
    TCLAP::ValueArg<double> x("x", "x-component", "x component of the normal",
                              false, 0, "floating point value");
    cmd.add(x);
    TCLAP::ValueArg<double> y("y", "y-component", "y component of the normal",
                              false, 0, "floating point value");
    cmd.add(y);
    TCLAP::ValueArg<double> z("z", "z-component", "z component of the normal",
                              false, -1.0, "floating point value");
    cmd.add(z);

    TCLAP::ValueArg<std::string> node_prop_name(
        "n", "node-property-name",
        "the name of the data array the subsurface/bulk node id will be stored "
        "to",
        false, "bulk_node_ids", "string");
    cmd.add(node_prop_name);
    TCLAP::ValueArg<std::string> element_prop_name(
        "e", "element-property-name",
        "the name of the data array the subsurface/bulk element id will be "
        "stored to",
        false, "OriginalSubsurfaceElementIDs", "string");
    cmd.add(element_prop_name);
    TCLAP::ValueArg<std::string> face_prop_name(
        "f", "face-property-name",
        "the name of the data array the surface face id of the subsurface/bulk "
        "element will be stored to",
        false, "OriginalFaceIDs", "string");
    cmd.add(face_prop_name);

    TCLAP::ValueArg<double> angle_arg(
        "a", "angle", "angle between given normal and element normal", false,
        90, "floating point value");
    cmd.add(angle_arg);

    cmd.parse(argc, argv);

    std::unique_ptr<MeshLib::Mesh const> mesh(
        MeshLib::IO::readMeshFromFile(mesh_in.getValue()));
    INFO("Mesh read: %u nodes, %u elements.", mesh->getNumberOfNodes(),
         mesh->getNumberOfElements());

    // extract surface
    MathLib::Vector3 const dir(x.getValue(), y.getValue(), z.getValue());
    double const angle(angle_arg.getValue());
    std::unique_ptr<MeshLib::Mesh> surface_mesh(
        MeshLib::MeshSurfaceExtraction::getMeshSurface(
            *mesh, dir, angle, node_prop_name.getValue(),
            element_prop_name.getValue(), face_prop_name.getValue()));

    std::string out_fname(mesh_out.getValue());
    if (out_fname.empty())
        out_fname = BaseLib::dropFileExtension(mesh_in.getValue()) + "_sfc.vtu";
    MeshLib::IO::writeMeshToFile(*surface_mesh, out_fname);

    return EXIT_SUCCESS;
}
