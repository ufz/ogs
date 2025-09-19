/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>

#include <array>
#include <string>

#include "BaseLib/Logging.h"
#include "BaseLib/MPI.h"
#include "BaseLib/RunTime.h"
#include "BaseLib/TCLAPArguments.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/Utils/addPropertyToMesh.h"
#include "MeshToolsLib/MeshQuality/ElementQualityInterface.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Add element quality as a mesh property.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2025, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::ValueArg<std::string> mesh_out_arg("o", "output_mesh_file",
                                              "Output (.vtu | .msh) mesh file ",
                                              true, "", "OUTPUT_FILE");
    cmd.add(mesh_out_arg);
    std::vector<std::string> allowed_element_criterions{
        "ElementSize", "EdgeRatio", "EquiAngleSkew", "RadiusEdgeRatio",
        "SizeDifference"};
    TCLAP::ValuesConstraint<std::string> element_criterions{
        allowed_element_criterions};
    TCLAP::ValueArg<std::string> criterion_arg{
        "c", "quality_criterion", "quality criterion", true,
        "",  &element_criterions};
    cmd.add(criterion_arg);
    TCLAP::ValueArg<std::string> mesh_in_arg("i", "input_mesh_file",
                                             "Input (.vtu | .msh) mesh file",
                                             true, "", "INPUT_FILE");
    cmd.add(mesh_in_arg);
    auto log_level_arg = BaseLib::makeLogLevelArg();
    cmd.add(log_level_arg);
    cmd.parse(argc, argv);

    BaseLib::MPI::Setup mpi_setup(argc, argv);
    BaseLib::initOGSLogger(log_level_arg.getValue());

    // read the mesh file
    BaseLib::RunTime run_time;
    run_time.start();
    std::unique_ptr<MeshLib::Mesh> mesh(MeshLib::IO::readMeshFromFile(
        mesh_in_arg.getValue(), true /* compute_element_neighbors */));
    if (!mesh)
    {
        return EXIT_FAILURE;
    }
    INFO("Time for reading: {:g} s", run_time.elapsed());

    // Geometric information
    MeshLib::MeshQualityType const type =
        MeshLib::String2MeshQualityType(criterion_arg.getValue());
    MeshToolsLib::ElementQualityInterface element_quality(*mesh, type);
    auto const element_quality_vector = element_quality.getQualityVector();
    MeshLib::addPropertyToMesh<double>(*mesh, criterion_arg.getValue(),
                                       MeshLib::MeshItemType::Cell, 1,
                                       {element_quality_vector});
    INFO("Writing mesh '{:s}' ... ", mesh_out_arg.getValue());
    MeshLib::IO::writeMeshToFile(*mesh, mesh_out_arg.getValue());
    INFO("done.");
    return EXIT_SUCCESS;
}
