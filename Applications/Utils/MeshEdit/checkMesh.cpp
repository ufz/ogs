// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <spdlog/fmt/ranges.h>
#include <tclap/CmdLine.h>

#include <array>
#include <string>

#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/MPI.h"
#include "BaseLib/MemWatch.h"
#include "BaseLib/RunTime.h"
#include "BaseLib/StringTools.h"
#include "BaseLib/TCLAPArguments.h"
#include "GeoLib/AABB.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshToolsLib/MeshInformation.h"
#include "MeshToolsLib/MeshQuality/MeshValidation.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Checks mesh properties.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2025, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::UnlabeledValueArg<std::string> mesh_arg(
        "mesh-file", "Input (.vtu | .msh). Mesh file", true, "", "INPUT_FILE");
    cmd.add(mesh_arg);
    TCLAP::SwitchArg valid_arg("v", "validation", "validate the mesh");
    cmd.add(valid_arg);
    TCLAP::SwitchArg print_properties_arg(
        "p", "print_properties", "print properties stored in the mesh");
    cmd.add(print_properties_arg);
    auto log_level_arg = BaseLib::makeLogLevelArg();
    cmd.add(log_level_arg);

    cmd.parse(argc, argv);

    BaseLib::MPI::Setup mpi_setup(argc, argv);
    BaseLib::initOGSLogger(log_level_arg.getValue());

    // read the mesh file
    BaseLib::MemWatch mem_watch;
    const unsigned long mem_without_mesh(mem_watch.getVirtMemUsage());
    BaseLib::RunTime run_time;
    run_time.start();
    std::unique_ptr<MeshLib::Mesh> mesh(MeshLib::IO::readMeshFromFile(
        mesh_arg.getValue(), true /* compute_element_neighbors */));
    if (!mesh)
    {
        return EXIT_FAILURE;
    }

    const unsigned long mem_with_mesh(mem_watch.getVirtMemUsage());
    if (mem_with_mesh > 0)
    {
        INFO("Memory size: {} MiB",
             (mem_with_mesh - mem_without_mesh) / (1024 * 1024));
        (void)mem_with_mesh;
    }
    INFO("Time for reading: {:g} s", run_time.elapsed());

    // Geometric information
    const GeoLib::AABB aabb(
        MeshToolsLib::MeshInformation::getBoundingBox(*mesh));
    INFO("Axis aligned bounding box: {}", aabb);

    auto const [min, max] = minMaxEdgeLength(mesh->getElements());
    INFO("Min/max edge lengths: [{:g}, {:g}]", min, max);

    // Element information

    MeshToolsLib::MeshInformation::writeAllNumbersOfElementTypes(*mesh);

    if (print_properties_arg.isSet())
    {
        MeshToolsLib::MeshInformation::writePropertyVectorInformation(*mesh);
        INFO("MaterialID-list: [{}]",
             fmt::join(MeshToolsLib::MeshInformation::getMaterialIDs(*mesh),
                       ", "));
    }

    if (valid_arg.isSet())
    {
        // MeshValidation outputs error messages
        // Remark: MeshValidation can modify the original mesh
        MeshToolsLib::MeshInformation::writeMeshValidationResults(*mesh);
    }
    return EXIT_SUCCESS;
}
