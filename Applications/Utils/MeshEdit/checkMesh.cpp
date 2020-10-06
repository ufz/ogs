/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <array>
#include <string>

#include <tclap/CmdLine.h>

#include "InfoLib/GitInfo.h"
#include "BaseLib/StringTools.h"
#include "BaseLib/MemWatch.h"
#include "BaseLib/RunTime.h"
#include "BaseLib/FileTools.h"

#include "GeoLib/AABB.h"

#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshInformation.h"
#include "MeshLib/MeshQuality/MeshValidation.h"

#include "MeshLib/IO/readMeshFromFile.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Checks mesh properties.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2020, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::UnlabeledValueArg<std::string> mesh_arg(
        "mesh-file", "input mesh file", true, "", "string");
    cmd.add(mesh_arg);
    TCLAP::SwitchArg valid_arg("v", "validation", "validate the mesh");
    cmd.add(valid_arg);
    TCLAP::SwitchArg print_properties_arg(
        "p", "print_properties", "print properties stored in the mesh");
    cmd.add(print_properties_arg);

    cmd.parse(argc, argv);

    // read the mesh file
    BaseLib::MemWatch mem_watch;
    const unsigned long mem_without_mesh(mem_watch.getVirtMemUsage());
    BaseLib::RunTime run_time;
    run_time.start();
    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshLib::IO::readMeshFromFile(mesh_arg.getValue()));
    if (!mesh)
    {
        return EXIT_FAILURE;
    }

    const unsigned long mem_with_mesh(mem_watch.getVirtMemUsage());
    if (mem_with_mesh > 0)
    {
        INFO("Memory size: {:i} MB",
             (mem_with_mesh - mem_without_mesh) / (1024 * 1024));
        (void)mem_with_mesh;
    }
    INFO("Time for reading: {:g} s", run_time.elapsed());

    // Geometric information
    const GeoLib::AABB aabb(MeshLib::MeshInformation::getBoundingBox(*mesh));
    auto minPt(aabb.getMinPoint());
    auto maxPt(aabb.getMaxPoint());
    INFO("Node coordinates:");
    INFO("\tx [{:g}, {:g}] (extent {:g})", minPt[0], maxPt[0],
         maxPt[0] - minPt[0]);
    INFO("\ty [{:g}, {:g}] (extent {:g})", minPt[1], maxPt[1],
         maxPt[1] - minPt[1]);
    INFO("\tz [{:g}, {:g}] (extent {:g})", minPt[2], maxPt[2],
         maxPt[2] - minPt[2]);

    INFO("Edge length: [{:g}, {:g}]", mesh->getMinEdgeLength(),
         mesh->getMaxEdgeLength());

    // Element information

    MeshLib::MeshInformation::writeAllNumbersOfElementTypes(*mesh);

    MeshLib::MeshInformation::writePropertyVectorInformation(*mesh);

    if (valid_arg.isSet())
    {
        // MeshValidation outputs error messages
        // Remark: MeshValidation can modify the original mesh
        MeshLib::MeshInformation::writeMeshValidationResults(*mesh);
    }
}
