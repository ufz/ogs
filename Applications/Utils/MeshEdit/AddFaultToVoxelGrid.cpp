/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

// ThirdParty
#include <tclap/CmdLine.h>

#ifdef USE_PETSC
#include <mpi.h>
#endif

#include "InfoLib/GitInfo.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/Mesh.h"
#include "MeshToolsLib/MeshGenerators/AddFaultToVoxelGrid.h"

int main(int argc, char* argv[])
{
    constexpr int mat_not_set = std::numeric_limits<int>::max();

    TCLAP::CmdLine cmd(
        "Marks all elements in a voxel grid (i.e. a structured hex grid, for "
        "instance created with Layers2Grid or Vtu2Grid) that are intersected "
        "by a triangulated 2D mesh representing a fault or some other "
        "significant structure. The material group for those intersected "
        "elements can be explicitly specified, otherwise the largest existing "
        "MaterialID will be increased by one.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2024, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<int> id_arg("m", "material",
                                "material id for cells intersected by fault",
                                false, mat_not_set, "non-negative integer");
    cmd.add(id_arg);

    TCLAP::ValueArg<std::string> output_arg(
        "o", "output", "name of output mesh (*.vtu)", true, "", "string");
    cmd.add(output_arg);

    TCLAP::ValueArg<std::string> fault_arg(
        "f", "fault", "name of mesh representing fault (*.vtu)", true, "",
        "string");
    cmd.add(fault_arg);

    TCLAP::ValueArg<std::string> input_arg(
        "i", "input",
        "name of the input file list containing the paths the all input layers "
        "in correct order from top to bottom",
        true, "", "string");
    cmd.add(input_arg);
    cmd.parse(argc, argv);

#ifdef USE_PETSC
    MPI_Init(&argc, &argv);
#endif

    std::string const input_name = input_arg.getValue();
    std::string const fault_name = fault_arg.getValue();
    std::string const output_name = output_arg.getValue();

    using namespace MeshToolsLib::MeshGenerator::AddFaultToVoxelGrid;

    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshLib::IO::readMeshFromFile(input_name));
    std::unique_ptr<MeshLib::Mesh> fault(
        MeshLib::IO::readMeshFromFile(fault_name));
    if (mesh == nullptr)
    {
        ERR("Input mesh not found...");
#ifdef USE_PETSC
        MPI_Finalize();
#endif
        return EXIT_FAILURE;
    }
    auto const& mat_ids = MeshLib::materialIDs(*mesh);
    if (!mat_ids)
    {
#ifdef USE_PETSC
        MPI_Finalize();
#endif
        ERR("Input mesh has no material IDs");
        return EXIT_FAILURE;
    }
    int fault_id = id_arg.getValue();
    if (!id_arg.isSet())
    {
        auto it = std::max_element(mat_ids->cbegin(), mat_ids->cend());
        fault_id = *it + 1;
    }
    if (addFaultToVoxelGrid(mesh.get(), fault.get(), fault_id))
    {
        MeshLib::IO::VtuInterface vtu(mesh.get());
        vtu.writeToFile(output_name);
#ifdef USE_PETSC
        MPI_Finalize();
#endif
        INFO("The fault was successfully added.");
        return EXIT_SUCCESS;
    }
#ifdef USE_PETSC
    MPI_Finalize();
#endif
    ERR("No fault could be added.");
    return EXIT_FAILURE;
}
