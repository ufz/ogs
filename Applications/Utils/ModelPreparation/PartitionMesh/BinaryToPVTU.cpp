// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <tclap/CmdLine.h>
#include <vtkMPIController.h>
#include <vtkSmartPointer.h>

#include "BaseLib/CPUTime.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/MPI.h"
#include "BaseLib/RunTime.h"
#include "BaseLib/TCLAPArguments.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/NodePartitionedMesh.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Reads the binary mesh format and writes the data as vtus/pvtu.\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2025, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<std::string> mesh_input(
        "i", "mesh-input-file-base",
        "Input (.bin). The base name of the files containing the input mesh, "
        "i.e., the file "
        "name without the string beginning with '_partitioned' and ending "
        "with ",
        true, "", "BASE_FILENAME_INPUT_MESH");
    cmd.add(mesh_input);

    TCLAP::ValueArg<std::string> output_directory_arg(
        "o", "output",
        "output directory name and output base file name without extensions",
        true, "", "directory/base_file_name_without_extensions");
    cmd.add(output_directory_arg);

    auto log_level_arg = BaseLib::makeLogLevelArg();
    cmd.add(log_level_arg);

    cmd.parse(argc, argv);

    BaseLib::MPI::Setup mpi_setup(argc, argv);
    BaseLib::initOGSLogger(log_level_arg.getValue());

    // start the timer
    BaseLib::RunTime run_timer;
    run_timer.start();
    BaseLib::CPUTime CPU_timer;
    CPU_timer.start();

    // init vtkMPI
    vtkSmartPointer<vtkMPIController> controller =
        vtkSmartPointer<vtkMPIController>::New();
    controller->Initialize(&argc, &argv, 1);
    vtkMPIController::SetGlobalController(controller);

    MeshLib::NodePartitionedMesh* mesh =
        dynamic_cast<MeshLib::NodePartitionedMesh*>(
            MeshLib::IO::readMeshFromFile(mesh_input.getValue()));
    INFO("Mesh '{:s}' read: {:d} nodes, {:d} elements.",
         mesh->getName(),
         mesh->getNumberOfNodes(),
         mesh->getNumberOfElements());
    INFO("Mesh read runtime: {:g} s.", run_timer.elapsed());
    INFO("Mesh read CPU time: {:g} s.", CPU_timer.elapsed());

    // Write output file
    auto const file_name = output_directory_arg.getValue() + ".vtu";
    DBUG("Writing output to '{:s}'.", file_name);

    MeshLib::IO::VtuInterface vtu_interface(mesh);
    vtu_interface.writeToFile(file_name);

    INFO("Total runtime: {:g} s.", run_timer.elapsed());
    INFO("Total CPU time: {:g} s.", CPU_timer.elapsed());
    return EXIT_SUCCESS;
}
