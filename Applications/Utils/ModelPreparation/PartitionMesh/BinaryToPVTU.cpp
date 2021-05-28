/*!
  \file
  \brief  A tool for debugging the mesh partitioning.

  \copyright
  Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/

#include <spdlog/spdlog.h>
#include <tclap/CmdLine.h>
#include <mpi.h>

#include <vtkMPIController.h>
#include <vtkSmartPointer.h>

#include "BaseLib/CPUTime.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/RunTime.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/IO/readMeshFromFile.h"

#include "MeshLib/IO/MPI_IO/NodePartitionedMeshReader.h"
#include "MeshLib/NodePartitionedMesh.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Reads the binary mesh format and writes the data as vtus/pvtu.\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2021, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<std::string> mesh_input(
        "i", "mesh-input-file-base",
        "the base name of the files containing the input mesh, i.e., the file "
        "name without the string beginning with '_partitioned' and ending with "
        "'.bin'",
        true, "", "base_file_name_of_input_mesh");
    cmd.add(mesh_input);

    TCLAP::ValueArg<std::string> output_directory_arg(
        "o", "output",
        "output directory name and output base file name without extensions",
        true, "", "directory/base_file_name_without_extensions");
    cmd.add(output_directory_arg);

    TCLAP::ValueArg<std::string> log_level_arg(
        "l", "log-level",
        "the verbosity of logging messages: none, error, warn, info, debug, "
        "all",
        false,
#ifdef NDEBUG
        "info",
#else
        "all",
#endif
        "LOG_LEVEL");
    cmd.add(log_level_arg);

    cmd.parse(argc, argv);

    BaseLib::setConsoleLogLevel(log_level_arg.getValue());
    spdlog::set_pattern("%^%l:%$ %v");
    spdlog::set_error_handler([](const std::string& msg) {
        std::cerr << "spdlog error: " << msg << std::endl;
        OGS_FATAL("spdlog logger error occurred.");
    });

    // init MPI
    MPI_Init(&argc, &argv);

    // start the timer
    BaseLib::RunTime run_timer;
    run_timer.start();
    BaseLib::CPUTime CPU_timer;
    CPU_timer.start();

    // add mpi_rank to logger output
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    spdlog::set_pattern(fmt::format("[{}] %^%l:%$ %v", mpi_rank));

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

    MPI_Finalize();

    return EXIT_SUCCESS;
}
