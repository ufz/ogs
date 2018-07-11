/*!
  \file PartitionMesh.cpp
  \date   2016.05

  \brief  A tool for mesh partitioning.

  \copyright
  Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/

#include <tclap/CmdLine.h>

#ifdef WIN32
#include <windows.h>
#else
#include <climits>
#endif

#include "Applications/ApplicationsLib/LogogSetup.h"
#include "BaseLib/CPUTime.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/RunTime.h"

#include "MeshLib/IO/readMeshFromFile.h"

#include "NodeWiseMeshPartitioner.h"
#include "Metis.h"

using namespace ApplicationUtils;

int main(int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    const std::string m_str =
        "Partition a mesh for parallel computing."
        "The tasks of this tool are in twofold:\n"
        "1. Convert mesh file to the input file of the partitioning tool,\n"
        "2. Partition a mesh using the partitioning tool,\n"
        "\tcreate the mesh data of each partition,\n"
        "\trenumber the node indices of each partition,\n"
        "\tand output the results for parallel computing.\n"
        "Note: If this tool is installed as a system command,\n"
        "\tthe command must be run with its full path.";

    TCLAP::CmdLine cmd(m_str, ' ', "0.1");
    TCLAP::ValueArg<std::string> mesh_input(
        "i", "mesh-input-file",
        "the name of the file containing the input mesh", true, "",
        "file name of input mesh");
    cmd.add(mesh_input);

    TCLAP::ValueArg<std::string> output_directory_arg(
        "o", "output", "directory name for the output files", false, "",
        "directory");
    cmd.add(output_directory_arg);

    TCLAP::ValueArg<int> nparts("n", "np", "the number of partitions", false, 2,
                                "integer");
    cmd.add(nparts);

    TCLAP::SwitchArg ogs2metis_flag(
        "s", "ogs2metis",
        "Indicator to convert the ogs mesh file to METIS input file", cmd,
        false);

    TCLAP::SwitchArg exe_metis_flag(
        "m", "exe_metis", "Call mpmetis inside the programme via system().",
        false);
    cmd.add(exe_metis_flag);

    TCLAP::SwitchArg lh_elems_flag(
        "q", "lh_elements", "Mixed linear and high order elements.", false);
    cmd.add(lh_elems_flag);

    TCLAP::SwitchArg ascii_flag("a", "ascii", "Enable ASCII output.", false);
    cmd.add(ascii_flag);

    cmd.parse(argc, argv);

    BaseLib::RunTime run_timer;
    run_timer.start();
    BaseLib::CPUTime CPU_timer;
    CPU_timer.start();

    const std::string input_file_name_wo_extension =
        BaseLib::dropFileExtension(mesh_input.getValue());
    std::unique_ptr<MeshLib::Mesh> mesh_ptr(
        MeshLib::IO::readMeshFromFile(input_file_name_wo_extension + ".vtu"));
    INFO("Mesh read: %d nodes, %d elements.",
         mesh_ptr->getNumberOfNodes(),
         mesh_ptr->getNumberOfElements());

    if (ogs2metis_flag.getValue())
    {
        INFO("Write the mesh into METIS input file.");
        ApplicationUtils::writeMETIS(mesh_ptr->getElements(),
                                     input_file_name_wo_extension + ".mesh");
        INFO("Total runtime: %g s.", run_timer.elapsed());
        INFO("Total CPU time: %g s.", CPU_timer.elapsed());

        return EXIT_SUCCESS;
    }

    ApplicationUtils::NodeWiseMeshPartitioner mesh_partitioner(
        nparts.getValue(), std::move(mesh_ptr));

    std::string const output_file_name_wo_extension = BaseLib::joinPaths(
        output_directory_arg.getValue(),
        BaseLib::extractBaseNameWithoutExtension(mesh_input.getValue()));
    const int num_partitions = nparts.getValue();

    if (num_partitions < 1)
    {
        OGS_FATAL("Number of partitions must be positive.");
    }

    if (num_partitions == 1)
    {
        OGS_FATAL(
            "Partitioning the mesh into one domain is unnecessary because OGS "
            "reads vtu mesh data directly when called with 'mpirun bin/ogs "
            "-np=1'.");
    }

    // Execute mpmetis via system(...)
    if (exe_metis_flag.getValue())
    {
        INFO("METIS is running ...");
        const std::string exe_name = argv[0];
        const std::string exe_path = BaseLib::extractPath(exe_name);
        INFO("Path to mpmetis is: \n\t%s", exe_path.c_str());

        const std::string mpmetis_com =
            exe_path + "/mpmetis " + " -gtype=nodal " +
            input_file_name_wo_extension + ".mesh " +
            std::to_string(nparts.getValue());

        const int status = system(mpmetis_com.c_str());
        if (status != 0)
        {
            INFO("Failed in system calling.");
            INFO("Return value of system call %d ", status);
            return EXIT_FAILURE;
        }
    }
    mesh_partitioner.resetPartitionIdsForNodes(
        readMetisData(input_file_name_wo_extension, num_partitions,
                      mesh_partitioner.mesh().getNumberOfNodes()));

    removeMetisPartitioningFiles(input_file_name_wo_extension, num_partitions);

    INFO("Partitioning the mesh in the node wise way ...");
    mesh_partitioner.partitionByMETIS(lh_elems_flag.getValue());
    if (ascii_flag.getValue())
    {
        INFO("Write the data of partitions into ASCII files ...");
        mesh_partitioner.writeASCII(output_file_name_wo_extension);
    }
    else
    {
        INFO("Write the data of partitions into binary files ...");
        mesh_partitioner.writeBinary(output_file_name_wo_extension);
    }

    INFO("Total runtime: %g s.", run_timer.elapsed());
    INFO("Total CPU time: %g s.", CPU_timer.elapsed());

    return EXIT_SUCCESS;
}
