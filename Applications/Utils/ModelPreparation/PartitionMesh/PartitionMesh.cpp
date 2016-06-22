/*!
  \file PartitionMesh.cpp
  \date   2016.05

  \brief  A tool for mesh partitioning.

  \copyright
  Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/

#include <sstream>
#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"
#include "BaseLib/FileTools.h"

#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"

#include "BaseLib/CPUTime.h"
#include "BaseLib/RunTime.h"

#include "MeshPartitioning.h"

int main (int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    std::string m_str = "Partition a mesh for parallel computing."
                        "The tasks of this tool are in twofold:\n"
                        "1. Convert mesh file to the input file of the partitioning tool,\n"
                        "2. Partition a mesh using the partitioning tool,\n"
                        "\tcreate the mesh data of each partition,\n"
                        "\trenumber the node indicies of each partition,\n"
                        "\tand output the results for parallel computing.";

    TCLAP::CmdLine cmd(m_str, ' ', "0.1");
    TCLAP::ValueArg<std::string> mesh_input("i", "mesh-input-file",
                                            "the name of the file containing the input mesh", true,
                                            "", "file name of input mesh");
    cmd.add(mesh_input);
    TCLAP::ValueArg<int> nparts("n", "np", "the number of partitions", false, 2, "integer");
    cmd.add(nparts);

    TCLAP::SwitchArg ogs2metis_flag("1", "ogs2metis",
                                    "Indicator to convert the ogs mesh file to METIS input file", cmd, false);

    TCLAP::SwitchArg exe_metis_flag("m","exe_metis","Call mpmetis inside the programme via system().", false);
    cmd.add(exe_metis_flag);
    TCLAP::SwitchArg asci_flag("a","asci","Enable ASCII output.", false);
    cmd.add(asci_flag);
    TCLAP::SwitchArg elem_part_flag("e","element_wise","Enable element wise partitioing.", false);
    cmd.add(elem_part_flag);

    cmd.parse(argc, argv);

    BaseLib::RunTime run_timer;
    run_timer.start();
    BaseLib::CPUTime CPU_timer;
    CPU_timer.start();

    const std::string ifile_name = mesh_input.getValue();
    const std::string file_name_base = BaseLib::dropFileExtension(ifile_name);
    MeshLib::MeshPartitioning* mesh = static_cast<MeshLib::MeshPartitioning*>
                                      (MeshLib::IO::readMeshFromFile(file_name_base + ".vtu"));
    INFO("Mesh read: %d nodes, %d elements.", mesh->getNumberOfNodes(), mesh->getNumberOfElements());

    if (ogs2metis_flag.getValue())
    {
        INFO("Write the mesh into METIS input file.");
        mesh->write2METIS(BaseLib::dropFileExtension(ifile_name) + ".mesh");
    }
    else
    {
        if ( elem_part_flag.getValue() )
        {
            INFO("Partitioning the mesh in the element wise way ...");
            mesh->partitionByElementMETIS( file_name_base, nparts.getValue());
        }
        else
        {
            const int num_partitions = nparts.getValue();
            std::string str_nparts = std::to_string(num_partitions);

            // Execute mpmetis via system(...)
            if (num_partitions > 1 && exe_metis_flag.getValue())
            {
                INFO("METIS is running ...");
                std::stringstream ss;
                ss << "mpmetis " << " -gtype=nodal "
                   << file_name_base + ".mesh "
                   << nparts.getValue();

                int status = system(ss.str().c_str());
                INFO("Return value of system calling %d ", status);
            }

            INFO("Partitioning the mesh in the node wise way ...");
            // Binary output is default.
            mesh->partitionByNodeMETIS(file_name_base, num_partitions,
                                       !asci_flag.getValue());

            INFO("Write the mesh with renumbered node indicies into VTU");
            MeshLib::IO::VtuInterface writer(mesh);
            writer.writeToFile(file_name_base + "_node_id_renumbered_partitions_"
                               + str_nparts + ".vtu");
        }
    }

    INFO( "Total runtime: %g s.\n", run_timer.elapsed() );
    INFO( "Total CPU time: %g s.\n", CPU_timer.elapsed() );

    delete mesh;

    return EXIT_SUCCESS;
}
