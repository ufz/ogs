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

#include "tclap/CmdLine.h"

#include "logog/include/logog.hpp"

#include "LogogSimpleFormatter.h"
#include "BaseLib/FileTools.h"

#include "FileIO/VtkIO/VtuInterface.h"
#include "Legacy/MeshIO.h"

#include "MeshPartitioning.h"

int main (int argc, char* argv[])
{
    LOGOG_INITIALIZE();
    logog::Cout* logog_cout (new logog::Cout);
    BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
    logog_cout->SetFormatter(*custom_format);

    std::string m_str = "Partition a mesh for parallel computing."
                        "The tasks of this tool are in twofold:\n"
                        "1. Convert mesh file to the input file of the partitioning tool,\n"
                        "2. Partition a mesh using the partitioning tool,\n"
                        "\tcreate the mesh data of each partition,\n"
                        "\trenumber the node indicies of each partition,\n"
                        "\tand output the results for parallel computing.";

    TCLAP::CmdLine cmd(m_str, ' ', "0.1");
    TCLAP::ValueArg<std::string> mesh_in("i", "mesh-input-file",
                                         "the name of the file containing the input mesh", true,
                                         "", "file name of input mesh");

    TCLAP::SwitchArg ogs2metis_flag("1", "ogs2metis",
                                    "Indicator to convert the ogs mesh file to METIS input file", cmd, false);

    cmd.parse(argc, argv);

    const std::string ifile_name = mesh_in.getValue();

    MeshLib::MeshPartitioning* mesh = dynamic_cast<MeshLib::MeshPartitioning*>
                                      (FileIO::VtuInterface::readVTUFile(ifile_name));
    INFO("Mesh read: %d nodes, %d elements.", mesh->getNNodes(), mesh->getNElements());

    if ( ogs2metis_flag.getValue() )
    {
        mesh->write2METIS(BaseLib::dropFileExtension(ifile_name) + ".mesh");
    }
    else
    {

    }

    delete custom_format;
    delete logog_cout;
    LOGOG_SHUTDOWN();

    return 0;
}
