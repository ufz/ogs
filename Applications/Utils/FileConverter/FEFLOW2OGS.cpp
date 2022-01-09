/**
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <memory>
#include <string>

// ThirdParty
#include <tclap/CmdLine.h>

// BaseLib
#include "BaseLib/FileTools.h"
#include "BaseLib/RunTime.h"
#include "InfoLib/GitInfo.h"
#ifndef WIN32
#include "BaseLib/MemWatch.h"
#endif

// FileIO
#include "Applications/FileIO/FEFLOW/FEFLOWMeshInterface.h"
#include "MeshLib/IO/Legacy/MeshIO.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Converting a mesh in FEFLOW file format (ASCII, version 5.4) to a vtk "
        "unstructured grid file (new OGS file format) or to the old OGS file "
        "format - see options.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2022, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::ValueArg<std::string> ogs_mesh_arg(
        "o",
        "out",
        "filename for output mesh (if extension is msh, old OGS fileformat is "
        "written)",
        true,
        "",
        "filename as string");
    cmd.add(ogs_mesh_arg);

    TCLAP::ValueArg<std::string> feflow_mesh_arg(
        "i", "in", "FEFLOW input file (*.fem)", true, "", "filename as string");
    cmd.add(feflow_mesh_arg);

    cmd.parse(argc, argv);

    // *** read mesh
    INFO("Reading {:s}.", feflow_mesh_arg.getValue());
#ifndef WIN32
    BaseLib::MemWatch mem_watch;
    unsigned long mem_without_mesh(mem_watch.getVirtMemUsage());
#endif
    BaseLib::RunTime run_time;
    run_time.start();
    FileIO::FEFLOWMeshInterface feflowIO;
    std::unique_ptr<MeshLib::Mesh const> mesh(
        feflowIO.readFEFLOWFile(feflow_mesh_arg.getValue()));

    if (mesh == nullptr)
    {
        INFO("Could not read mesh from {:s}.", feflow_mesh_arg.getValue());
        return EXIT_FAILURE;
    }
#ifndef WIN32
    unsigned long mem_with_mesh(mem_watch.getVirtMemUsage());
    INFO("Mem for mesh: {} MiB",
         (mem_with_mesh - mem_without_mesh) / (1024 * 1024));
#endif
    INFO("Time for reading: {:f} seconds.", run_time.elapsed());
    INFO("Read {:d} nodes and {:d} elements.", mesh->getNumberOfNodes(),
         mesh->getNumberOfElements());

    std::string ogs_mesh_fname(ogs_mesh_arg.getValue());
    INFO("Writing {:s}.", ogs_mesh_fname);
    MeshLib::IO::writeMeshToFile(*mesh, ogs_mesh_fname);
    INFO("\tDone.");
    return EXIT_SUCCESS;
}
