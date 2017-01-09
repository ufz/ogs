/**
 * \file
 * \author Karsten Rink
 * \date   2012/05/09
 * \brief  Test for reading meshes.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <tclap/CmdLine.h>
#include <logog/include/logog.hpp>

#include "BaseLib/LogogSimpleFormatter.h"
#include "BaseLib/MemWatch.h"
#include "BaseLib/RunTime.h"
#include "BaseLib/StringTools.h"
#include "BaseLib/FileTools.h"

#include "MeshLib/IO/readMeshFromFile.h"

#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"


int main(int argc, char *argv[])
{
    LOGOG_INITIALIZE();
    BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
    logog::Cout *logogCout(new logog::Cout);
    logogCout->SetFormatter(*custom_format);

    TCLAP::CmdLine cmd("Simple mesh loading test", ' ', "0.1");

    // Define a value argument and add it to the command line.
    // A value arg defines a flag and a type of value that it expects,
    // such as "-m meshfile".
    TCLAP::ValueArg<std::string> mesh_arg("m","mesh","input mesh file",true,"homer","string");

    // Add the argument mesh_arg to the CmdLine object. The CmdLine object
    // uses this Arg to parse the command line.
    cmd.add( mesh_arg );

    cmd.parse( argc, argv );

    std::string fname (mesh_arg.getValue());

#ifndef WIN32
    BaseLib::MemWatch mem_watch;
    unsigned long mem_without_mesh (mem_watch.getVirtMemUsage());
#endif
    BaseLib::RunTime run_time;
    run_time.start();
    MeshLib::Mesh* mesh = MeshLib::IO::readMeshFromFile(fname);
#ifndef WIN32
    unsigned long mem_with_mesh (mem_watch.getVirtMemUsage());
//    std::cout << "mem for mesh: " << (mem_with_mesh - mem_without_mesh)/(1024*1024) << " MB" << std::endl;
    INFO ("mem for mesh: %i MB", (mem_with_mesh - mem_without_mesh)/(1024*1024));
#endif

//    std::cout << "time for reading: " << run_time.elapsed() << " s" << std::endl;
    INFO ("time for reading: %f s", run_time.elapsed());

/*
    unsigned elem_id = 25000;
    const MeshLib::Element* e = mesh->getElement(elem_id);
    const std::size_t nElems = mesh->getNumberOfElements();
    for (unsigned i=0; i< e->getNumberOfNeighbors(); i++)
    {
        for (unsigned j=0; j< nElems; j++)
            if (mesh->getElement(j) == e->getNeighbor(i))
                std::cout << "neighbour of element " << elem_id << " : " << j << std::endl;
    }
*/

    delete mesh;
    delete logogCout;
    delete custom_format;
    LOGOG_SHUTDOWN();
}

