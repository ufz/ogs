/*
 * MeshRead.cpp
 *
 *  Created on: 2012/05/09
 *      Author: KR
 */

// BaseLib
#include "MemWatch.h"
#include "RunTime.h"
#include "tclap/CmdLine.h"

// BaseLib/logog
#include "logog.hpp"


// MeshLib
#include "Node.h"
#include "Elements/Element.h"
#include "Mesh.h"
#include "MeshIO.h"

int main(int argc, char *argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logogCout = new logog::Cout;

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

	FileIO::MeshIO mesh_io;
#ifndef WIN32
	BaseLib::MemWatch mem_watch;
	unsigned long mem_without_mesh (mem_watch.getVirtMemUsage());
	BaseLib::RunTime run_time;
	run_time.start();
#endif
	MeshLib::Mesh* mesh = mesh_io.loadMeshFromFile(fname);
#ifndef WIN32
	unsigned long mem_with_mesh (mem_watch.getVirtMemUsage());
//	std::cout << "mem for mesh: " << (mem_with_mesh - mem_without_mesh)/(1024*1024) << " MB" << std::endl;
	INFO ("mem for mesh: %i MB", (mem_with_mesh - mem_without_mesh)/(1024*1024));
	run_time.stop();
//	std::cout << "time for reading: " << run_time.elapsed() << " s" << std::endl;
	INFO ("time for reading: %f s", run_time.elapsed());
#endif
	delete mesh;
	delete logogCout;
	LOGOG_SHUTDOWN();
}

