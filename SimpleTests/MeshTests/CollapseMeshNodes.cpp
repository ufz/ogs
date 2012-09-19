/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file CollapseMeshNodes.cpp
 *
 *  Created on  Jul 31, 2012 by Thomas Fischer
 */

// BaseLib
#include "MemWatch.h"
#include "RunTime.h"
#include "tclap/CmdLine.h"
#include "LogogSimpleFormatter.h"

// BaseLib/logog
#include "logog.hpp"

// MeshLib
#include "Node.h"
#include "Elements/Element.h"
#include "Mesh.h"
#include "Legacy/MeshIO.h"
#include "MeshCoarsener.h"

int main(int argc, char *argv[])
{
	LOGOG_INITIALIZE();
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog::Cout *logogCout(new logog::Cout);
	logogCout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Collapse mesh nodes and, if necessary, remove elements", ' ', "0.1");

	// Define a value argument and add it to the command line.
	// A value arg defines a flag and a type of value that it expects,
	// such as "-m meshfile".
	TCLAP::ValueArg<std::string> input_mesh_arg("m","mesh","input mesh file name",true,"","string");
	// Add the argument mesh_arg to the CmdLine object. The CmdLine object
	// uses this Arg to parse the command line.
	cmd.add( input_mesh_arg );

	TCLAP::ValueArg<std::string> output_mesh_arg("","out-mesh","mesh file name for output",false,"","string");
	cmd.add( output_mesh_arg );

	TCLAP::ValueArg<double> distance_arg("d","collapse-distance","maximal distance two nodes are collapsed",false,0.01,"for example you can set this parameter to 10^{-6} times maximal area length");
	cmd.add( distance_arg );

	cmd.parse( argc, argv );

	std::string fname (input_mesh_arg.getValue());

	FileIO::MeshIO mesh_io;
#ifndef WIN32
	BaseLib::MemWatch mem_watch;
	unsigned long mem_without_mesh (mem_watch.getVirtMemUsage());
	BaseLib::RunTime run_time;
	run_time.start();
#endif
	MeshLib::Mesh* mesh = mesh_io.loadMeshFromFile(fname);
#ifndef WIN32
	if (mesh) {
		unsigned long mem_with_mesh (mem_watch.getVirtMemUsage());
		INFO ("mem for mesh: %i MB", (mem_with_mesh - mem_without_mesh)/(1024*1024));
	}
	run_time.stop();
	if (mesh) {
		INFO ("time for reading: %f s", run_time.elapsed());
	}
#endif

#ifndef WIN32
	unsigned long mem_without_meshgrid (mem_watch.getVirtMemUsage());
	run_time.start();
#endif
	MeshLib::MeshCoarsener mesh_coarsener(mesh);
	MeshLib::Mesh *collapsed_mesh(mesh_coarsener (distance_arg.getValue()));

#ifndef WIN32
	run_time.stop();
	unsigned long mem_with_meshgrid (mem_watch.getVirtMemUsage());
	INFO ("mem for meshgrid: %i MB", (mem_with_meshgrid - mem_without_meshgrid)/(1024*1024));
	INFO ("time for collapsing: %f s", run_time.elapsed());
#endif

	mesh_io.setMesh(collapsed_mesh);
	std::string out_fname (output_mesh_arg.getValue());
	if (out_fname.empty()) {
		out_fname = "/home/fischeth/workspace/OGS-6/Build/CollapsedMesh.msh";
	}
	INFO ("writing collapsed mesh to %s", out_fname.c_str());
	mesh_io.writeToFile(out_fname);
	INFO ("done");

	delete mesh;
	delete collapsed_mesh;
	delete custom_format;
	delete logogCout;
	LOGOG_SHUTDOWN();

	return 0;
}
