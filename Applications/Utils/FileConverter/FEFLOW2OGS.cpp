/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// STL
#include <string>

// ThirdParty
#include "tclap/CmdLine.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "FileTools.h"
#include "RunTime.h"
#ifndef WIN32
#include "MemWatch.h"
#endif
#include "LogogSimpleFormatter.h"

// FileIO
#include "Legacy/MeshIO.h"
#include "FEFLOWInterface.h"
#include "XmlIO/Boost/BoostVtuInterface.h"

// MeshLib
#include "Mesh.h"

int main (int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logog_cout (new logog::Cout);
	BaseLib::LogogSimpleFormatter *const custom_format (new BaseLib::LogogSimpleFormatter);
	logog_cout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Converting a mesh in FEFLOW file format (ASCII, version 5.4) to a vtk unstructured grid file (new OGS file format) or to the old OGS file format - see options.", ' ', "0.1");

	TCLAP::ValueArg<std::string> ogs_mesh_arg(
		"o",
		"out",
		"filename for output mesh (if extension is msh, old OGS fileformat is written)",
		true,
		"",
		"filename as string");
	cmd.add(ogs_mesh_arg);

	TCLAP::ValueArg<std::string> feflow_mesh_arg(
		"i",
		"in",
		"FEFLOW input file (*.fem)",
		true,
		"",
		"filename as string");
	cmd.add(feflow_mesh_arg);

	cmd.parse(argc, argv);

	// *** read mesh
	INFO("Reading %s.", feflow_mesh_arg.getValue().c_str());
#ifndef WIN32
	BaseLib::MemWatch mem_watch;
	unsigned long mem_without_mesh (mem_watch.getVirtMemUsage());
#endif
	BaseLib::RunTime run_time;
	run_time.start();
	FileIO::FEFLOWInterface feflowIO(nullptr);
	MeshLib::Mesh const*const mesh(feflowIO.readFEFLOWFile(feflow_mesh_arg.getValue()));

	if (mesh == nullptr) {
		INFO("Could not read mesh from %s.", feflow_mesh_arg.getValue().c_str());
		return -1;
	}
#ifndef WIN32
	unsigned long mem_with_mesh (mem_watch.getVirtMemUsage());
	INFO("Mem for mesh: %i MB", (mem_with_mesh - mem_without_mesh)/(1024*1024));
#endif
	INFO("Time for reading: %f seconds.", run_time.elapsed());
	INFO("Read %d nodes and %d elements.", mesh->getNNodes(), mesh->getNElements());

	// *** write mesh in new format
	std::string ogs_mesh_fname(ogs_mesh_arg.getValue());

	if (BaseLib::getFileExtension(ogs_mesh_fname).compare("msh") == 0) {
		INFO("Writing %s.", ogs_mesh_fname.c_str());
		FileIO::Legacy::MeshIO mesh_io;
		mesh_io.setMesh(mesh);
		mesh_io.writeToFile(ogs_mesh_fname);
	} else {
		if (BaseLib::getFileExtension(ogs_mesh_fname).compare("vtu") != 0) {
			ogs_mesh_fname += ".vtu";
		}
		INFO("Writing %s.", ogs_mesh_fname.c_str());
		FileIO::BoostVtuInterface mesh_io;
		mesh_io.setMesh(mesh);
		mesh_io.writeToFile(ogs_mesh_fname);
	}
	INFO("\tDone.");

	delete mesh;
}

