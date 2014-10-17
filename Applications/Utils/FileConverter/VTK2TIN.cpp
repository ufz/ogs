/**
 * @copyright
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

// STL
#include <memory>
#include <string>
#include <fstream>

// TCLAP
#include "tclap/CmdLine.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "BaseLib/BuildInfo.h"
#include "BaseLib/LogogSimpleFormatter.h"

// GeoLib
#include "GeoLib/GEOObjects.h"
#include "GeoLib/Surface.h"

// MeshLib
#include "MeshLib/Mesh.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"
#include "MeshLib/convertMeshToGeo.h"

// FileIO
#include "FileIO/XmlIO/Boost/BoostVtuInterface.h"
#include "FileIO/TINInterface.h"


int main (int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logog_cout (new logog::Cout);
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog_cout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Converts VTK mesh into TIN file.", ' ', BaseLib::BuildInfo::git_version_sha1);
	TCLAP::ValueArg<std::string> mesh_in("i", "mesh-input-file",
	                                     "the name of the file containing the input mesh", true,
	                                     "", "file name of input mesh");
	cmd.add(mesh_in);
	TCLAP::ValueArg<std::string> mesh_out("o", "TIN-output-file",
	                                      "the name of the file the TIN will be written to", true,
	                                      "", "file name of output TIN");
	cmd.add(mesh_out);
	cmd.parse(argc, argv);

	std::unique_ptr<MeshLib::Mesh> mesh (FileIO::BoostVtuInterface::readVTUFile(mesh_in.getValue()));
	INFO("Mesh read: %d nodes, %d elements.", mesh->getNNodes(), mesh->getNElements());

	INFO("Converting the mesh to TIN");
	GeoLib::GEOObjects geo_objects;
	if (MeshLib::convertMeshToGeo(*mesh, geo_objects)) {
		INFO("Writing TIN into the file");
		FileIO::TINInterface::writeSurfaceAsTIN(*(*geo_objects.getSurfaceVec(mesh->getName()))[0], mesh_out.getValue());
	}

	delete custom_format;
	delete logog_cout;
	LOGOG_SHUTDOWN();

	return 0;
}
