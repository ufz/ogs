/**
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

// STL
#include <memory>
#include <string>
#include <vector>

// TCLAP
#include "tclap/CmdLine.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "BaseLib/BuildInfo.h"
#include "BaseLib/LogogSimpleFormatter.h"
#include "BaseLib/FileTools.h"

// GeoLib
#include "GeoLib/Point.h"
#include "GeoLib/Surface.h"

// FileIO
#include "FileIO/XmlIO/Boost/BoostVtuInterface.h"
#include "FileIO/TINInterface.h"

// MeshLib
#include "MeshLib/Mesh.h"
#include "MeshLib/convertMeshToGeo.h"


int main (int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logog_cout (new logog::Cout);
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog_cout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Converts TIN file into VTU file.", ' ', BaseLib::BuildInfo::git_version_sha1);
	TCLAP::ValueArg<std::string> inArg("i", "input-tin-file",
	                                     "the name of the file containing the input TIN", true,
	                                     "", "string");
	cmd.add(inArg);
	TCLAP::ValueArg<std::string> outArg("o", "output-vtu-file",
	                                      "the name of the file the mesh will be written to", true,
	                                      "", "string");
	cmd.add(outArg);
	cmd.parse(argc, argv);

	INFO("reading the TIN file...");
	const std::string tinFileName(inArg.getValue());
	std::vector<GeoLib::Point*> pnt_vec;
	std::unique_ptr<GeoLib::Surface> sfc(FileIO::TINInterface::readTIN(tinFileName, pnt_vec));
	if (!sfc)
		return 1;
	INFO("TIN read:  %d points, %d triangles", pnt_vec.size(), sfc->getNTriangles());

	INFO("converting to mesh data");
	std::unique_ptr<MeshLib::Mesh> mesh(MeshLib::convertSurfaceToMesh(*sfc, BaseLib::extractBaseNameWithoutExtension(tinFileName), std::numeric_limits<double>::epsilon()));
	INFO("Mesh created: %d nodes, %d elements.", mesh->getNNodes(), mesh->getNElements());

	INFO("Write it into VTU");
	FileIO::BoostVtuInterface writer;
	writer.setMesh(mesh.get());
	writer.writeToFile(outArg.getValue());

	for (auto p : pnt_vec) delete p;
	delete custom_format;
	delete logog_cout;
	LOGOG_SHUTDOWN();

	return 0;
}
