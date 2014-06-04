/**
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

// TCLAP
#include "tclap/CmdLine.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "LogogSimpleFormatter.h"
#include "FileTools.h"

// GeoLib
#include "GEOObjects.h"
#include "PolylineVec.h"

// FileIO
#include "Legacy/MeshIO.h"
#include "readMeshFromFile.h"
#include "XmlIO/Boost/BoostXmlGmlInterface.h"

// MeshLib
#include "Mesh.h"

// MeshGeoToolsLib
#include "MeshGeoToolsLib/AppendLinesAlongPolyline.h"


int main (int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logog_cout (new logog::Cout);
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog_cout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Append line elements into a mesh.", ' ', "0.1");
	TCLAP::ValueArg<std::string> mesh_in("i", "mesh-input-file",
	                                     "the name of the file containing the input mesh", true,
	                                     "", "file name of input mesh");
	cmd.add(mesh_in);
	TCLAP::ValueArg<std::string> mesh_out("o", "mesh-output-file",
	                                      "the name of the file the mesh will be written to", true,
	                                      "", "file name of output mesh");
	cmd.add(mesh_out);
	TCLAP::ValueArg<std::string> geoFileArg("g", "geo-file",
	                                      "the name of the geometry file which contains polylines", true, "", "the name of the geometry file");
	cmd.add(geoFileArg);

	// parse arguments
	cmd.parse(argc, argv);

	// read GEO objects
	GeoLib::GEOObjects geo_objs;
	FileIO::BoostXmlGmlInterface xml(geo_objs);
	xml.readFile(geoFileArg.getValue());

	std::vector<std::string> geo_names;
	geo_objs.getGeometryNames (geo_names);
	if (geo_names.empty ())
	{
		std::cout << "no geometries found" << std::endl;
		return -1;
	}
	const GeoLib::PolylineVec* ply_vec (geo_objs.getPolylineVecObj(geo_names[0]));
	if (!ply_vec)
	{
		std::cout << "could not found polylines" << std::endl;
		return -1;
	}

	// read a mesh
	MeshLib::Mesh const*const mesh (FileIO::readMeshFromFile(mesh_in.getValue()));
	if (!mesh)
	{
		ERR("Mesh file %s not found", mesh_in.getValue().c_str());
		return 1;
	}
	INFO("Mesh read: %d nodes, %d elements.", mesh->getNNodes(), mesh->getNElements());

	// add line elements
	MeshLib::Mesh* new_mesh = MeshGeoToolsLib::appendLinesAlongPolylines(*mesh, *ply_vec);
	INFO("Mesh created: %d nodes, %d elements.", new_mesh->getNNodes(), new_mesh->getNElements());

	// write into a file
	FileIO::Legacy::MeshIO meshIO;
	meshIO.setMesh(new_mesh);
	meshIO.writeToFile(mesh_out.getValue());

	delete custom_format;
	delete logog_cout;
	LOGOG_SHUTDOWN();

	return 1;
}

