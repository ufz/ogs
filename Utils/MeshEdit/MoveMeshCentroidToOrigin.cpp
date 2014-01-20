/**
 * @file MoveMeshCentroidToOrigin.cpp
 * @date Jan 17, 2014
 * @brief 
 *
 * @copyright
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

// ThirdParty
#include "tclap/CmdLine.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "StringTools.h"
#include "FileTools.h"
#include "LogogSimpleFormatter.h"

// FileIO
#include "readMeshFromFile.h"
#include "FileIO/RapidXmlIO/BoostVtuInterface.h"

// GeoLib
#include "AABB.h"

// MeshLib
#include "Node.h"
#include "Elements/Element.h"
#include "Mesh.h"

int main(int argc, char *argv[])
{
	LOGOG_INITIALIZE();
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog::Cout *logogCout(new logog::Cout);
	logogCout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Moves the centroid of the mesh to the origin", ' ', "0.1");

	// Define a value argument and add it to the command line.
	// A value arg defines a flag and a type of value that it expects,
	// such as "-m meshfile".
	TCLAP::ValueArg<std::string> mesh_arg("m","mesh","input mesh file",true,"homer","string");

	// Add the argument mesh_arg to the CmdLine object. The CmdLine object
	// uses this Arg to parse the command line.
	cmd.add( mesh_arg );

	cmd.parse( argc, argv );

	std::string fname (mesh_arg.getValue());

	MeshLib::Mesh* mesh = FileIO::readMeshFromFile(fname);

	GeoLib::AABB<MeshLib::Node> aabb(mesh->getNodes().begin(), mesh->getNodes().end());
	MeshLib::Node center(
			(aabb.getMaxPoint()[0] + aabb.getMinPoint()[0])/2.0,
			(aabb.getMaxPoint()[1] + aabb.getMinPoint()[1])/2.0,
			0.0);
	INFO("translate model (-%f, -%f, -%f).", center[0], center[1], center[2]);
	std::for_each(mesh->getNodes().begin(), mesh->getNodes().end(),
					[&center](MeshLib::Node* node)
					{
							(*node)[0] -= center[0];
							(*node)[1] -= center[1];
					}
	);

	std::string out_fname(BaseLib::dropFileExtension(mesh_arg.getValue()));
	out_fname += "_transformed.vtu";
	FileIO::BoostVtuInterface mesh_io;
	mesh_io.setMesh(mesh);
	mesh_io.writeToFile(out_fname);

	delete mesh;
	delete logogCout;
	delete custom_format;
	LOGOG_SHUTDOWN();
}


