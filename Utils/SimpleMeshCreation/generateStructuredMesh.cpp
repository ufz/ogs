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

// FileIO
#include "Legacy/MeshIO.h"

// MeshLib
#include "Mesh.h"
#include "Node.h"
#include "Elements/Element.h"
#include "MeshEnums.h"
#include "MeshGenerators/MeshGenerator.h"


int main (int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logog_cout (new logog::Cout);
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog_cout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Generate a structured mesh.", ' ', "0.1");
	TCLAP::ValueArg<std::string> mesh_out("o", "mesh-output-file",
	                                      "the name of the file the mesh will be written to", true,
	                                      "", "file name of output mesh");
	cmd.add(mesh_out);
	TCLAP::ValueArg<std::string> eleTypeArg("e", "element-type to be created",
	                                      "element type to be removed", true, "line", "element type");
	cmd.add(eleTypeArg);
	TCLAP::ValueArg<double> lengthArg("l", "length",
	                                      "length of a domain", true, 10.0, "length of a domain");
	cmd.add(lengthArg);
	TCLAP::ValueArg<unsigned> nsubdivArg("n", "nr-subdivision",
	                                      "the number of subdivision", true, 10, "the number of subdivision");
	cmd.add(nsubdivArg);

	// parse arguments
	cmd.parse(argc, argv);
	const std::string eleTypeName(eleTypeArg.getValue());
	const MeshElemType eleType = String2MeshElemType(eleTypeName);
	const double length = lengthArg.getValue();
	const unsigned n_subdivision = nsubdivArg.getValue();

	// generate a mesh
	MeshLib::Mesh* mesh = nullptr;
	switch (eleType)
	{
	case MeshElemType::LINE:
		mesh = MeshLib::MeshGenerator::generateLineMesh(length, n_subdivision);
		break;
	case MeshElemType::TRIANGLE:
		mesh = MeshLib::MeshGenerator::generateRegularTriMesh(length, n_subdivision);
		break;
	case MeshElemType::QUAD:
		mesh = MeshLib::MeshGenerator::generateRegularQuadMesh(length, n_subdivision);
		break;
	case MeshElemType::HEXAHEDRON:
		mesh = MeshLib::MeshGenerator::generateRegularHexMesh(length, n_subdivision);
		break;
	default:
		ERR("Given element type is not supported.");
		break;
	}

	if (mesh)
	{
		INFO("Mesh created: %d nodes, %d elements.", mesh->getNNodes(), mesh->getNElements());

		// write into a file
		FileIO::Legacy::MeshIO meshIO;
		meshIO.setMesh(mesh);
		meshIO.writeToFile(mesh_out.getValue());

		delete mesh;
	}

	delete custom_format;
	delete logog_cout;
	LOGOG_SHUTDOWN();

	return 1;
}

