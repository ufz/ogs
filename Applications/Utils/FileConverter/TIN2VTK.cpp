/**
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

// STL
#include <string>
#include <fstream>
#include <memory>

// TCLAP
#include "tclap/CmdLine.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "BaseLib/LogogSimpleFormatter.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"

// FileIO
#include "FileIO/XmlIO/Boost/BoostVtuInterface.h"
#include "FileIO/Legacy/MeshIO.h"

// MeshLib
#include "MeshLib/Mesh.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Tri.h"
#include "MeshLib/Node.h"
#include "MeshLib/MeshEditing/MeshRevision.h"

MeshLib::Mesh* readTIN(const std::string &tinFile)
{
	std::ifstream is(tinFile.c_str());
	if (!is) {
		ERR("Error: cannot open a file %s", tinFile.c_str());
		return nullptr;
	}

	std::vector<MeshLib::Node*> nodes;
	std::vector<MeshLib::Element*> elements;

	std::string line;
	std::size_t eId;
	double coords[3];
	std::size_t nodeId = 0;

	while (!is.eof())
	{
		getline(is, line);
		BaseLib::simplify(line);
		if (line.empty())
			continue;

		std::stringstream ss(line);
		ss >> eId;
		MeshLib::Node** tri_nodes = new MeshLib::Node*[3];
		for (unsigned i=0; i<3; i++) {
			for (unsigned j=0; j<3; j++)
				ss >> coords[j];
			tri_nodes[i] = new MeshLib::Node(coords, nodeId++);
		}
		elements.push_back(new MeshLib::Tri(tri_nodes, 0, eId));
		for (unsigned i=0; i<3; i++)
			nodes.push_back(tri_nodes[i]);
	}

	is.close();

	return new MeshLib::Mesh(BaseLib::extractBaseNameWithoutExtension(tinFile), nodes, elements);
}

int main (int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logog_cout (new logog::Cout);
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog_cout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Converts TIN file into VTK mesh.", ' ', "0.1");
	TCLAP::ValueArg<std::string> inArg("i", "input-file",
	                                     "the name of the file containing the input TIN", true,
	                                     "", "string");
	cmd.add(inArg);
	TCLAP::ValueArg<std::string> outArg("o", "output-file",
	                                      "the name of the file the mesh will be written to", true,
	                                      "", "string");
	cmd.add(outArg);
	cmd.parse(argc, argv);

	INFO("reading the TIN file...");
	std::unique_ptr<MeshLib::Mesh> mesh_with_duplicated_nodes(readTIN(inArg.getValue()));
	if (!mesh_with_duplicated_nodes)
		return 1;
	INFO("TIN read:  %d points, %d triangles", mesh_with_duplicated_nodes->getNNodes(), mesh_with_duplicated_nodes->getNElements());

	INFO("removing duplicated nodes");
	MeshLib::MeshRevision rev(*mesh_with_duplicated_nodes);
	std::unique_ptr<MeshLib::Mesh> mesh(rev.simplifyMesh(mesh_with_duplicated_nodes->getName(), std::numeric_limits<double>::epsilon()));
	INFO("Mesh created: %d nodes, %d elements.", mesh->getNNodes(), mesh->getNElements());

	INFO("Write it into VTK");
	FileIO::BoostVtuInterface writer;
	writer.setMesh(mesh.get());
	writer.writeToFile(outArg.getValue());

	delete custom_format;
	delete logog_cout;
	LOGOG_SHUTDOWN();

	return 0;
}
