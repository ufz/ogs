/**
 * @file   reverseMeshNodeOrdering.cpp
 * @author Norihiro Watanabe
 * @date   2013/10/15
 * @brief  Reverse element node ordering
 *
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
#include "readMeshFromFile.h"

// MeshLib
#include "Mesh.h"
#include "Elements/Element.h"

void reverseNodeOrdering(std::vector<MeshLib::Element*> & ele_vec)
{
	for (MeshLib::Element* ele : ele_vec) {
		unsigned nElemNodes (ele->getNNodes());
		std::vector<MeshLib::Node*> originalNodes(ele->getNodes(), ele->getNodes() + nElemNodes);
		if (ele->getGeomType() == MeshElemType::PRISM)
		{
			for(size_t j = 0; j < 3; ++j)
			{
				ele->setNode(j, originalNodes[j+3]);
				ele->setNode(j+3, originalNodes[j]);
			}
		}
		else if (ele->getGeomType() == MeshElemType::TETRAHEDRON)
		{
			for(size_t j = 0; j < 4; ++j)
				ele->setNode(j, originalNodes[(j+1)%4]);
		}
		else 
		{
			for(size_t j = 0; j < nElemNodes; ++j)
				ele->setNode(j, originalNodes[nElemNodes - j - 1]);
		}
	}
}

int main (int argc, char* argv[])
{
	LOGOG_INITIALIZE();
	logog::Cout* logog_cout (new logog::Cout);
	BaseLib::LogogSimpleFormatter *custom_format (new BaseLib::LogogSimpleFormatter);
	logog_cout->SetFormatter(*custom_format);

	TCLAP::CmdLine cmd("Reverse the node ordering of mesh elements.", ' ', "0.1");
	TCLAP::ValueArg<std::string> mesh_in("i", "mesh-input-file",
	                                     "the name of the file containing the input mesh", true,
	                                     "", "file name of input mesh");
	cmd.add(mesh_in);
	TCLAP::ValueArg<std::string> mesh_out("o", "mesh-output-file",
	                                      "the name of the file the mesh will be written to", true,
	                                      "", "file name of output mesh");
	cmd.add(mesh_out);
	cmd.parse(argc, argv);

	MeshLib::Mesh* mesh (FileIO::readMeshFromFile(mesh_in.getValue()));
	INFO("Mesh read: %d nodes, %d elements.", mesh->getNNodes(), mesh->getNElements());

	INFO("Reversing the node ordering...");
	reverseNodeOrdering(const_cast<std::vector<MeshLib::Element*>&>(mesh->getElements()));

	FileIO::MeshIO meshIO;
	meshIO.setMesh(mesh);
	meshIO.writeToFile(mesh_out.getValue());

	delete custom_format;
	delete logog_cout;
	LOGOG_SHUTDOWN();

	return 0;
}



