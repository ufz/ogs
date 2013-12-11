/**
 * \file DataExplorer5NodeReordering.cpp
 * 2013/13/06 KR Initial implementation
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <array>
#include <algorithm>
#include <vector>

#include "logog/include/logog.hpp"
#include "LogogSimpleFormatter.h"

// FileIO
#include "readMeshFromFile.h"
#include "RapidXmlIO/BoostVtuInterface.h"

// MeshLib
#include "Mesh.h"
#include "Elements/Element.h"

void reorderNodes(std::vector<MeshLib::Element*> &elements)
{
	std::size_t nElements (elements.size());
	for (std::size_t i=0; i<nElements; ++i) 
	{
		const unsigned nElemNodes (elements[i]->getNNodes());
		std::vector<MeshLib::Node*> nodes(elements[i]->getNodes(), elements[i]->getNodes() + nElemNodes);

		switch (elements[i]->getGeomType())
		{
			case MeshElemType::TETRAHEDRON:
				for(size_t j = 0; j < 4; ++j)
					elements[i]->setNode(j, nodes[(j+1)%4]);
				break;
			case MeshElemType::PYRAMID:
				for(size_t j = 0; j < 5; ++j)
					elements[i]->setNode(j, nodes[(j+1)%5]);
				break;
			case MeshElemType::PRISM:
				for(size_t j = 0; j < 3; ++j)
				{
					elements[i]->setNode(j, nodes[j+3]);
					elements[i]->setNode(j+3, nodes[j]);
				}
				break;
			case MeshElemType::HEXAHEDRON:
				for(size_t j = 0; j < 4; ++j)
				{
					elements[i]->setNode(j, nodes[j+4]);
					elements[i]->setNode(j+4, nodes[j]);
				}
				break;
			default:
				for(size_t j = 0; j < nElemNodes; ++j)
					elements[i]->setNode(j, nodes[nElemNodes - j - 1]);
		}
	}
}

int main (int argc, char* argv[])
{
	if (argc != 3)
	{
		std::cout << "Reordering of mesh nodes to make OGS Data Explorer 5 meshes compatible with OGS6." << std::endl;
		std::cout << std::endl;
		std::cout << "Usage: " << argv[0] << " <oldmesh.msh> <newmesh.vtu>" << std::endl;
		return 0;
	}
	
	LOGOG_INITIALIZE();
	logog::Cout* logogCout = new logog::Cout;
	BaseLib::LogogSimpleFormatter* formatter = new BaseLib::LogogSimpleFormatter;
	logogCout->SetFormatter(*formatter);

	MeshLib::Mesh* mesh (FileIO::readMeshFromFile(argv[1]));

	INFO("Reordering nodes... ");
	reorderNodes(const_cast<std::vector<MeshLib::Element*>&>(mesh->getElements()));
	
	FileIO::BoostVtuInterface writer;
	writer.setMesh(mesh);
	writer.writeToFile(argv[2]);

	INFO("VTU file written."); 

	delete formatter;
	delete logogCout;
	LOGOG_SHUTDOWN();

	return 1;
}



