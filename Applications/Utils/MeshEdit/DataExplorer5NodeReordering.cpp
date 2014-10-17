/**
 * \file DataExplorer5NodeReordering.cpp
 * 2013/13/06 KR Initial implementation
 *
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <array>
#include <algorithm>
#include <vector>

#include "logog/include/logog.hpp"
#include "LogogSimpleFormatter.h"

// TCLAP
#include "tclap/CmdLine.h"

// FileIO
#include "readMeshFromFile.h"
#include "FileIO/VtkIO/VtuInterface.h"

// MeshLib
#include "Mesh.h"
#include "Elements/Element.h"

void reorderNodes(std::vector<MeshLib::Element*> &elements)
{
	std::size_t nElements (elements.size());
	for (std::size_t i=0; i<nElements; ++i)
	{
		const unsigned nElemNodes (elements[i]->getNBaseNodes());
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
	LOGOG_INITIALIZE();
	logog::Cout* logogCout = new logog::Cout;
	BaseLib::LogogSimpleFormatter* formatter = new BaseLib::LogogSimpleFormatter;
	logogCout->SetFormatter(*formatter);

	TCLAP::CmdLine cmd("Reordering of mesh nodes to make OGS Data Explorer 5 meshes compatible with OGS6.",
			' ', "0.1");
	TCLAP::UnlabeledValueArg<std::string> input_mesh_arg("OGS5_file_input_mesh",
	                                           "the name of the mesh file used for input containing elements in OGS5 node ordering",
	                                           true,
	                                           "",
	                                           "oldmesh.msh");
	cmd.add(input_mesh_arg);
	TCLAP::UnlabeledValueArg<std::string> output_mesh_arg("OGS6_file_output_mesh",
	                                           "the name of the mesh file used for output with node ordering consistent to OGS-6",
	                                           true,
	                                           "",
	                                           "newmesh.vtu");
	cmd.add(output_mesh_arg);
	cmd.parse(argc, argv);

	MeshLib::Mesh* mesh (FileIO::readMeshFromFile(input_mesh_arg.getValue().c_str()));

	INFO("Reordering nodes... ");
	reorderNodes(const_cast<std::vector<MeshLib::Element*>&>(mesh->getElements()));

	FileIO::VtuInterface writer;
	writer.setMesh(mesh);
	writer.writeToFile(output_mesh_arg.getValue().c_str());

	INFO("VTU file written.");

	delete formatter;
	delete logogCout;
	LOGOG_SHUTDOWN();

	return 0;
}



