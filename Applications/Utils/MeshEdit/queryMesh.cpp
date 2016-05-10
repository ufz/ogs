/**
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 */

#include <array>
#include <memory>
#include <string>

#include "logog/include/logog.hpp"
#include "tclap/CmdLine.h"

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "BaseLib/BuildInfo.h"
#include "BaseLib/StringTools.h"
#include "BaseLib/FileTools.h"

#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"

#include "MeshLib/IO/readMeshFromFile.h"

int main(int argc, char *argv[])
{
	ApplicationsLib::LogogSetup logog_setup;

	TCLAP::CmdLine cmd("Query mesh information", ' ', BaseLib::BuildInfo::git_describe);
	TCLAP::UnlabeledValueArg<std::string> mesh_arg("mesh-file","input mesh file",true,"","string");
	cmd.add( mesh_arg );
	TCLAP::MultiArg<std::size_t> eleId_arg("e","element-id","element ID",false,"number");
	cmd.add( eleId_arg );
	TCLAP::MultiArg<std::size_t> nodeId_arg("n","node-id","node ID",false,"number");
	cmd.add( nodeId_arg );

	cmd.parse( argc, argv );

	const std::string filename(mesh_arg.getValue());

	// read the mesh file
	auto const mesh = std::unique_ptr<MeshLib::Mesh>(
		MeshLib::IO::readMeshFromFile(filename));
	if (!mesh)
		return EXIT_FAILURE;

	auto materialIds = mesh->getProperties().getPropertyVector<int>("MaterialIDs");

	std::cout << std::scientific << std::setprecision(12);
	for (auto ele_id : eleId_arg.getValue())
	{
		std::cout << "--------------------------------------------------------" << std::endl;
		auto* ele = mesh->getElement(ele_id);
		std::cout << "# Element " << ele->getID() << std::endl;
		std::cout << "Type : " << CellType2String(ele->getCellType()) << std::endl;
		if (materialIds)
			std::cout << "Mat ID : " << (*materialIds)[ele_id] << std::endl;
		std::cout << "Nodes: " << std::endl;
		for (unsigned i=0; i<ele->getNNodes(); i++)
			std::cout <<  ele->getNode(i)->getID() << " " << *ele->getNode(i) << std::endl;
		std::cout << "Content: " << ele->getContent() << std::endl;
		std::cout << "Neighbors: ";
		for (unsigned i=0; i<ele->getNNeighbors(); i++)
		{
			if (ele->getNeighbor(i))
				std::cout << ele->getNeighbor(i)->getID() << " ";
			else
				std::cout << "none ";
		}
		std::cout << std::endl;
	}

	for (auto node_id : nodeId_arg.getValue())
	{
		std::cout << "--------------------------------------------------------" << std::endl;
		auto* node = mesh->getNode(node_id);
		std::cout << "# Node" << node->getID() << std::endl;
		std::cout << "Coordinates: " << *node << std::endl;
		std::cout << "Connected elements: " ;
		for (unsigned i=0; i<node->getNElements(); i++)
			std::cout << node->getElement(i)->getID() << " ";
		std::cout << std::endl;
	}
}
