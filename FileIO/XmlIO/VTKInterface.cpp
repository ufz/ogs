/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file VTKInterface.cpp
 *
 *  Created on 2012-08-30 by Karsten Rink
 */

#include "VTKInterface.h"
#include <iostream>
#include <fstream>

#include "StringTools.h"

// MSH
#include "Mesh.h"
#include "Node.h"
#include "Elements/Edge.h"
#include "Elements/Tri.h"
#include "Elements/Quad.h"
#include "Elements/Tet.h"
#include "Elements/Hex.h"
#include "Elements/Pyramid.h"
#include "Elements/Prism.h"


#include "RapidXML/rapidxml.hpp"

namespace FileIO {

VTKInterface::VTKInterface()
{
}

VTKInterface::~VTKInterface()
{
}

MeshLib::Mesh* VTKInterface::readVTUFile(const std::string &file_name)
{
	std::ifstream in(file_name.c_str());
	if (in.fail())
	{
		std::cout << "VTKInterface::readVTUFile() - Can't open xml-file." << std::endl;
		return NULL;
	}

	in.seekg(0, std::ios::end);
	size_t length = in.tellg();
	in.seekg(0, std::ios::beg);
	char* buffer = new char[length+1];
	in.read(buffer, length);
	buffer[in.gcount()] = '\0';
	in.close();

	// build DOM tree
	rapidxml::xml_document<> doc;
	doc.parse<0>(buffer);

	// error management
	const rapidxml::xml_node<>* vtk_root (doc.first_node());
	if (vtk_root == NULL || std::string(vtk_root->name()).compare("VTKFile"))
	{
		std::cout << "Error in VTKInterface::readVTUFile() - Not a VTK File." << std::endl;
		return NULL;
	}

	if (std::string(vtk_root->first_attribute("version")->value()).compare("0.1"))
	{
		std::cout << "Error in VTKInterface::readVTUFile() - Unsupported file format version." << std::endl;
		return NULL;
	}

	if (std::string(vtk_root->first_attribute("byte_order")->value()).compare("LittleEndian"))
	{
		std::cout << "Error in VTKInterface::readVTUFile() - Only little endian files are supported." << std::endl;
		return NULL;
	}
	
	const rapidxml::xml_attribute<>* compressor (vtk_root->first_attribute("compressor"));
	const rapidxml::xml_node<>* grid_root (vtk_root->first_node());
	if (std::string(grid_root->name()).compare("UnstructuredGrid"))
	{
		std::cout << "Error in VTKInterface::readVTUFile() - Not an unstructured grid." << std::endl;
		return NULL;
	}
	
	//parse content
	const rapidxml::xml_node<>* piece_node (grid_root->first_node("Piece"));
	if (piece_node)
	{
		const unsigned nNodes = static_cast<unsigned>(atoi(piece_node->first_attribute("NumberOfPoints")->value()));
		const unsigned nElems = static_cast<unsigned>(atoi(piece_node->first_attribute("NumberOfCells")->value()));
		std::vector<MeshLib::Node*> nodes(nNodes);
		std::vector<MeshLib::Element*> elements(nElems);
		std::vector<unsigned> mat_ids(nElems, 0);
		std::vector<unsigned> cell_types(nElems);

		const rapidxml::xml_node<>* mat_id_node (piece_node->first_node("CellData")->first_node("DataArray"));
		if (mat_id_node && (std::string(mat_id_node->first_attribute("Name")->value()).compare("MaterialIDs") == 0))
		{
			if (std::string(mat_id_node->first_attribute("format")->value()).compare("ascii") == 0)
			{
				std::stringstream iss (mat_id_node->value());
				for(unsigned i=0; i<nElems; i++)
					iss >> mat_ids[i];
			}
		}
		else
			std::cout << "Warning in VTKInterface::readVTUFile() - MaterialID array not found." << std::endl;


		const rapidxml::xml_node<>* points_node (piece_node->first_node("Points")->first_node("DataArray"));
		if (points_node && (std::string(points_node->first_attribute("Name")->value()).compare("Points") == 0))
		{
			if (std::string(points_node->first_attribute("format")->value()).compare("ascii") == 0)
			{
				std::stringstream iss (points_node->value());
				double x,y,z;
				for(unsigned i=0; i<nNodes; i++)
				{
					iss >> x >> y >> z;
					nodes[i] = new MeshLib::Node(x,y,z,i);
				}
			}
		}
		else
		{
			std::cout << "Error in VTKInterface::readVTUFile() - Points array not found." << std::endl;
			return NULL;
		}

		rapidxml::xml_node<>* cells_node (piece_node->first_node("Cells"));
		rapidxml::xml_node<>* connectivity_node (NULL);
		rapidxml::xml_node<>* celltype_node (NULL);
		for (rapidxml::xml_node<>* cells_data = cells_node->first_node("DataArray"); cells_data; cells_data = cells_data->next_sibling("DataArray"))
		{
			if (std::string(cells_data->first_attribute("Name")->value()).compare("connectivity") == 0)
				connectivity_node = cells_data;
			if (std::string(cells_data->first_attribute("Name")->value()).compare("types") == 0)
				celltype_node = cells_data;
		}

		if (connectivity_node && celltype_node)
		{
			if (std::string(celltype_node->first_attribute("format")->value()).compare("ascii") == 0)
			{
				std::stringstream iss (celltype_node->value());
				for(unsigned i=0; i<nElems; i++)
					iss >> cell_types[i];
			}
			if (std::string(connectivity_node->first_attribute("format")->value()).compare("ascii") == 0)
			{
				std::stringstream iss (connectivity_node->value());
				for(unsigned i=0; i<nElems; i++)
					elements[i] = readElement(iss, nodes, mat_ids[i], cell_types[i]);
			}
		}
		else
		{
			std::cout << "Error in VTKInterface::readVTUFile() - Cell data not found." << std::endl;
			return NULL;
		}

		return new MeshLib::Mesh(BaseLib::getFileNameFromPath(file_name), nodes, elements);
	}
	else
		std::cout << "Error in VTKInterface::readVTUFile() - Number of nodes and elements not specified." << std::endl;
	
	return NULL;
}

MeshLib::Element* VTKInterface::readElement(std::stringstream &iss, const std::vector<MeshLib::Node*> &nodes, unsigned material, unsigned type)
{
	unsigned node_ids[8];
	switch (type)
	{
	case 3: //line
		for (unsigned i(0); i<2; i++) iss >> node_ids[i];
		return new MeshLib::Edge(nodes[node_ids[0]], nodes[node_ids[1]], material);
		break;
	case 5: //triangle
		for (unsigned i(0); i<3; i++) iss >> node_ids[i];
		return new MeshLib::Tri(nodes[node_ids[0]], nodes[node_ids[1]], nodes[node_ids[2]], material);
		break;
	case 9: //quad
		for (unsigned i(0); i<4; i++) iss >> node_ids[i];
		return new MeshLib::Quad(nodes[node_ids[0]], nodes[node_ids[1]], nodes[node_ids[2]], nodes[node_ids[3]], material);
		break;
	case 8: //pixel
		for (unsigned i(0); i<4; i++) iss >> node_ids[i];
		return new MeshLib::Quad(nodes[node_ids[0]], nodes[node_ids[1]], nodes[node_ids[3]], nodes[node_ids[2]], material);
		break;
	case 10:
		for (unsigned i(0); i<4; i++) iss >> node_ids[i];
		return new MeshLib::Tet(nodes[node_ids[0]], nodes[node_ids[1]], nodes[node_ids[2]], nodes[node_ids[3]], material);
		break;
	case 12: //hexahedron
		for (unsigned i(0); i<8; i++) iss >> node_ids[i];
		return new MeshLib::Hex(nodes[node_ids[0]], nodes[node_ids[1]], nodes[node_ids[2]], nodes[node_ids[3]],
				                nodes[node_ids[4]], nodes[node_ids[5]], nodes[node_ids[6]], nodes[node_ids[7]], material);
		break;
	case 11: //voxel
		for (unsigned i(0); i<8; i++) iss >> node_ids[i];
		return new MeshLib::Hex(nodes[node_ids[0]], nodes[node_ids[1]], nodes[node_ids[3]], nodes[node_ids[2]],
				                nodes[node_ids[4]], nodes[node_ids[5]], nodes[node_ids[7]], nodes[node_ids[6]], material);
		break;
	case 14: //pyramid
		for (unsigned i(0); i<5; i++) iss >> node_ids[i];
		return new MeshLib::Pyramid(nodes[node_ids[0]], nodes[node_ids[1]], nodes[node_ids[2]],
				                    nodes[node_ids[3]], nodes[node_ids[4]], material);
		break;
	case 13: //wedge
		for (unsigned i(0); i<6; i++) iss >> node_ids[i];
		return new MeshLib::Prism(nodes[node_ids[0]], nodes[node_ids[1]], nodes[node_ids[2]],
				                    nodes[node_ids[3]], nodes[node_ids[4]], nodes[node_ids[5]], material);
		break;
	default:
		std::cout << "Error in VTKInterface::readElement() - Unknown mesh element type \"" << type << "\" ..." << std::endl;
		return NULL;
	}

}

int VTKInterface::write(std::ostream& stream)
{
	return 0;
}


} // end namespace FileIO
