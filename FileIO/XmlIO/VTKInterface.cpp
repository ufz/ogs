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

#include "RapidXML/rapidxml_print.hpp"
#include "StringTools.h"
#include "ProjectData.h"

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


namespace FileIO {

using namespace rapidxml;

VTKInterface::VTKInterface()
: _export_name(""), _mesh(NULL), _doc(new xml_document<>), _use_compressor(false)
{
}

VTKInterface::~VTKInterface()
{
	delete _doc;
}

MeshLib::Mesh* VTKInterface::readVTUFile(const std::string &file_name)
{
	std::cout << "Reading OGS mesh ... " << std::endl;
	std::ifstream in(file_name.c_str());
	if (in.fail())
	{
		std::cout << "\nVTKInterface::readVTUFile() - Can't open xml-file." << std::endl;
		return NULL;
	}

	in.seekg(0, std::ios::end);
	const size_t length = in.tellg();
	in.seekg(0, std::ios::beg);
	char* buffer = new char[length+1];
	in.read(buffer, length);
	buffer[in.gcount()] = '\0';
	in.close();

	// build DOM tree
	rapidxml::xml_document<> doc;
	doc.parse<0>(buffer);

	rapidxml::xml_node<>* root_node (doc.first_node());
	if (isVTKUnstructuredGrid(root_node))
	{
		bool is_compressed(false);
		//check if content is compressed
		const rapidxml::xml_attribute<>* compressor (root_node->first_attribute("compressor"));
		if (compressor )
		{
			if (std::string(compressor->value()).compare("vtkZLibDataCompressor") == 0)
			{
				is_compressed = true;
				uncompressData(root_node);
			}
			else
			{
				std::cout << "VTKInterface::readVTUFile() - Unknown compression method." << std::endl;
				return NULL;
			}
		}

		//skip to <Piece>-tag and start parsing content
		const rapidxml::xml_node<>* piece_node (doc.first_node()->first_node()->first_node("Piece"));
		if (piece_node)
		{
			const unsigned nNodes = static_cast<unsigned>(atoi(piece_node->first_attribute("NumberOfPoints")->value()));
			const unsigned nElems = static_cast<unsigned>(atoi(piece_node->first_attribute("NumberOfCells")->value()));
			std::vector<MeshLib::Node*> nodes(nNodes);
			std::vector<MeshLib::Element*> elements(nElems);
			std::vector<unsigned> mat_ids(nElems, 0);
			std::vector<unsigned> cell_types(nElems);

			const rapidxml::xml_node<>* mat_id_node (piece_node->first_node("CellData")->first_node("DataArray"));
			if (mat_id_node &&
				((std::string(mat_id_node->first_attribute("Name")->value()).compare("MaterialIDs") == 0) ||
				 (std::string(mat_id_node->first_attribute("Name")->value()).compare("MatGroup") == 0)))
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
			// This _may_ have an attribute "Name" with the value "Points" but you cannot count on it.
			// However, there shouldn't be any other DataArray nodes so most likely not checking the name isn't a problem.
			if (points_node)
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
				const std::string format = std::string(celltype_node->first_attribute("format")->value());
				if (format.compare("ascii") == 0)
				{
					std::stringstream iss (celltype_node->value());
					for(unsigned i=0; i<nElems; i++)
						iss >> cell_types[i];
				}
				else if (format.compare("appended") == 0)
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

			std::cout << "finished." << std::endl;
			std::cout << "Nr. Nodes: " << nodes.size() << std::endl;
			std::cout << "Nr. Elements: " << elements.size() << std::endl;
			delete [] buffer;
			return new MeshLib::Mesh(BaseLib::getFileNameFromPath(file_name), nodes, elements);
		}
		else {
			std::cout << "Error in VTKInterface::readVTUFile() - Number of nodes and elements not specified." << std::endl;
			delete [] buffer;
		}
	}
	delete [] buffer;
	return NULL;
}

MeshLib::Element* VTKInterface::readElement(std::stringstream &iss, const std::vector<MeshLib::Node*> &nodes, unsigned material, unsigned type)
{
	unsigned node_ids[8];
	switch (type)
	{
	case 3: { //line
		for (unsigned i(0); i<2; i++) iss >> node_ids[i];
		MeshLib::Node** edge_nodes = new MeshLib::Node*[2];
		edge_nodes[0] = nodes[node_ids[0]];
		edge_nodes[1] = nodes[node_ids[1]];
		return new MeshLib::Edge(edge_nodes, material);
		break;
	}
	case 5: { //triangle
		for (unsigned i(0); i<3; i++) iss >> node_ids[i];
		MeshLib::Node** tri_nodes = new MeshLib::Node*[3];
		tri_nodes[0] = nodes[node_ids[0]];
		tri_nodes[1] = nodes[node_ids[1]];
		tri_nodes[2] = nodes[node_ids[2]];
		return new MeshLib::Tri(tri_nodes, material);
		break;
	}
	case 9: { //quad
		for (unsigned i(0); i<4; i++) iss >> node_ids[i];
		MeshLib::Node** quad_nodes = new MeshLib::Node*[4];
		for (unsigned k(0); k<4; k++)
			quad_nodes[k] = nodes[node_ids[k]];
		return new MeshLib::Quad(quad_nodes, material);
		break;
	}
	case 8: { //pixel
		for (unsigned i(0); i<4; i++) iss >> node_ids[i];
		MeshLib::Node** quad_nodes = new MeshLib::Node*[4];
		quad_nodes[0] = nodes[node_ids[0]];
		quad_nodes[1] = nodes[node_ids[1]];
		quad_nodes[2] = nodes[node_ids[3]];
		quad_nodes[3] = nodes[node_ids[2]];
		return new MeshLib::Quad(quad_nodes, material);
		break;
	}
	case 10: {
		for (unsigned i(0); i<4; i++) iss >> node_ids[i];
		MeshLib::Node** tet_nodes = new MeshLib::Node*[4];
		for (unsigned k(0); k<4; k++)
			tet_nodes[k] = nodes[node_ids[k]];
		return new MeshLib::Tet(tet_nodes, material);
		break;
	}
	case 12: { //hexahedron
		for (unsigned i(0); i<8; i++) iss >> node_ids[i];
		MeshLib::Node** hex_nodes = new MeshLib::Node*[8];
		for (unsigned k(0); k<8; k++)
			hex_nodes[k] = nodes[node_ids[k]];
		return new MeshLib::Hex(hex_nodes, material);
		break;
	}
	case 11: { //voxel
		for (unsigned i(0); i<8; i++) iss >> node_ids[i];
		MeshLib::Node** voxel_nodes = new MeshLib::Node*[8];
		voxel_nodes[0] = nodes[node_ids[0]];
		voxel_nodes[1] = nodes[node_ids[1]];
		voxel_nodes[2] = nodes[node_ids[3]];
		voxel_nodes[3] = nodes[node_ids[2]];
		voxel_nodes[4] = nodes[node_ids[4]];
		voxel_nodes[5] = nodes[node_ids[5]];
		voxel_nodes[6] = nodes[node_ids[7]];
		voxel_nodes[7] = nodes[node_ids[6]];
		return new MeshLib::Hex(voxel_nodes, material);
		break;
	}
	case 14: { //pyramid
		for (unsigned i(0); i<5; i++) iss >> node_ids[i];
		MeshLib::Node** pyramid_nodes = new MeshLib::Node*[5];
		for (unsigned k(0); k<5; k++)
			pyramid_nodes[k] = nodes[node_ids[k]];
		return new MeshLib::Pyramid(pyramid_nodes, material);
		break;
	}
	case 13: { //wedge
		for (unsigned i(0); i<6; i++) iss >> node_ids[i];
		MeshLib::Node** prism_nodes = new MeshLib::Node*[6];
		for (unsigned k(0); k<6; k++)
			prism_nodes[k] = nodes[node_ids[k]];
		return new MeshLib::Prism(prism_nodes, material);
		break;
	}
	default:
		std::cout << "Error in VTKInterface::readElement() - Unknown mesh element type \"" << type << "\" ..." << std::endl;
		return NULL;
	}

}

bool VTKInterface::isVTKFile(const rapidxml::xml_node<>* vtk_root)
{
	if (vtk_root == NULL || std::string(vtk_root->name()).compare("VTKFile"))
	{
		std::cout << "Error in VTKInterface::readVTUFile() - Not a VTK File." << std::endl;
		return false;
	}
	if (std::string(vtk_root->first_attribute("version")->value()).compare("0.1"))
	{
		std::cout << "Error in VTKInterface::readVTUFile() - Unsupported file format version." << std::endl;
		return false;
	}
	if (std::string(vtk_root->first_attribute("byte_order")->value()).compare("LittleEndian"))
	{
		std::cout << "Error in VTKInterface::readVTUFile() - Only little endian files are supported." << std::endl;
		return false;
	}
	return true;
}

bool VTKInterface::isVTKUnstructuredGrid(const rapidxml::xml_node<>* vtk_root)
{
	if (isVTKFile(vtk_root))
	{
		if (std::string(vtk_root->first_node()->name()).compare("UnstructuredGrid") == 0)
			return true;
		std::cout << "Error in VTKInterface::readVTUFile() - Not an unstructured grid." << std::endl;
	}
	return false;
}

unsigned char* VTKInterface::uncompressData(const rapidxml::xml_node<>* node)
{
	rapidxml::xml_node<>* data_node = node->first_node("AppendedData");
	char* compressed_data (NULL);
	if (data_node)
		compressed_data = data_node->value();

	return NULL;
}



int VTKInterface::write(std::ostream& stream)
{
	//if (this->_export_name.empty())
	if (!_mesh)
	{
		std::cout << "Error in XmlStnInterface::write() - No station list specified..." << std::endl;
		return 0;
	}

	const size_t nNodes (_mesh->getNNodes());
	const size_t nElems (_mesh->getNElements());
	const std::vector<MeshLib::Node*> &nodes (_mesh->getNodes());
	const std::vector<MeshLib::Element*> &elements (_mesh->getElements());

	const std::string data_array_close("\t\t\t\t");
	const std::string data_array_indent("\t\t\t\t  ");

	stream << "<?xml version=\"1.0\"?>" << std::endl;

	xml_node<> *root_node (_doc->allocate_node(node_element, "VTKFile"));
	_doc->append_node(root_node);
	root_node->append_attribute(_doc->allocate_attribute("type", "UnstructuredGrid"));
	root_node->append_attribute(_doc->allocate_attribute("version", "0.1"));
	root_node->append_attribute(_doc->allocate_attribute("byte_order", "LittleEndian"));
	if (_use_compressor)
		root_node->append_attribute(_doc->allocate_attribute("compressor", "vtkZLibDataCompressor"));

	xml_node<> *grid_node (_doc->allocate_node(node_element, "UnstructuredGrid"));
	root_node->append_node(grid_node);
	xml_node<> *piece_node (_doc->allocate_node(node_element, "Piece"));
	grid_node->append_node(piece_node);
	const std::string str_nNodes (number2str(nNodes));
	piece_node->append_attribute (_doc->allocate_attribute("NumberOfPoints", str_nNodes.c_str()));
	const std::string str_nElems (number2str(nElems));
	piece_node->append_attribute(_doc->allocate_attribute("NumberOfCells", str_nElems.c_str()));

	// scalar arrays for point- and cell-data
	xml_node<> *pointdata_node (_doc->allocate_node(node_element, "PointData", "\n\t\t\t"));
	piece_node->append_node(pointdata_node);
	// add node_area array here if necessary!
	xml_node<> *celldata_node (_doc->allocate_node(node_element, "CellData"));
	piece_node->append_node(celldata_node);
	celldata_node->append_attribute(_doc->allocate_attribute("Scalars", "MaterialIDs"));

	std::stringstream oss(std::stringstream::out);
	oss << std::endl << data_array_indent;
	for (unsigned i=0; i<nElems; i++)
		oss << elements[i]->getValue() << " ";
	oss << std::endl << data_array_close;
	celldata_node->append_node(this->addDataArray("MaterialIDs", "Int32", oss.str()));
	oss.str(std::string());
	oss.clear();

	// point coordinates
	xml_node<> *points_node (_doc->allocate_node(node_element, "Points"));
	piece_node->append_node(points_node);
	oss << std::endl;
	for (unsigned i=0; i<nNodes; i++)
		oss << data_array_indent << (*nodes[i])[0] << " " << (*nodes[i])[1] << " " << (*nodes[i])[2] << std::endl;
	oss << data_array_close;
	points_node->append_node(this->addDataArray("Points", "Float32", oss.str(), 3));
	oss.str(std::string());
	oss.clear();

	// cells with node ids
	xml_node<> *cells_node (_doc->allocate_node(node_element, "Cells"));
	piece_node->append_node(cells_node);
	std::stringstream offstream(std::stringstream::out);
	std::stringstream typestream(std::stringstream::out);
	oss << std::endl;
	offstream << std::endl << data_array_indent;
	typestream << std::endl << data_array_indent;

	unsigned offset_count(0);
	for (unsigned i=0; i<nElems; i++)
	{
		MeshLib::Element* element (elements[i]);
		const unsigned nElemNodes (element->getNNodes());
		oss << data_array_indent;
		for (unsigned j=0; j<nElemNodes; j++)
			oss << element->getNode(j)->getID() << " ";
		oss << std::endl;
		offset_count += nElemNodes;
		offstream << offset_count << " ";
		typestream << this->getVTKElementID(element->getGeoType()) << " ";
	}
	oss << data_array_close;
	offstream << std::endl << data_array_close;
	typestream << std::endl << data_array_close;

	// connectivity attributes
	cells_node->append_node(this->addDataArray("connectivity", "Int32", oss.str()));
	cells_node->append_node(this->addDataArray("offsets", "Int32", offstream.str()));
	cells_node->append_node(this->addDataArray("types", "UInt8", typestream.str()));

	stream << *_doc;
	return 1;
}

unsigned VTKInterface::getVTKElementID(MshElemType::type type) const
{
	switch (type)
	{
	case MshElemType::EDGE:
		return 3;
	case MshElemType::TRIANGLE:
		return 5;
	case MshElemType::QUAD:
		return 9;
	case MshElemType::TETRAHEDRON:
		return 10;
	case MshElemType::HEXAHEDRON:
		return 12;
	case MshElemType::PYRAMID:
		return 14;
	case MshElemType::PRISM:
		return 13;
	default:
		return std::numeric_limits<unsigned>::max();
	}
}

xml_node<>* VTKInterface::addDataArray(const std::string &name, const std::string &data_type, const std::string &data, unsigned nComponents)
{

	xml_attribute<> *attr (NULL);
	xml_node<> *dataarray_node (_doc->allocate_node(node_element, _doc->allocate_string("DataArray"), _doc->allocate_string(data.c_str())));
	attr = _doc->allocate_attribute(_doc->allocate_string("type"), _doc->allocate_string(data_type.c_str()));
	dataarray_node->append_attribute(attr);
	attr = _doc->allocate_attribute(_doc->allocate_string("Name"), _doc->allocate_string(name.c_str()));
	dataarray_node->append_attribute(attr);
	if (nComponents > 1)
	{
		attr = _doc->allocate_attribute(_doc->allocate_string("NumberOfComponents"), _doc->allocate_string(number2str(nComponents).c_str()));
		dataarray_node->append_attribute(attr);
	}
	std::string comp_type = (_use_compressor) ? "appended" : "ascii";
	attr = _doc->allocate_attribute(_doc->allocate_string("format"), _doc->allocate_string(comp_type.c_str()));
	dataarray_node->append_attribute(attr);
	// ---- offset attribute for compressed data! ----
	return dataarray_node;
}

} // end namespace FileIO
