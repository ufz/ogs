/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file BoostVtuInterface.cpp
 *
 *  Created on 2012-12-05 by Karsten Rink
 */

#include "BoostVtuInterface.h"
#include <iostream>
#include <fstream>

#include <boost/foreach.hpp>

#include "StringTools.h"
#include "FileTools.h"
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

using namespace boost;

BoostVtuInterface::BoostVtuInterface()
: _export_name(""), _mesh(nullptr), /*_doc(),*/ _use_compressor(false)
{
}

BoostVtuInterface::~BoostVtuInterface()
{
//	delete _doc;
}

MeshLib::Mesh* BoostVtuInterface::readVTUFile(const std::string &file_name)
{
	std::cout << "Reading OGS mesh ... " << std::endl;
	std::ifstream in(file_name.c_str());
	if (in.fail())
	{
		std::cout << "\nRapidVtuInterface::readVTUFile() - Can't open xml-file." << std::endl;
		return nullptr;
	}

	// build DOM tree
	using boost::property_tree::ptree;
	ptree doc;
	read_xml(in, doc);

	//const ptree& v = doc.get_child("VTKFile");

	if (isVTKUnstructuredGrid(doc))
	{
		bool is_compressed = doc.get("VTKFile.<xmlattr>.compressor", false);
		if (is_compressed)
		{
			if (doc.get("VTKFile.<xmlattr>.compressor", "missing") == "vtkZLibDataCompressor")
			{
				//uncompressData(root_node);
			}
			else
			{
				std::cout << "BoostVtuInterface::readVTUFile() - Unknown compression method." << std::endl;
				return nullptr;
			}
		}

		//skip to <Piece>-tag and start parsing content
		optional<property_tree::ptree&> piece_node = doc.get_child_optional("VTKFile.UnstructuredGrid.Piece");
		if (piece_node)
		{
				
			const unsigned nNodes = static_cast<unsigned>(doc.get("VTKFile.UnstructuredGrid.Piece.<xmlattr>.NumberOfPoints", 0));
			const unsigned nElems = static_cast<unsigned>(doc.get("VTKFile.UnstructuredGrid.Piece.<xmlattr>.NumberOfCells", 0));

			if ((nNodes == 0) || (nElems == 0))
			{
				std::cout << "BoostVtuInterface::readVTUFile() - Number of nodes is " << nNodes 
					      << ", number of elements is " << nElems << "." << std::endl;
				return nullptr;
			}

			std::vector<MeshLib::Node*> nodes(nNodes);
			std::vector<MeshLib::Element*> elements(nElems);
			std::vector<unsigned> mat_ids(nElems, 0);
			std::vector<unsigned> cell_types(nElems);

			BOOST_FOREACH( ptree::value_type const& v, doc.get_child("VTKFile.UnstructuredGrid.Piece") ) 
			{
				if (v.first == "CellData")
				{
					std::string name (v.second.get("DataArray.<xmlattr>.Name", ""));
					std::string format (v.second.get("DataArray.<xmlattr>.format", "ascii"));
					if ((name.compare("MaterialIDs") == 0) || (name.compare("MatGroups") == 0))
					{
						std::stringstream iss (v.second.get<std::string>("DataArray"));
						if (format.compare("ascii") == 0)
						{
							for(unsigned i=0; i<nElems; i++)
								iss >> mat_ids[i];
						}
						else if (format.compare("appended") == 0)
						{
							//uncompress
						}
					}
				}
				
				if (v.first == "Points")
				{
					// This node may or may not have an attribute "Name" with the value "Points".
					// However, there shouldn't be any other DataArray nodes so most likely not checking the name isn't a problem.
					std::string format (v.second.get("DataArray.<xmlattr>.format", ""));
					if (format.compare("ascii") == 0)
					{
						std::stringstream iss (v.second.get<std::string>("DataArray"));
						double x,y,z;
						for(unsigned i=0; i<nNodes; i++)
						{
							iss >> x >> y >> z;
							nodes[i] = new MeshLib::Node(x,y,z,i);
						}
					}
					else if (format.compare("appended") == 0)
					{
						//uncompress
					}
				}

				if (v.first == "Cells")
				{
					std::string conn_string ("");
					BOOST_FOREACH( ptree::value_type const& c, doc.get_child("VTKFile.UnstructuredGrid.Piece.Cells") ) 
					{
						std::string attr_name (c.second.get<std::string>("<xmlattr>.Name"));
						if (attr_name.compare("connectivity") == 0)
						{
							conn_string = c.second.get<std::string>("DataArray");
							std::string format (c.second.get("DataArray.<xmlattr>.format", ""));
							if (format.compare("appended") == 0)
							{
								//uncompress
							}
						}
						if (attr_name.compare("types") == 0)
						{
							std::stringstream iss (c.second.get<std::string>("DataArray"));
							std::string format (c.second.get("DataArray.<xmlattr>.format", ""));
							if (format.compare("ascii") == 0)
							{
								for(unsigned i=0; i<nElems; i++)
									iss >> cell_types[i];
							}
							else if (format.compare("appended") == 0)
							{
								//uncompress
							}
						}
					}
					for(unsigned i=0; i<nElems; i++)
					{
						if (!conn_string.empty())
						{
							std::stringstream iss (conn_string);
							for(unsigned i=0; i<nElems; i++)
								elements[i] = readElement(iss, nodes, mat_ids[i], cell_types[i]);
						}
					}
				}


			}

			std::cout << "finished." << std::endl;
			std::cout << "Nr. Nodes: " << nodes.size() << std::endl;
			std::cout << "Nr. Elements: " << elements.size() << std::endl;
			return new MeshLib::Mesh(BaseLib::extractBaseNameWithoutExtension(file_name), nodes, elements);

		} // piece
	} // unstructured grid

	return nullptr;
}

MeshLib::Element* BoostVtuInterface::readElement(std::stringstream &iss, const std::vector<MeshLib::Node*> &nodes, unsigned material, unsigned type)
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
		std::cout << "Error in RapidVtuInterface::readElement() - Unknown mesh element type \"" << type << "\" ..." << std::endl;
		return nullptr;
	}

}

bool BoostVtuInterface::isVTKFile(const property_tree::ptree &vtk_root)
{
	if (!vtk_root.get_child_optional("VTKFile"))
	{
		std::cout << "Error in BoostVtuInterface::readVTUFile() - Not a VTK File." << std::endl;
		return false;
	}
	if (vtk_root.get("VTKFile.<xmlattr>.version", "missing") != "0.1")
	{
		std::cout << "Error in BoostVtuInterface::readVTUFile() - Unsupported file format version." << std::endl;
		return false;
	}
	if (vtk_root.get("VTKFile.<xmlattr>.byte_order", "missing") != "LittleEndian")
	{
		std::cout << "Error in BoostVtuInterface::readVTUFile() - Only little endian files are supported." << std::endl;
		return false;
	}
	return true;
}

bool BoostVtuInterface::isVTKUnstructuredGrid(const property_tree::ptree &vtk_root)
{
	if (isVTKFile(vtk_root))
	{
		optional<const property_tree::ptree&> u_grid_node = vtk_root.get_child_optional("VTKFile.UnstructuredGrid");
		if (u_grid_node)
			return true;
		std::cout << "Error in BoostVtuInterface::readVTUFile() - Not an unstructured grid." << std::endl;
	}
	return false;
}
/*
unsigned char* BoostVtuInterface::uncompressData(const rapidxml::xml_node<>* node)
{
	rapidxml::xml_node<>* data_node = node->first_node("AppendedData");
	char* compressed_data (NULL);
	if (data_node)
		compressed_data = data_node->value();

	return NULL;
}



int BoostVtuInterface::write(std::ostream& stream)
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
	const std::string str_nNodes (BaseLib::number2str(nNodes));
	piece_node->append_attribute (_doc->allocate_attribute("NumberOfPoints", str_nNodes.c_str()));
	const std::string str_nElems (BaseLib::number2str(nElems));
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
		typestream << this->getVTKElementID(element->getGeomType()) << " ";
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

unsigned BoostVtuInterface::getVTKElementID(MshElemType::type type) const
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

xml_node<>* BoostVtuInterface::addDataArray(const std::string &name, const std::string &data_type, const std::string &data, unsigned nComponents)
{

	xml_attribute<> *attr (NULL);
	xml_node<> *dataarray_node (_doc->allocate_node(node_element, _doc->allocate_string("DataArray"), _doc->allocate_string(data.c_str())));
	attr = _doc->allocate_attribute(_doc->allocate_string("type"), _doc->allocate_string(data_type.c_str()));
	dataarray_node->append_attribute(attr);
	attr = _doc->allocate_attribute(_doc->allocate_string("Name"), _doc->allocate_string(name.c_str()));
	dataarray_node->append_attribute(attr);
	if (nComponents > 1)
	{
		attr = _doc->allocate_attribute(_doc->allocate_string("NumberOfComponents"), _doc->allocate_string(BaseLib::number2str(nComponents).c_str()));
		dataarray_node->append_attribute(attr);
	}
	std::string comp_type = (_use_compressor) ? "appended" : "ascii";
	attr = _doc->allocate_attribute(_doc->allocate_string("format"), _doc->allocate_string(comp_type.c_str()));
	dataarray_node->append_attribute(attr);
	// ---- offset attribute for compressed data! ----
	return dataarray_node;
}
*/

} // end namespace FileIO

