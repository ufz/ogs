/**
 * \file
 * \author Karsten Rink
 * \date   2012-12-05
 * \brief  Implementation of the BoostVtuInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file BoostVtuInterface.cpp
 *  @date 2012-12-05
 *  @author Karsten Rink
 *  @brief Read VTU files employing boost.
 */

#include "BoostVtuInterface.h"
#include "zLibDataCompressor.h"
#include <fstream>

#include <boost/version.hpp>
#include <boost/foreach.hpp>
#include <boost/property_tree/xml_parser.hpp>

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "FileTools.h"

// MSH
#include "Elements/Line.h"
#include "Elements/Hex.h"
#include "Elements/Prism.h"
#include "Elements/Pyramid.h"
#include "Elements/Quad.h"
#include "Elements/Tet.h"
#include "Elements/Tri.h"
#include "Mesh.h"
#include "MeshLib/Node.h"

namespace FileIO
{
using namespace boost;

BoostVtuInterface::BoostVtuInterface() :
	_mesh(nullptr), _use_compressor(false), _doc()
{
}

BoostVtuInterface::~BoostVtuInterface()
{}

MeshLib::Mesh* BoostVtuInterface::readVTUFile(const std::string &file_name)
{
	std::ifstream in(file_name.c_str());
	if (in.fail())
	{
		ERR("BoostVtuInterface::readVTUFile(): Can't open xml-file %s.", file_name.c_str());
		return nullptr;
	}

	// build DOM tree
	using boost::property_tree::ptree;
	ptree doc;
	read_xml(in, doc);

	if (isVTKUnstructuredGrid(doc))
	{
		ptree const& root_node = doc.get_child("VTKFile");
		optional<std::string> const& compressor (getXmlAttribute("compressor", root_node));
		bool is_compressed = static_cast<bool>(compressor);
		if (is_compressed)
		{
			if (*compressor != "vtkZLibDataCompressor")
			{
				ERR("BoostVtuInterface::readVTUFile(): Unknown compression method.");
				return nullptr;
			}

			// TODO: remove this once compressed data can be handled!!
			INFO("Handling of compressed meshes not yet implemented.");
			return nullptr;
		}

		//skip to <Piece>-tag and start parsing content
		OptionalPtree const& piece_node = root_node.get_child_optional(
		        "UnstructuredGrid.Piece");
		if (piece_node)
		{
			const unsigned nNodes =
			        static_cast<unsigned>(piece_node->get("<xmlattr>.NumberOfPoints", 0));
			const unsigned nElems =
			        static_cast<unsigned>(piece_node->get("<xmlattr>.NumberOfCells", 0));

			if ((nNodes == 0) || (nElems == 0))
			{
				ERR("BoostVtuInterface::readVTUFile() - Number of nodes is %d, number of elements is %d.",
				    nNodes, nElems);
				return nullptr;
			}

			std::vector<MeshLib::Node*> nodes(nNodes);
			std::vector<MeshLib::Element*> elements(nElems);
			std::vector<unsigned> mat_ids(nElems, 0);
			std::vector<unsigned> cell_types(nElems);

			BOOST_FOREACH( ptree::value_type const & grid_piece, *piece_node )
			{
				if (grid_piece.first == "CellData")
				{
					const OptionalPtree& cell_data_node = findDataArray(
					        "MaterialIDs",
					        grid_piece.second);
					if (cell_data_node)
					{
						optional<std::string> const& format =
						        getXmlAttribute("format", *cell_data_node);
						std::stringstream iss (cell_data_node->data()); //v.second.get<std::string>("DataArray"));
						if (format)
						{
							if (*format == "ascii")
								for(unsigned i = 0; i < nElems; i++)
									iss >> mat_ids[i];
							else if (*format == "appended")
							{
								//uncompress
							}
						}
					}
					else
						WARN("BoostVtuInterface::readVTUFile(): MaterialIDs not found, setting every cell to 0.");
				}

				if (grid_piece.first == "Points")
				{
					// This node may or may not have an attribute "Name" with the value "Points".
					// However, there shouldn't be any other DataArray nodes so most likely not checking the name isn't a problem.
					ptree const& data_array_node = grid_piece.second.get_child(
					        "DataArray");
					optional<std::string> const& format = getXmlAttribute(
					        "format",
					        data_array_node);

					if (format)
					{
						if (*format == "ascii")
						{
							std::stringstream iss (data_array_node.data());
							double x,y,z;
							for(unsigned i = 0; i < nNodes; i++)
							{
								iss >> x >> y >> z;
								nodes[i] = new MeshLib::Node(x, y, z, i);
							}
						}
						else if (*format == "appended")
						{
							//uncompress
						}
					}
				}

				if (grid_piece.first == "Cells")
				{
					ptree const& cells = grid_piece.second;

					// cell types
					OptionalPtree const& types = findDataArray("types", cells);
					if (!types)
						ERR("BoostVtuInterface::readVTUFile(): Cannot find \"types\" data array.");

					std::stringstream iss (types->data());
					optional<std::string> const& format = getXmlAttribute("format", *types);
					if (*format == "ascii")
					{
						for(unsigned i = 0; i < nElems; i++)
							iss >> cell_types[i];
					}
					else if (*format == "appended")
					{
						//uncompress
					}

					// connectivity / element nodes
					OptionalPtree const& connectivity = findDataArray("connectivity", cells);
					if (!connectivity)
						ERR("BoostVtuInterface::readVTUFile(): Cannot find \"connectivity\" data array.");

					std::string conn_string = connectivity->data();

					if (!conn_string.empty())
					{
						optional<std::string> const& format =
						        getXmlAttribute("format", *connectivity);
						if (*format == "appended")
						{
							//uncompress
						}

						std::stringstream iss (conn_string);
						for(unsigned i = 0; i < nElems; i++)
							elements[i] = readElement(iss, nodes, mat_ids[i],cell_types[i]);
					}
				}
			}

			INFO("Reading VTU mesh finished.");
			INFO("Name: %s", BaseLib::extractBaseNameWithoutExtension(file_name).c_str());
			INFO("Nr. Nodes: %d", nodes.size());
			INFO("Nr. Elements: %d", elements.size());
			return new MeshLib::Mesh(BaseLib::extractBaseNameWithoutExtension(file_name), nodes,
			                         elements);

		} // piece
	} // unstructured grid

	return nullptr;
}

MeshLib::Element* BoostVtuInterface::readElement(std::stringstream &iss,
                                                 const std::vector<MeshLib::Node*> &nodes,
                                                 unsigned material, unsigned type)
{
	unsigned node_ids[8];
	switch (type)
	{
	case 3: { //line
		for (unsigned i(0); i < 2; i++)
			iss >> node_ids[i];
		MeshLib::Node** edge_nodes = new MeshLib::Node*[2];
		edge_nodes[0] = nodes[node_ids[0]];
		edge_nodes[1] = nodes[node_ids[1]];
		return new MeshLib::Line(edge_nodes, material);
		break;
	}
	case 5: { //triangle
		for (unsigned i(0); i < 3; i++)
			iss >> node_ids[i];
		MeshLib::Node** tri_nodes = new MeshLib::Node*[3];
		tri_nodes[0] = nodes[node_ids[0]];
		tri_nodes[1] = nodes[node_ids[1]];
		tri_nodes[2] = nodes[node_ids[2]];
		return new MeshLib::Tri(tri_nodes, material);
		break;
	}
	case 9: { //quad
		for (unsigned i(0); i < 4; i++)
			iss >> node_ids[i];
		MeshLib::Node** quad_nodes = new MeshLib::Node*[4];
		for (unsigned k(0); k < 4; k++)
			quad_nodes[k] = nodes[node_ids[k]];
		return new MeshLib::Quad(quad_nodes, material);
		break;
	}
	case 8: { //pixel
		for (unsigned i(0); i < 4; i++)
			iss >> node_ids[i];
		MeshLib::Node** quad_nodes = new MeshLib::Node*[4];
		quad_nodes[0] = nodes[node_ids[0]];
		quad_nodes[1] = nodes[node_ids[1]];
		quad_nodes[2] = nodes[node_ids[3]];
		quad_nodes[3] = nodes[node_ids[2]];
		return new MeshLib::Quad(quad_nodes, material);
		break;
	}
	case 10: {
		for (unsigned i(0); i < 4; i++)
			iss >> node_ids[i];
		MeshLib::Node** tet_nodes = new MeshLib::Node*[4];
		for (unsigned k(0); k < 4; k++)
			tet_nodes[k] = nodes[node_ids[k]];
		return new MeshLib::Tet(tet_nodes, material);
		break;
	}
	case 12: { //hexahedron
		for (unsigned i(0); i < 8; i++)
			iss >> node_ids[i];
		MeshLib::Node** hex_nodes = new MeshLib::Node*[8];
		for (unsigned k(0); k < 8; k++)
			hex_nodes[k] = nodes[node_ids[k]];
		return new MeshLib::Hex(hex_nodes, material);
		break;
	}
	case 11: { //voxel
		for (unsigned i(0); i < 8; i++)
			iss >> node_ids[i];
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
		for (unsigned i(0); i < 5; i++)
			iss >> node_ids[i];
		MeshLib::Node** pyramid_nodes = new MeshLib::Node*[5];
		for (unsigned k(0); k < 5; k++)
			pyramid_nodes[k] = nodes[node_ids[k]];
		return new MeshLib::Pyramid(pyramid_nodes, material);
		break;
	}
	case 13: { //wedge
		for (unsigned i(0); i < 6; i++)
			iss >> node_ids[i];
		MeshLib::Node** prism_nodes = new MeshLib::Node*[6];
		for (unsigned k(0); k < 6; k++)
			prism_nodes[k] = nodes[node_ids[k]];
		return new MeshLib::Prism(prism_nodes, material);
		break;
	}
	default:
		ERR("BoostVtuInterface::readElement(): Unknown mesh element type \"%d\".", type);
		return nullptr;
	}
}

bool BoostVtuInterface::isVTKFile(const property_tree::ptree &vtk_root)
{
	if (!vtk_root.get_child_optional("VTKFile"))
	{
		ERR("BoostVtuInterface::isVTKFile(): Not a VTK file.");
		return false;
	}
	optional<std::string> const& att_version (getXmlAttribute("version", vtk_root));
	if (att_version && *att_version == "0.1")
	{
		ERR("BoostVtuInterface::isVTKFile(): Unsupported file format version.");
		return false;
	}
	optional<std::string> const& att_order (getXmlAttribute("byte_order", vtk_root));
	if (att_order && *att_order == "LittleEndian")
	{
		ERR("BoostVtuInterface::isVTKFile(): Only little endian files are supported.");
		return false;
	}
	return true;
}

bool BoostVtuInterface::isVTKUnstructuredGrid(const property_tree::ptree &vtk_root)
{
	if (isVTKFile(vtk_root))
	{
		const OptionalPtree &u_grid_node = vtk_root.get_child_optional(
		        "VTKFile.UnstructuredGrid");
		if (u_grid_node)
			return true;
		ERR("Error in BoostVtuInterface::isVTKUnstructuredGrid(): Not an unstructured grid.");
	}
	return false;
}

unsigned char* BoostVtuInterface::uncompressData(property_tree::ptree const& compressed_data_node)
{
	const unsigned char* compressed_data = reinterpret_cast<const unsigned char*>(compressed_data_node.data().c_str());
	unsigned long compressed_size = strlen(compressed_data_node.data().c_str());
	unsigned char* uncompressed_data = NULL;
	unsigned long uncompressed_size = 0;
	// unsigned long result =
	zLibDataCompressor::UncompressBuffer(compressed_data, compressed_size, uncompressed_data, uncompressed_size);
	return uncompressed_data;
}

const optional<std::string> BoostVtuInterface::getXmlAttribute(std::string const& key,
                                                               property_tree::ptree const& tree)
{
	for (property_tree::ptree::const_iterator it = tree.begin(); it != tree.end(); ++it)
	{
		if (it->first != "<xmlattr>")
			continue;
		if (it->second.get_child_optional(key))
			return it->second.get_child(key).data();
	}

	return optional<std::string>();
}

const OptionalPtree BoostVtuInterface::findDataArray(std::string const& array_name,
                                                     property_tree::ptree const& tree)
{
	// Loop over all "DataArray" children.
	for (property_tree::ptree::const_iterator it = tree.begin(); it != tree.end(); ++it)
		if (it->first == "DataArray")
		{
			optional<std::string> const& value = getXmlAttribute("Name", it->second);
			if (value && *value == array_name)
				return it->second;
		}

	return OptionalPtree();
}

void BoostVtuInterface::setMesh(const MeshLib::Mesh* mesh)
{
	if (!mesh)
	{
		ERR("BoostVtuInterface::write(): No mesh specified.");
		return;
	}
	this->_mesh = const_cast<MeshLib::Mesh*>(mesh);
	buildPropertyTree();
};

void BoostVtuInterface::addScalarPointProperty(std::string const& name,
	std::vector<double> const& prop_vals)
{
	if (_doc.empty()) {
		WARN("BoostVtuInterface::addPointProperty(): propertry tree empty (no mesh set)");
		return;
	}

	if (_mesh->getNNodes() != prop_vals.size()) {
		WARN("BoostVtuInterface::addPointProperty(): number of values for propertry %s (%d) does not match the number of nodes (%d)", name.c_str(), prop_vals.size(), _mesh->getNNodes());
		return;
	}

#if BOOST_VERSION <= 105500
	const std::string data_array_close("\t\t\t\t");
	const std::string data_array_indent("\t\t\t\t  ");
#endif

	// go to the node where data should be inserted
	using boost::property_tree::ptree;
	ptree &root_node = _doc.get_child("VTKFile");
	ptree &piece_node = root_node.get_child("UnstructuredGrid.Piece");
	ptree &pointdata_node = piece_node.get_child("PointData");

	// prepare the data
	std::stringstream oss(std::stringstream::out);
	oss.precision(_out.precision());
#if BOOST_VERSION <= 105500
	oss << "\n" << data_array_indent;
	std::copy(prop_vals.cbegin(), prop_vals.cend(), std::ostream_iterator<double>(oss, " "));
	oss << "\n" << data_array_close;
#else
	std::copy(prop_vals.cbegin(), prop_vals.cend(), std::ostream_iterator<double>(oss, " "));
#endif
	this->addDataArray(pointdata_node, name, "Float64", oss.str());
	oss.str(std::string());
}

void BoostVtuInterface::buildPropertyTree()
{
	_doc.clear();
	const std::size_t nNodes (_mesh->getNNodes());
	const std::size_t nElems (_mesh->getNElements());
	const std::vector<MeshLib::Node*> &nodes (_mesh->getNodes());
	const std::vector<MeshLib::Element*> &elements (_mesh->getElements());

#if BOOST_VERSION <= 105500
	const std::string data_array_close("\t\t\t\t");
	const std::string data_array_indent("\t\t\t\t  ");
#endif

	using boost::property_tree::ptree;

	ptree &root_node = _doc.put("VTKFile", "");
	root_node.put("<xmlattr>.type", "UnstructuredGrid");
	root_node.put("<xmlattr>.version", "0.1");
	root_node.put("<xmlattr>.byte_order", "LittleEndian");

	if (_use_compressor)
		root_node.put("<xmlattr>.compressor", "vtkZLibDataCompressor");

	ptree &piece_node = root_node.put("UnstructuredGrid.Piece", "");
	const std::string str_nNodes (std::to_string(nNodes));
	const std::string str_nElems (std::to_string(nElems));
	piece_node.put("<xmlattr>.NumberOfPoints", str_nNodes.c_str());
	piece_node.put("<xmlattr>.NumberOfCells", str_nElems.c_str());

	// scalar arrays for point- and cell-data
	piece_node.add("PointData", "");
	// add node_area array here if necessary!
	ptree &celldata_node = piece_node.add("CellData", "");
	celldata_node.put("<xmlattr>.Scalars", "MaterialIDs");

	std::stringstream oss(std::stringstream::out);
	oss.precision(_out.precision());
#if BOOST_VERSION <= 105500
	oss << std::endl << data_array_indent;
	for (unsigned i = 0; i < nElems; i++)
		oss << elements[i]->getValue() << " ";
	oss << std::endl << data_array_close;
#else
	for (unsigned i = 0; i < nElems; i++)
		oss << elements[i]->getValue() << " ";
#endif
	this->addDataArray(celldata_node, "MaterialIDs", "Int32", oss.str());
	oss.str(std::string());
	oss.clear();

	// point coordinates
	ptree &points_node = piece_node.add("Points", "");
#if BOOST_VERSION <= 105500
	oss << std::endl;
	for (unsigned i = 0; i < nNodes; i++)
		oss << data_array_indent << (*nodes[i])[0] << " " << (*nodes[i])[1] << " " <<
		(*nodes[i])[2] << std::endl;
	oss << data_array_close;
#else
	for (unsigned i = 0; i < nNodes; i++)
		oss << (*nodes[i])[0] << " " << (*nodes[i])[1] << " " << (*nodes[i])[2] << " ";
#endif
	this->addDataArray(points_node, "Points", "Float64", oss.str(), 3);
	oss.str(std::string());
	oss.clear();

	// cells with node ids
	ptree &cells_node = piece_node.add("Cells", "");
	std::stringstream offstream(std::stringstream::out);
	std::stringstream typestream(std::stringstream::out);
#if BOOST_VERSION <= 105500
	oss << std::endl;
	offstream << std::endl << data_array_indent;
	typestream << std::endl << data_array_indent;
#endif
	unsigned offset_count(0);
	for (unsigned i = 0; i < nElems; i++)
	{
		MeshLib::Element* element (elements[i]);
		const unsigned nElemNodes (element->getNBaseNodes());
#if BOOST_VERSION <= 105500
		oss << data_array_indent;
#endif
		for (unsigned j = 0; j < nElemNodes; j++)
			oss << element->getNode(j)->getID() << " ";
#if BOOST_VERSION <= 105500
		oss << std::endl;
#endif
		offset_count += nElemNodes;
		offstream << offset_count << " ";
		typestream << this->getVTKElementID(element->getGeomType()) << " ";
	}
#if BOOST_VERSION <= 105500
	oss << data_array_close;
	offstream << std::endl << data_array_close;
	typestream << std::endl << data_array_close;
#endif

	// connectivity attributes
	this->addDataArray(cells_node, "connectivity", "Int32", oss.str());
	this->addDataArray(cells_node, "offsets", "Int32", offstream.str());
	this->addDataArray(cells_node, "types", "UInt8", typestream.str());
}

bool BoostVtuInterface::write()
{
	if (_doc.empty()) {
		ERR("BoostVtuInterface::write(): No mesh specified.");
		return false;
	}
#if BOOST_VERSION <= 105500
	property_tree::xml_writer_settings<char> settings('\t', 1);
#else
	property_tree::xml_writer_settings<std::string> settings('\t', 1);
#endif  // BOOST_VERSION
	write_xml(_out, _doc, settings);
	return true;
}

unsigned BoostVtuInterface::getVTKElementID(MeshElemType type) const
{
	switch (type)
	{
	case MeshElemType::LINE:
		return 3;
	case MeshElemType::TRIANGLE:
		return 5;
	case MeshElemType::QUAD:
		return 9;
	case MeshElemType::TETRAHEDRON:
		return 10;
	case MeshElemType::HEXAHEDRON:
		return 12;
	case MeshElemType::PYRAMID:
		return 14;
	case MeshElemType::PRISM:
		return 13;
	default:
		return std::numeric_limits<unsigned>::max();
	}
}

void BoostVtuInterface::addDataArray(property_tree::ptree &parent_node, const std::string &name,
                                     const std::string &data_type, const std::string &data,
                                     unsigned nComponents)
{
	property_tree::ptree &dataarray_node = parent_node.add("DataArray", data.c_str());
	dataarray_node.put("<xmlattr>.type", data_type.c_str());
	dataarray_node.put("<xmlattr>.Name", name.c_str());
	if (nComponents > 1)
		dataarray_node.put("<xmlattr>.NumberOfComponents", std::to_string(nComponents).c_str());
	std::string comp_type = (_use_compressor) ? "appended" : "ascii";
	dataarray_node.put("<xmlattr>.format", comp_type.c_str());
	// ---- offset attribute for compressed data! ----
}
} // end namespace FileIO
