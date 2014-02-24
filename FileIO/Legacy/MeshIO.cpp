/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-08
 * \brief  Implementation of the MeshIO class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * @file MeshIO.cpp
 * @date 2012-05-08
 * @author Karsten Rink
 */

#include <iomanip>
#include <sstream>

// ThirdParty/logog
#include "logog/include/logog.hpp"

// GeoLib
#include "GEOObjects.h"

// MeshLib
#include "Elements/Line.h"
#include "Elements/Hex.h"
#include "Elements/Prism.h"
#include "Elements/Pyramid.h"
#include "Elements/Quad.h"
#include "Elements/Tet.h"
#include "Elements/Tri.h"
#include "MeshIO.h"
#include "Node.h"

// BaseLib
#include "FileTools.h"
#include "StringTools.h"

namespace FileIO
{
namespace Legacy {

MeshIO::MeshIO()
	: _mesh(NULL)
{
}

MeshLib::Mesh* MeshIO::loadMeshFromFile(const std::string& file_name)
{
	INFO("Reading OGS legacy mesh ... ");

	std::ifstream in (file_name.c_str(),std::ios::in);
	if (!in.is_open())
	{
		WARN("MeshIO::loadMeshFromFile() - Could not open file %s.", file_name.c_str());
		return nullptr;
	}

	std::string line_string ("");
	getline(in, line_string);

	std::vector<MeshLib::Node*> nodes;
	std::vector<MeshLib::Element*> elements;

	if(line_string.find("#FEM_MSH") != std::string::npos) // OGS mesh file
	{
		while (!in.eof())
		{
			getline(in, line_string);

			// check keywords
			if (line_string.find("#STOP") != std::string::npos)
				break;
			else if (line_string.find("$NODES") != std::string::npos)
			{
				double x, y, z, double_dummy;
				unsigned idx;
				getline(in, line_string);
				BaseLib::trim(line_string);
				unsigned nNodes = atoi(line_string.c_str());
				std::string s;
				for (unsigned i = 0; i < nNodes; ++i)
				{
					getline(in, line_string);
					std::stringstream iss(line_string);
					iss >> idx >> x >> y >> z;
					MeshLib::Node* node(new MeshLib::Node(x, y, z, idx));
					nodes.push_back(node);
					iss >> s;
					if (s.find("$AREA") != std::string::npos)
						iss >> double_dummy;
				}
			}
			else if (line_string.find("$ELEMENTS") != std::string::npos)
			{
				getline(in, line_string);
				BaseLib::trim(line_string);
				unsigned nElements = atoi(line_string.c_str());
				for (unsigned i = 0; i < nElements; ++i)
				{
					getline(in, line_string);
					elements.push_back(readElement(line_string, nodes));
				}
			}
		}

		if (elements.empty())
		{
			ERR ("MeshIO::loadMeshFromFile() - File did not contain element information.");
			for (auto it = nodes.begin(); it!=nodes.end(); ++it)
				delete *it;
			return nullptr;
		}

		MeshLib::Mesh* mesh (new MeshLib::Mesh(BaseLib::extractBaseNameWithoutExtension(
		                                               file_name), nodes, elements));

		INFO("\t... finished.");
		INFO("Nr. Nodes: %d.", nodes.size());
		INFO("Nr. Elements: %d.", elements.size());

		in.close();
		return mesh;
	}
	else
	{
		in.close();
		return nullptr;
	}
}

MeshLib::Element* MeshIO::readElement(const std::string& line,
                                      const std::vector<MeshLib::Node*> &nodes)
{
	std::stringstream ss (line);
	std::string elem_type_str("");
	MeshElemType elem_type (MeshElemType::INVALID);
	unsigned index, patch_index;
	ss >> index >> patch_index;

	do {
		ss >> elem_type_str;
		if (ss.fail())
			return NULL;
		elem_type = String2MeshElemType(elem_type_str);
	} while (elem_type == MeshElemType::INVALID);

	unsigned* idx = new unsigned[8];
	MeshLib::Element* elem;

	switch(elem_type)
	{
	case MeshElemType::LINE: {
		for (int i = 0; i < 2; ++i)
			ss >> idx[i];
		// edge_nodes array will be deleted from Line object
		MeshLib::Node** edge_nodes = new MeshLib::Node*[2];
		for (unsigned k(0); k < 2; ++k)
			edge_nodes[k] = nodes[idx[k]];
		elem = new MeshLib::Line(edge_nodes, patch_index);
		break;
	}
	case MeshElemType::TRIANGLE: {
		for (int i = 0; i < 3; ++i)
			ss >> idx[i];
		MeshLib::Node** tri_nodes = new MeshLib::Node*[3];
		for (unsigned k(0); k < 3; ++k)
			tri_nodes[k] = nodes[idx[k]];
		elem = new MeshLib::Tri(tri_nodes, patch_index);
		break;
	}
	case MeshElemType::QUAD: {
		for (int i = 0; i < 4; ++i)
			ss >> idx[i];
		MeshLib::Node** quad_nodes = new MeshLib::Node*[4];
		for (unsigned k(0); k < 4; ++k)
			quad_nodes[k] = nodes[idx[k]];
		elem = new MeshLib::Quad(quad_nodes, patch_index);
		break;
	}
	case MeshElemType::TETRAHEDRON: {
		for (int i = 0; i < 4; ++i)
			ss >> idx[i];
		MeshLib::Node** tet_nodes = new MeshLib::Node*[4];
		for (unsigned k(0); k < 4; ++k)
			tet_nodes[k] = nodes[idx[k]];
		elem = new MeshLib::Tet(tet_nodes, patch_index);
		break;
	}
	case MeshElemType::HEXAHEDRON: {
		for (int i = 0; i < 8; ++i)
			ss >> idx[i];
		MeshLib::Node** hex_nodes = new MeshLib::Node*[8];
		for (unsigned k(0); k < 8; ++k)
			hex_nodes[k] = nodes[idx[k]];
		elem = new MeshLib::Hex(hex_nodes, patch_index);
		break;
	}
	case MeshElemType::PYRAMID: {
		for (int i = 0; i < 5; ++i)
			ss >> idx[i];
		MeshLib::Node** pyramid_nodes = new MeshLib::Node*[5];
		for (unsigned k(0); k < 5; ++k)
			pyramid_nodes[k] = nodes[idx[k]];
		elem = new MeshLib::Pyramid(pyramid_nodes, patch_index);
		break;
	}
	case MeshElemType::PRISM: {
		for (int i = 0; i < 6; ++i)
			ss >> idx[i];
		MeshLib::Node** prism_nodes = new MeshLib::Node*[6];
		for (unsigned k(0); k < 6; ++k)
			prism_nodes[k] = nodes[idx[k]];
		elem = new MeshLib::Prism(prism_nodes, patch_index);
		break;
	}
	default:
		elem = NULL;
		break;
	}

	delete [] idx;

	return elem;
}

bool MeshIO::write(std::ostream &out)
{
	if(!_mesh) {
		WARN("MeshIO::write(): Cannot write: no mesh object specified.");
		return false;
	}

	out << "#FEM_MSH\n"
		<< "$PCS_TYPE\n"
		<< "  NO_PCS\n"
		<< "$NODES\n"
		<< "  ";
	const size_t n_nodes(_mesh->getNNodes());
	out << n_nodes << "\n";
	for (size_t i(0); i < n_nodes; ++i) {
		out << i << " " << *(_mesh->getNode(i)) << "\n";
	}

	out << "$ELEMENTS\n"
		<< "  ";

	writeElements(_mesh->getElements(), out);

	out << " $LAYER\n"
		<< "  0\n"
		<< "#STOP\n";

	return true;
}

void MeshIO::setMesh(const MeshLib::Mesh* mesh)
{
	_mesh = mesh;
}

void MeshIO::writeElements(std::vector<MeshLib::Element*> const& ele_vec,
                                      std::ostream &out)
{
	const size_t ele_vector_size (ele_vec.size());

	out << ele_vector_size << "\n";
	for (size_t i(0); i < ele_vector_size; ++i) {
		out << i << " " << ele_vec[i]->getValue() << " " << MeshElemType2String(ele_vec[i]->getGeomType()) << " ";
		unsigned nElemNodes (ele_vec[i]->getNNodes());
		for(size_t j = 0; j < nElemNodes; ++j)
			out << ele_vec[i]->getNode(j)->getID() << " ";
		out << "\n";
	}
}

}
} // end namespace FileIO
