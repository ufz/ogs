/**
 * \file
 * \author Thomas Fischer
 * \date   2011-09-12
 * \brief  Implementation of the TetGenInterface class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <cstddef>
#include <string>

// BaseLib
#include "FileTools.h"
#include "StringTools.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// FileIO
#include "TetGenInterface.h"

// MeshLib
#include "Mesh.h"
#include "Node.h"
#include "Elements/Element.h"
#include "Elements/Tet.h"

namespace FileIO
{
TetGenInterface::TetGenInterface() :
	_zero_based_idx (false)
{
}

TetGenInterface::~TetGenInterface()
{
}

MeshLib::Mesh* TetGenInterface::readTetGenMesh (std::string const& nodes_fname,
                                                std::string const& ele_fname)
{
	std::ifstream ins_nodes (nodes_fname.c_str());
	std::ifstream ins_ele (ele_fname.c_str());

	if (!ins_nodes || !ins_ele)
	{
		if (!ins_nodes)
			ERR ("TetGenInterface::readTetGenMesh failed to open %s", nodes_fname.c_str());
		if (!ins_ele)
			ERR ("TetGenInterface::readTetGenMesh failed to open %s", ele_fname.c_str());
		return nullptr;
	}

	std::vector<MeshLib::Node*> nodes;
	if (!readNodesFromStream (ins_nodes, nodes)) {
		// remove nodes read until now
		for (std::size_t k(0); k<nodes.size(); k++) {
			delete nodes[k];
		}
		return nullptr;
	}

	std::vector<MeshLib::Element*> elements;
	if (!readElementsFromStream (ins_ele, elements, nodes)) {
		// remove elements read until now
		for (std::size_t k(0); k<elements.size(); k++) {
			delete elements[k];
		}
		// remove nodes
		for (std::size_t k(0); k<nodes.size(); k++) {
			delete nodes[k];
		}
		return nullptr;
	}

	const std::string mesh_name (BaseLib::extractBaseNameWithoutExtension(nodes_fname));
	return new MeshLib::Mesh(mesh_name, nodes, elements);
}

bool TetGenInterface::readNodesFromStream (std::ifstream &ins, 
                                           std::vector<MeshLib::Node*> &nodes)
{
	std::string line;
	getline (ins, line);
	size_t pos_beg (line.find_first_not_of(" "));
	size_t n_nodes, dim, n_attributes;
	bool boundary_markers;

	while (!ins.fail())
	{
		line = line.substr(pos_beg);
		if (line.compare(0,1,"#") == 0)
		{
			// this line is a comment - skip
			getline (ins, line);
			pos_beg = line.find_first_not_of(" ");
			continue;
		}
		// read header line
		bool header_okay = parseNodesFileHeader(line, n_nodes, dim, n_attributes, boundary_markers);
		if (!header_okay)
			return false;
		if (!parseNodes(ins, nodes, n_nodes, dim))
			return false;
		return true;
	}
	return false;	
}

bool TetGenInterface::parseNodesFileHeader(std::string &line, 
                                           size_t &n_nodes, 
                                           size_t &dim,
                                           size_t &n_attributes, 
                                           bool &boundary_markers) const
{
	size_t pos_beg, pos_end;

	// number of nodes
	pos_beg = line.find_first_not_of (" ");
	pos_end = line.find_first_of(" ", pos_beg);
	if (pos_beg != std::string::npos && pos_end != std::string::npos)
		n_nodes = BaseLib::str2number<size_t> (line.substr(pos_beg, pos_end - pos_beg));
	else
	{
		ERR("TetGenInterface::parseNodesFileHeader(): could not number of nodes specified in header.");
		return false;
	}
	// dimension
	pos_beg = line.find_first_not_of (" ", pos_end);
	pos_end = line.find_first_of(" ", pos_beg);
	dim = BaseLib::str2number<size_t> (line.substr(pos_beg, pos_end - pos_beg));
	// number of attributes
	pos_beg = line.find_first_not_of (" ", pos_end);
	pos_end = line.find_first_of(" ", pos_beg);
	n_attributes = BaseLib::str2number<size_t> (line.substr(pos_beg, pos_end - pos_beg));
	// boundary marker at nodes?
	pos_beg = line.find_first_not_of (" ", pos_end);
	pos_end = line.find_first_of(" ", pos_beg);
	if (pos_end == std::string::npos)
		pos_end = line.size();
	if ((line.substr(pos_beg, pos_end - pos_beg)).compare("1") == 0)
		boundary_markers = true;
	else
		boundary_markers = false;

	return true;
}

bool TetGenInterface::parseNodes(std::ifstream &ins, 
                                 std::vector<MeshLib::Node*> &nodes, 
                                 size_t n_nodes, 
                                 size_t dim)
{
	std::size_t pos_beg, pos_end;
	std::string line;
	double* coordinates (static_cast<double*> (alloca (sizeof(double) * dim)));
	nodes.reserve(n_nodes);

	for (std::size_t k(0); k < n_nodes && !ins.fail(); k++) {
		getline(ins, line);
		if (ins.fail()) 
		{
			ERR("TetGenInterface::parseNodes(): Error reading node %d.", k);
			return false;
		}
		if (line.empty())
			continue;

		pos_end = 0;
		// read id
		size_t id;
		pos_beg = line.find_first_not_of(" ", pos_end);
		pos_end = line.find_first_of(" \n", pos_beg);
		if (pos_beg != std::string::npos && pos_end != std::string::npos) {
			id = BaseLib::str2number<size_t> (line.substr(pos_beg, pos_end - pos_beg));
			if (k == 0 && id == 0)
				_zero_based_idx = true;
		} else {
			ERR("TetGenInterface::parseNodes(): Error reading ID of node %d.", k);
			return false;
		}
		// read coordinates
		const unsigned offset = (_zero_based_idx) ? 0 : 1;
		for (size_t i(0); i < dim; i++) {
			pos_beg = line.find_first_not_of(" ", pos_end);
			pos_end = line.find_first_of(" \n", pos_beg);
			if (pos_end == std::string::npos) pos_end = line.size();
			if (pos_beg != std::string::npos)
				coordinates[i] = BaseLib::str2number<double>(line.substr(pos_beg, pos_end-pos_beg));
			else {
				ERR("TetGenInterface::parseNodes(): error reading coordinate %d of node %d.", i, k);
				return false;
			}
		}

		nodes.push_back(new MeshLib::Node(coordinates, id-offset));
		// read attributes and boundary markers ... - at the moment we do not use this information
	}

	return true;
}

bool TetGenInterface::readElementsFromStream(std::ifstream &ins, 
                                             std::vector<MeshLib::Element*> &elements, 
                                             const std::vector<MeshLib::Node*> &nodes)
{
	std::string line;
	getline (ins, line);
	size_t pos_beg (line.find_first_not_of(" "));
	size_t n_tets, n_nodes_per_tet;
	bool region_attributes;

	while (!ins.fail())
	{
		line = line.substr(pos_beg);
		if (line.compare(0,1,"#") == 0)
		{
			// this line is a comment - skip
			getline (ins, line);
			pos_beg = line.find_first_not_of(" ");
			continue;
		}
		
		// read header line
		bool header_okay = parseElementsFileHeader(line, n_tets, n_nodes_per_tet, region_attributes);
		if (!header_okay)
			return false;
		if (!parseElements(ins, elements, nodes, n_tets, n_nodes_per_tet, region_attributes))
			return false;
		return true;
	}
	return false;
}

bool TetGenInterface::parseElementsFileHeader(std::string &line,
                                              size_t& n_tets,
                                              size_t& n_nodes_per_tet,
                                              bool& region_attribute) const
{
	size_t pos_beg, pos_end;

	// number of tetrahedras
	pos_beg = line.find_first_not_of (" ");
	pos_end = line.find_first_of(" ", pos_beg);
	if (pos_beg != std::string::npos && pos_end != std::string::npos)
		n_tets = BaseLib::str2number<size_t> (line.substr(pos_beg, pos_end - pos_beg));
	else {
		ERR("TetGenInterface::parseElementsFileHeader(): Could not read number of tetrahedra specified in header.");
		return false;
	}
	// nodes per tet - either 4 or 10
	pos_beg = line.find_first_not_of (" \t", pos_end);
	pos_end = line.find_first_of(" \t", pos_beg);
	n_nodes_per_tet = BaseLib::str2number<size_t> (line.substr(pos_beg, pos_end - pos_beg));
	// region attribute at tetrahedra?
	pos_beg = line.find_first_not_of (" \t", pos_end);
	pos_end = line.find_first_of(" \t\n", pos_beg);
	if (pos_end == std::string::npos)
		pos_end = line.size();
	if ((line.substr(pos_beg, pos_end - pos_beg)).compare("1") == 0)
		region_attribute = true;
	else
		region_attribute = false;

	return true;
}

bool TetGenInterface::parseElements(std::ifstream& ins, 
                                    std::vector<MeshLib::Element*> &elements, 
                                    const std::vector<MeshLib::Node*> &nodes, 
                                    size_t n_tets, 
                                    size_t n_nodes_per_tet,
                                    bool region_attribute)
{
	size_t pos_beg, pos_end;
	std::string line;
	size_t* ids (static_cast<size_t*>(alloca (sizeof (size_t) * n_nodes_per_tet)));
	elements.reserve(n_tets);

	const unsigned offset = (_zero_based_idx) ? 0 : 1;
	for (size_t k(0); k < n_tets && !ins.fail(); k++)
	{
		getline (ins, line);
		if (ins.fail())
		{
			ERR("TetGenInterface::parseElements(): Error reading node %d.", k);
			return false;
		}
		if (line.empty())
			continue;
				
		pos_end = 0;
		// read id
		size_t id;
		pos_beg = line.find_first_not_of(" ", pos_end);
		pos_end = line.find_first_of(" \n", pos_beg);
		if (pos_beg != std::string::npos && pos_end != std::string::npos)
			id = BaseLib::str2number<size_t>(line.substr(pos_beg, pos_end - pos_beg));
		else {
			ERR("TetGenInterface::parseElements(): Error reading id of tetrahedron %d.", k);
			return false;
		}
		// read node ids
		for (size_t i(0); i < n_nodes_per_tet; i++)
		{
			pos_beg = line.find_first_not_of(" ", pos_end);
			pos_end = line.find_first_of(" ", pos_beg);
			if (pos_end == std::string::npos)
				pos_end = line.size();
			if (pos_beg != std::string::npos && pos_end != std::string::npos)
				ids[i] = BaseLib::str2number<std::size_t>(line.substr(pos_beg, pos_end - pos_beg)) - offset;
			else
			{
				ERR("TetGenInterface::parseElements(): Error reading node %d of tetrahedron %d.", i, k);
				return false;
			}
		}

		// read region attribute - this is something like material group
		unsigned region (0);
		if (region_attribute) {
			pos_beg = line.find_first_not_of(" ", pos_end);
			pos_end = line.find_first_of(" ", pos_beg);
			if (pos_end == std::string::npos) pos_end = line.size();
			if (pos_beg != std::string::npos && pos_end != std::string::npos)
				region = BaseLib::str2number<unsigned> (line.substr(pos_beg, pos_end - pos_beg));
			else {
				ERR("TetGenInterface::parseElements(): Error reading region attribute of tetrahedron %d.", k);
				return false;
			}
		}
		// insert new element into vector
		MeshLib::Node** tet_nodes = new MeshLib::Node*[4];
		for (unsigned k(0); k<4; k++) {
			tet_nodes[k] = nodes[ids[k]];
		}
		elements.push_back (new MeshLib::Tet(tet_nodes, region));
	}
	return true;
}

} // end namespace FileIO
