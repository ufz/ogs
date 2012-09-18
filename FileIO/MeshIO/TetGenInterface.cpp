/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file TetGenInterface.cpp
 *
 *  Created on 2011-09-12 by Thomas Fischer
 */

#include <cstddef>
#include <string>

// BaseLib
#include "StringTools.h"

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
	_nodes(), _elements(), _zero_based_idx (false)
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
			std::cout << "TetGenInterface::readTetGenMesh failed to open " <<
			nodes_fname << std::endl;
		if (!ins_ele)
			std::cout << "TetGenInterface::readTetGenMesh failed to open " <<
			ele_fname << std::endl;
		return NULL;
	}

	if (!readNodesFromStream (ins_nodes)) {
		// remove nodes read until now
		for (std::size_t k(0); k<_nodes.size(); k++) {
			delete _nodes[k];
		}
		return NULL;
	}

	if (!readElementsFromStream (ins_ele)) {
		// remove elements read until now
		for (std::size_t k(0); k<_elements.size(); k++) {
			delete _elements[k];
		}
		// remove nodes
		for (std::size_t k(0); k<_nodes.size(); k++) {
			delete _nodes[k];
		}
		return NULL;
	}

	return new MeshLib::Mesh(nodes_fname, _nodes, _elements);
}

bool TetGenInterface::readNodesFromStream (std::ifstream &ins)
{
	std::string line;
	getline (ins, line);
	size_t pos_beg (line.find_first_not_of(" "));
	size_t n_nodes, dim, n_attributes;
	bool boundary_markers;
	bool not_read_header (true);

	while (!ins.fail() && not_read_header)
	{
		line = line.substr(pos_beg);
		if (line.compare(0,1,"#") == 0)
		{
			// this line is a comment - skip
			getline (ins, line);
			pos_beg = line.find_first_not_of(" ");
		}
		else
			// read header line
			not_read_header = !parseNodesFileHeader(line,
			                                        n_nodes,
			                                        dim,
			                                        n_attributes,
			                                        boundary_markers);
	}
	if (not_read_header)
		return false;
	if (!parseNodes(ins, n_nodes, dim))
		return false;

	return true;
}

bool TetGenInterface::parseNodesFileHeader(std::string &line,
                                           size_t& n_nodes,
                                           size_t& dim,
                                           size_t& n_attributes,
                                           bool& boundary_markers) const
{
	size_t pos_beg, pos_end;

	// number of nodes
	pos_beg = line.find_first_not_of (" ");
	pos_end = line.find_first_of(" ", pos_beg);
	if (pos_beg != std::string::npos && pos_end != std::string::npos)
		n_nodes = str2number<size_t> (line.substr(pos_beg, pos_end - pos_beg));
	else
	{
		std::cout << "TetGenInterface::parseNodesFileHeader could not correct read TetGen mesh header - number of nodes"
		          << std::endl;
		return false;
	}
	// dimension
	pos_beg = line.find_first_not_of (" ", pos_end);
	pos_end = line.find_first_of(" ", pos_beg);
	dim = str2number<size_t> (line.substr(pos_beg, pos_end - pos_beg));
	// number of attributes
	pos_beg = line.find_first_not_of (" ", pos_end);
	pos_end = line.find_first_of(" ", pos_beg);
	n_attributes = str2number<size_t> (line.substr(pos_beg, pos_end - pos_beg));
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

bool TetGenInterface::parseNodes(std::ifstream& ins, size_t n_nodes, size_t dim)
{
	std::size_t pos_beg, pos_end;
	std::string line;
	double* coordinates (static_cast<double*> (alloca (sizeof(double) * dim)));

	for (std::size_t k(0); k < n_nodes && !ins.fail(); k++) {
		getline(ins, line);
		if (!ins.fail()) {
			if (!line.empty()) {
				pos_end = 0;
				// read id
				size_t id;
				pos_beg = line.find_first_not_of(" ", pos_end);
				pos_end = line.find_first_of(" \n", pos_beg);
				if (pos_beg != std::string::npos && pos_end != std::string::npos) {
					id = str2number<size_t> (line.substr(pos_beg, pos_end - pos_beg));
					if (k == 0 && id == 0)
						_zero_based_idx = true;
				} else {
					std::cout << "error reading id of node " << k
									<< " in TetGenInterface::parseNodes" << std::endl;
					return false;
				}
				// read coordinates
				for (size_t i(0); i < dim; i++) {
					pos_beg = line.find_first_not_of(" ", pos_end);
					pos_end = line.find_first_of(" \n", pos_beg);
					if (pos_end == std::string::npos) pos_end = line.size();
					if (pos_beg != std::string::npos)
						coordinates[i] = str2number<double> (
										line.substr(pos_beg, pos_end - pos_beg));
					else {
						std::cout << "error reading coordinate " << i << " of node " << k
										<< " in TetGenInterface::parseNodes" << std::endl;
						return false;
					}
				}
				if (!_zero_based_idx) id--;
				// since CFEMesh is our friend we can access private data of mesh
				_nodes.push_back(new MeshLib::Node(coordinates, id));
				// read attributes and boundary markers ... - at the moment we do not use this information
			}
		} else {
			std::cout << "error reading node " << k << " in TetGenInterface::parseNodes" << std::endl;
			return false;
		}
	}

	return true;
}

bool TetGenInterface::readElementsFromStream(std::ifstream &ins)
{
	std::string line;
	getline (ins, line);
	size_t pos_beg (line.find_first_not_of(" "));
	size_t n_tets, n_nodes_per_tet;
	bool region_attributes;
	bool not_read_header (true);

	while (!ins.fail() && not_read_header)
	{
		line = line.substr(pos_beg);
		if (line.compare(0,1,"#") == 0)
		{
			// this line is a comment - skip
			getline (ins, line);
			pos_beg = line.find_first_not_of(" ");
		}
		else
			// read header line
			not_read_header = !parseElementsFileHeader(line,
			                                           n_tets,
			                                           n_nodes_per_tet,
			                                           region_attributes);
	}
	if (not_read_header)
		return false;
	if (!parseElements(ins, n_tets, n_nodes_per_tet, region_attributes))
		return false;

	return true;
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
		n_tets = str2number<size_t> (line.substr(pos_beg, pos_end - pos_beg));
	else
	{
		std::cout <<
		"TetGenInterface::parseElementsFileHeader could not correct read TetGen mesh header - number of tetrahedras"
		          << std::endl;
		return false;
	}
	// nodes per tet - either 4 or 10
	pos_beg = line.find_first_not_of (" \t", pos_end);
	pos_end = line.find_first_of(" \t", pos_beg);
	n_nodes_per_tet = str2number<size_t> (line.substr(pos_beg, pos_end - pos_beg));
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

bool TetGenInterface::parseElements(std::ifstream& ins, size_t n_tets, size_t n_nodes_per_tet,
                                    bool region_attribute)
{
	size_t pos_beg, pos_end;
	std::string line;
	size_t* ids (static_cast<size_t*>(alloca (sizeof (size_t) * n_nodes_per_tet)));

	for (size_t k(0); k < n_tets && !ins.fail(); k++)
	{
		getline (ins, line);
		if (!ins.fail())
		{
			if (!line.empty())
			{
				pos_end = 0;
				// read id
				size_t id;
				pos_beg = line.find_first_not_of(" ", pos_end);
				pos_end = line.find_first_of(" \n", pos_beg);
				if (pos_beg != std::string::npos && pos_end != std::string::npos)
					id = str2number<size_t>(line.substr(pos_beg, pos_end - pos_beg));
				else {
					std::cout << "error reading id of tetrahedra " << k <<
					" in TetGenInterface::parseElements" << std::endl;
					return false;
				}
				// read node ids
				for (size_t i(0); i < n_nodes_per_tet; i++)
				{
					pos_beg = line.find_first_not_of(" ", pos_end);
					pos_end = line.find_first_of(" ", pos_beg);
					if (pos_end == std::string::npos)
						pos_end = line.size();
					if (pos_beg != std::string::npos && pos_end !=
					    std::string::npos)
						ids[i] = str2number<std::size_t>(line.substr(pos_beg, pos_end - pos_beg));
					else
					{
						std::cout << "error reading node " << i << " of tetrahedra " << k <<
							" in TetGenInterface::parseElements" << std::endl;
						return false;
					}
				}
				if (!_zero_based_idx) {
					id--;
					for (size_t i(0); i < n_nodes_per_tet; i++)
						ids[i]--;
				}

				// read region attribute - this is something like material group
				unsigned region (0);
				if (region_attribute) {
					pos_beg = line.find_first_not_of(" ", pos_end);
					pos_end = line.find_first_of(" ", pos_beg);
					if (pos_end == std::string::npos) pos_end = line.size();
					if (pos_beg != std::string::npos && pos_end != std::string::npos)
						region = str2number<unsigned> (line.substr(pos_beg, pos_end - pos_beg));
					else {
						std::cout << "error reading region attribute of tetrahedra " << k
										<< " in TetGenInterface::parseElements" << std::endl;
						return false;
					}
				}
				// insert new element into vector
				_elements.push_back (new MeshLib::Tet(_nodes[ids[0]], _nodes[ids[1]], _nodes[ids[2]], _nodes[ids[3]], region));

			}
		}
		else
		{
			std::cout << "error reading node " << k << " in TetGenInterface::parseElements" << std::endl;
			return false;
		}
	}
	return true;
}

} // end namespace FileIO
