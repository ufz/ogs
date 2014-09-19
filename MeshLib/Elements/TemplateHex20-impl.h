/**
 * \copyright
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "logog/include/logog.hpp"

#include "Node.h"
#include "Quad.h"
#include "Prism.h"

#include "AnalyticalGeometry.h"

namespace MeshLib {

template <unsigned NNODES, CellType CELLHEXTYPE>
const unsigned TemplateHex20<NNODES, CELLHEXTYPE>::n_all_nodes;

template <unsigned NNODES, CellType CELLHEXTYPE>
const unsigned TemplateHex20<NNODES, CELLHEXTYPE>::n_base_nodes;

template <unsigned NNODES, CellType CELLHEXTYPE>
const unsigned TemplateHex20<NNODES,CELLHEXTYPE>::_face_nodes[6][8] =
{
	{0, 3, 2, 1, 11, 10,  9,  8}, // Face 0
	{0, 1, 5, 4, 8,  17, 12, 16}, // Face 1
	{1, 2, 6, 5, 9,  18, 13, 17}, // Face 2
	{2, 3, 7, 6, 10, 19, 14, 18}, // Face 3
	{3, 0, 4, 7, 11, 16, 15, 12}, // Face 4
	{4, 5, 6, 7, 12, 13, 14, 15}  // Face 5
};

template <unsigned NNODES, CellType CELLHEXTYPE>
const unsigned TemplateHex20<NNODES,CELLHEXTYPE>::_edge_nodes[12][3] =
{
	{0, 1, 8}, // Edge 0
	{1, 2, 9}, // Edge 1
	{2, 3, 10}, // Edge 2
	{0, 3, 11}, // Edge 3
	{4, 5, 12}, // Edge 4
	{5, 6, 13}, // Edge 5
	{6, 7, 14}, // Edge 6
	{4, 7, 15}, // Edge 7
	{0, 4, 16}, // Edge 8
	{1, 5, 17}, // Edge 9
	{2, 6, 18}, // Edge 10
	{3, 7, 19}  // Edge 11
};

template <unsigned NNODES, CellType CELLHEXTYPE>
TemplateHex20<NNODES,CELLHEXTYPE>::TemplateHex20(Node* nodes[NNODES], unsigned value, std::size_t id)
	: Cell(value, id)
{
	_nodes = nodes;

	_neighbors = new Element*[6];
	std::fill(_neighbors, _neighbors + 6, nullptr);

	this->_volume = this->computeVolume();
}

template<unsigned NNODES, CellType CELLHEXTYPE>
TemplateHex20<NNODES,CELLHEXTYPE>::TemplateHex20(std::array<Node*, NNODES> const& nodes,
                                             unsigned value, std::size_t id)
	: Cell(value, id)
{
	_nodes = new Node*[NNODES];
	std::copy(nodes.begin(), nodes.end(), _nodes);

	_neighbors = new Element*[6];
	std::fill(_neighbors, _neighbors + 6, nullptr);

	this->_volume = this->computeVolume();
}

template <unsigned NNODES, CellType CELLHEXTYPE>
TemplateHex20<NNODES,CELLHEXTYPE>::TemplateHex20(const TemplateHex20<NNODES,CELLHEXTYPE> &hex)
	: Cell(hex.getValue(), hex.getID())
{
	_nodes = new Node*[NNODES];
	for (unsigned i=0; i<NNODES; i++)
		_nodes[i] = hex._nodes[i];

	_neighbors = new Element*[6];
	for (unsigned i=0; i<6; i++)
		_neighbors[i] = hex._neighbors[i];

	_volume = hex.getVolume();
}

template <unsigned NNODES, CellType CELLHEXTYPE>
TemplateHex20<NNODES,CELLHEXTYPE>::~TemplateHex20()
{
}

template <unsigned NNODES, CellType CELLHEXTYPE>
double TemplateHex20<NNODES,CELLHEXTYPE>::computeVolume()
{
	return GeoLib::calcTetrahedronVolume(_nodes[4]->getCoords(), _nodes[7]->getCoords(), _nodes[5]->getCoords(), _nodes[0]->getCoords())
		 + GeoLib::calcTetrahedronVolume(_nodes[5]->getCoords(), _nodes[3]->getCoords(), _nodes[1]->getCoords(), _nodes[0]->getCoords())
		 + GeoLib::calcTetrahedronVolume(_nodes[5]->getCoords(), _nodes[7]->getCoords(), _nodes[3]->getCoords(), _nodes[0]->getCoords())
		 + GeoLib::calcTetrahedronVolume(_nodes[5]->getCoords(), _nodes[7]->getCoords(), _nodes[6]->getCoords(), _nodes[2]->getCoords())
		 + GeoLib::calcTetrahedronVolume(_nodes[1]->getCoords(), _nodes[3]->getCoords(), _nodes[5]->getCoords(), _nodes[2]->getCoords())
		 + GeoLib::calcTetrahedronVolume(_nodes[3]->getCoords(), _nodes[7]->getCoords(), _nodes[5]->getCoords(), _nodes[2]->getCoords());
}

template <unsigned NNODES, CellType CELLHEXTYPE>
const Element* TemplateHex20<NNODES,CELLHEXTYPE>::getFace(unsigned i) const
{
	if (i<this->getNFaces())
	{
		unsigned nFaceNodes (this->getNFaceNodes(i));
		Node** nodes = new Node*[nFaceNodes];
		for (unsigned j=0; j<nFaceNodes; j++)
			nodes[j] = _nodes[_face_nodes[i][j]];
		return new Quad(nodes);
	}
	ERR("Error in MeshLib::Element::getFace() - Index %d does not exist.", i);
	return NULL;
}

template <unsigned NNODES, CellType CELLHEXTYPE>
bool TemplateHex20<NNODES,CELLHEXTYPE>::isEdge(unsigned idx1, unsigned idx2) const
{
	for (unsigned i(0); i<12; i++)
	{
		if (_edge_nodes[i][0]==idx1 && _edge_nodes[i][1]==idx2) return true;
		if (_edge_nodes[i][1]==idx1 && _edge_nodes[i][0]==idx2) return true;
	}
	return false;
}

template <unsigned NNODES, CellType CELLHEXTYPE>
Element* TemplateHex20<NNODES,CELLHEXTYPE>::clone() const
{
	return new TemplateHex20<NNODES,CELLHEXTYPE>(*this);
}

template <unsigned NNODES, CellType CELLHEXTYPE>
unsigned TemplateHex20<NNODES,CELLHEXTYPE>::identifyFace(Node* nodes[3]) const
{
	for (unsigned i=0; i<6; i++)
	{
		unsigned flag(0);
		for (unsigned j=0; j<4; j++)
			for (unsigned k=0; k<3; k++)
				if (_nodes[_face_nodes[i][j]] == nodes[k])
					flag++;
		if (flag==3)
			return i;
	}
	return std::numeric_limits<unsigned>::max();
}

template <unsigned NNODES, CellType CELLHEXTYPE>
ElementErrorCode TemplateHex20<NNODES,CELLHEXTYPE>::validate() const
{
	ElementErrorCode error_code;
	error_code[ElementErrorFlag::ZeroVolume] = this->hasZeroVolume();
		
	for (unsigned i=0; i<6; ++i)
	{
		if (error_code.all())
			break;

		const MeshLib::Element* quad (this->getFace(i));
		error_code |= quad->validate();
		delete quad;
	}
	error_code[ElementErrorFlag::NodeOrder]  = !this->testElementNodeOrder();
	return error_code;
}

} // end namespace MeshLib
