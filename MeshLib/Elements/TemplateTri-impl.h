/**
 * \file
 * \author Thomas Fischer
 * \date   Sep 27, 2012
 * \brief  Implementation of the TemplateTri class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

namespace MeshLib {

template <unsigned NNODES, CellType CELLTRITYPE, typename EDGENODES>
const unsigned TemplateTri<NNODES, CELLTRITYPE, EDGENODES>::n_all_nodes;

template <unsigned NNODES, CellType CELLTRITYPE, typename EDGENODES>
const unsigned TemplateTri<NNODES, CELLTRITYPE, EDGENODES>::n_base_nodes;

template <unsigned NNODES, CellType CELLTRITYPE, typename EDGENODES>
TemplateTri<NNODES,CELLTRITYPE,EDGENODES>::TemplateTri(Node* nodes[NNODES], unsigned value, std::size_t id) :
	Face(value, id)
{
	_nodes = nodes;
	_neighbors = new Element*[3];
	std::fill(_neighbors, _neighbors + 3, nullptr);
	this->_area = this->computeVolume();
}

template <unsigned NNODES, CellType CELLTRITYPE, typename EDGENODES>
TemplateTri<NNODES,CELLTRITYPE,EDGENODES>::TemplateTri(std::array<Node*, NNODES> const& nodes,
                                             unsigned value, std::size_t id)
	: Face(value, id)
{
	_nodes = new Node*[NNODES];
	std::copy(nodes.begin(), nodes.end(), _nodes);

	_neighbors = new Element*[3];
	std::fill(_neighbors, _neighbors + 3, nullptr);

	this->_area = this->computeVolume();
}

template <unsigned NNODES, CellType CELLTRITYPE, typename EDGENODES>
TemplateTri<NNODES,CELLTRITYPE,EDGENODES>::TemplateTri(const TemplateTri<NNODES,CELLTRITYPE,EDGENODES> &tri) :
	Face(tri.getValue(), tri.getID())
{
	_nodes = new Node*[NNODES];
	for (unsigned i=0; i<NNODES; i++)
	{
		_nodes[i] = tri._nodes[i];
	}

	_neighbors = new Element*[3];
	for (unsigned i=0; i<3; i++) {
		_neighbors[i] = tri._neighbors[i];
	}

	_area = tri.getArea();
}

template <unsigned NNODES, CellType CELLTRITYPE, typename EDGENODES>
TemplateTri<NNODES,CELLTRITYPE,EDGENODES>::~TemplateTri()
{}

template <unsigned NNODES, CellType CELLTRITYPE, typename EDGENODES>
bool TemplateTri<NNODES,CELLTRITYPE,EDGENODES>::isEdge(unsigned idx1, unsigned idx2) const
{
	for (unsigned i(0); i<3; i++)
	{
		if (_edge_nodes[i][0]==idx1 && _edge_nodes[i][1]==idx2) return true;
		if (_edge_nodes[i][1]==idx1 && _edge_nodes[i][0]==idx2) return true;
	}
	return false;
}

template <unsigned NNODES, CellType CELLTRITYPE, typename EDGENODES>
bool TemplateTri<NNODES,CELLTRITYPE,EDGENODES>::isPntInElement(GeoLib::Point const& pnt, double eps) const
{
	return GeoLib::isPointInTriangle(pnt, *_nodes[0], *_nodes[1], *_nodes[2], eps);
}

template <unsigned NNODES, CellType CELLTRITYPE, typename EDGENODES>
ElementErrorCode TemplateTri<NNODES,CELLTRITYPE,EDGENODES>::validate() const
{ 
	ElementErrorCode error_code;
	error_code[ElementErrorFlag::ZeroVolume] = this->hasZeroVolume();
	error_code[ElementErrorFlag::NodeOrder]  = !this->testElementNodeOrder();
	return error_code;
}

template <unsigned NNODES, CellType CELLTRITYPE, typename EDGENODES>
unsigned TemplateTri<NNODES,CELLTRITYPE,EDGENODES>::identifyFace(Node* nodes[3]) const
{
	for (unsigned i=0; i<3; i++)
	{
		unsigned flag(0);
		for (unsigned j=0; j<2; j++)
			for (unsigned k=0; k<2; k++)
				if (_nodes[_edge_nodes[i][j]] == nodes[k])
					flag++;
		if (flag==2)
			return i;
	}
	return std::numeric_limits<unsigned>::max();
}

} // namespace MeshLib

