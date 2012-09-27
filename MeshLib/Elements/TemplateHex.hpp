/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file TemplateHex.hpp
 *
 * Created on 2012-05-02 by Karsten Rink
 */

#include "Node.h"
#include "Quad.h"
#include "Prism.h"

#include "MathTools.h"

namespace MeshLib {

template <unsigned ORDER, unsigned NNODES>
const unsigned TemplateHex<ORDER,NNODES>::_face_nodes[6][4] =
{
	{0, 3, 2, 1}, // Face 0
	{0, 1, 5, 4}, // Face 1
	{1, 2, 6, 5}, // Face 2
	{2, 3, 7, 6}, // Face 3
	{3, 0, 4, 7}, // Face 4
	{4, 5, 6, 7}  // Face 5
};

template <unsigned ORDER, unsigned NNODES>
const unsigned TemplateHex<ORDER,NNODES>::_edge_nodes[12][2] =
{
	{0, 1}, // Edge 0
	{1, 2}, // Edge 1
	{2, 3}, // Edge 2
	{0, 3}, // Edge 3
	{0, 4}, // Edge 4
	{1, 5}, // Edge 5
	{2, 6}, // Edge 6
	{3, 7}, // Edge 7
	{4, 5}, // Edge 8
	{5, 6}, // Edge 9
	{6, 7}, // Edge 10
	{4, 7}  // Edge 11
};

template <unsigned ORDER, unsigned NNODES>
TemplateHex<ORDER,NNODES>::TemplateHex(Node* nodes[8], unsigned value)
	: Cell(value)
{
	_nodes = nodes;

	_neighbors = new Element*[6];
	for (unsigned i=0; i<6; i++)
		_neighbors[i] = NULL;

	this->_volume = this->computeVolume();
}

template <unsigned ORDER, unsigned NNODES>
TemplateHex<ORDER,NNODES>::TemplateHex(const TemplateHex<ORDER,NNODES> &hex)
	: Cell(hex.getValue())
{
	_nodes = new Node*[NNODES];
	for (unsigned i=0; i<NNODES; i++)
		_nodes[i] = hex._nodes[i];

	_neighbors = new Element*[6];
	for (unsigned i=0; i<6; i++)
		_neighbors[i] = hex._neighbors[i];

	_volume = hex.getVolume();
}

template <unsigned ORDER, unsigned NNODES>
TemplateHex<ORDER,NNODES>::~TemplateHex()
{
}

template <unsigned ORDER, unsigned NNODES>
double TemplateHex<ORDER,NNODES>::computeVolume()
{
	return MathLib::calcTetrahedronVolume(_nodes[4]->getCoords(), _nodes[7]->getCoords(), _nodes[5]->getCoords(), _nodes[0]->getCoords())
		 + MathLib::calcTetrahedronVolume(_nodes[5]->getCoords(), _nodes[3]->getCoords(), _nodes[1]->getCoords(), _nodes[0]->getCoords())
		 + MathLib::calcTetrahedronVolume(_nodes[5]->getCoords(), _nodes[7]->getCoords(), _nodes[3]->getCoords(), _nodes[0]->getCoords())
		 + MathLib::calcTetrahedronVolume(_nodes[5]->getCoords(), _nodes[7]->getCoords(), _nodes[6]->getCoords(), _nodes[2]->getCoords())
		 + MathLib::calcTetrahedronVolume(_nodes[1]->getCoords(), _nodes[3]->getCoords(), _nodes[5]->getCoords(), _nodes[2]->getCoords())
		 + MathLib::calcTetrahedronVolume(_nodes[3]->getCoords(), _nodes[7]->getCoords(), _nodes[5]->getCoords(), _nodes[2]->getCoords());
}

template <unsigned ORDER, unsigned NNODES>
const Element* TemplateHex<ORDER,NNODES>::getFace(unsigned i) const
{
	if (i<this->getNFaces())
	{
		unsigned nFaceNodes (this->getNFaceNodes(i));
		Node** nodes = new Node*[nFaceNodes];
		for (unsigned j=0; j<nFaceNodes; j++)
			nodes[j] = _nodes[_face_nodes[i][j]];
		return new Quad(nodes);
	}
	std::cerr << "Error in MeshLib::Element::getFace() - Index does not exist." << std::endl;
	return NULL;
}

template <unsigned ORDER, unsigned NNODES>
bool TemplateHex<ORDER,NNODES>::isEdge(unsigned idx1, unsigned idx2) const
{
	for (unsigned i(0); i<12; i++)
	{
		if (_edge_nodes[i][0]==idx1 && _edge_nodes[i][1]==idx2) return true;
		if (_edge_nodes[i][1]==idx1 && _edge_nodes[i][0]==idx2) return true;
	}
	return false;
}

template <unsigned ORDER, unsigned NNODES>
Element* TemplateHex<ORDER,NNODES>::clone() const
{
	return new TemplateHex<ORDER,NNODES>(*this);
}

template <unsigned ORDER, unsigned NNODES>
unsigned TemplateHex<ORDER,NNODES>::identifyFace(Node* nodes[3]) const
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

template <unsigned ORDER, unsigned NNODES>
Element* TemplateHex<ORDER,NNODES>::reviseElement() const
{
	std::vector<size_t> collapsed_edges;
	for (size_t edge(0); edge<getNEdges(); edge++) {
		if (_nodes[_edge_nodes[edge][0]] == _nodes[_edge_nodes[edge][1]]) {
			collapsed_edges.push_back(edge);
		}
	}

	if (collapsed_edges.size() == 1) {
		std::cerr << "[TemplateHex<ORDER,NNODES>::reviseElement()] collapsing of one edge in hexahedron not handled" << std::endl;
		return NULL;
	}

	if (collapsed_edges.size() == 2) {
		// try to create a prism out of the hex
		if (collapsed_edges[0] == 0 && collapsed_edges[1] == 2) {
			Node** prism_nodes(new Node*[6]);
			prism_nodes[0] = _nodes[0];
			prism_nodes[1] = _nodes[4];
			prism_nodes[2] = _nodes[5];
			prism_nodes[3] = _nodes[3];
			prism_nodes[4] = _nodes[7];
			prism_nodes[5] = _nodes[6];
			return new Prism(prism_nodes, _value);
		}
		if (collapsed_edges[0] == 1 && collapsed_edges[1] == 3) {
			Node** prism_nodes(new Node*[6]);
			prism_nodes[0] = _nodes[0];
			prism_nodes[1] = _nodes[4];
			prism_nodes[2] = _nodes[7];
			prism_nodes[3] = _nodes[1];
			prism_nodes[4] = _nodes[5];
			prism_nodes[5] = _nodes[6];
			return new Prism(prism_nodes, _value);
		}
		if (collapsed_edges[0] == 4 && collapsed_edges[1] == 5) {
			Node** prism_nodes(new Node*[6]);
			prism_nodes[0] = _nodes[0];
			prism_nodes[1] = _nodes[7];
			prism_nodes[2] = _nodes[3];
			prism_nodes[3] = _nodes[1];
			prism_nodes[4] = _nodes[6];
			prism_nodes[5] = _nodes[2];
			return new Prism(prism_nodes, _value);
		}
		if (collapsed_edges[0] == 5 && collapsed_edges[1] == 6) {
			Node** prism_nodes(new Node*[6]);
			prism_nodes[0] = _nodes[0];
			prism_nodes[1] = _nodes[1];
			prism_nodes[2] = _nodes[4];
			prism_nodes[3] = _nodes[3];
			prism_nodes[4] = _nodes[2];
			prism_nodes[5] = _nodes[7];
			return new Prism(prism_nodes, _value);
		}
		if (collapsed_edges[0] == 6 && collapsed_edges[1] == 7) {
			Node** prism_nodes(new Node*[6]);
			prism_nodes[0] = _nodes[0];
			prism_nodes[1] = _nodes[3];
			prism_nodes[2] = _nodes[4];
			prism_nodes[3] = _nodes[1];
			prism_nodes[4] = _nodes[2];
			prism_nodes[5] = _nodes[5];
			return new Prism(prism_nodes, _value);
		}
		if (collapsed_edges[0] == 7 && collapsed_edges[1] == 4) {
			Node** prism_nodes(new Node*[6]);
			prism_nodes[0] = _nodes[0];
			prism_nodes[1] = _nodes[1];
			prism_nodes[2] = _nodes[5];
			prism_nodes[3] = _nodes[3];
			prism_nodes[4] = _nodes[2];
			prism_nodes[5] = _nodes[6];
			return new Prism(prism_nodes, _value);
		}
		if (collapsed_edges[0] == 8 && collapsed_edges[1] == 10) {
			Node** prism_nodes(new Node*[6]);
			prism_nodes[0] = _nodes[0];
			prism_nodes[1] = _nodes[1];
			prism_nodes[2] = _nodes[4];
			prism_nodes[3] = _nodes[3];
			prism_nodes[4] = _nodes[2];
			prism_nodes[5] = _nodes[7];
			return new Prism(prism_nodes, _value);
		}
		if (collapsed_edges[0] == 9 && collapsed_edges[1] == 11) {
			Node** prism_nodes(new Node*[6]);
			prism_nodes[0] = _nodes[0];
			prism_nodes[1] = _nodes[3];
			prism_nodes[2] = _nodes[4];
			prism_nodes[3] = _nodes[1];
			prism_nodes[4] = _nodes[2];
			prism_nodes[5] = _nodes[5];
			return new Prism(prism_nodes, _value);
		}
		return NULL;
	}

	return NULL;
}

} // end namespace MeshLib
