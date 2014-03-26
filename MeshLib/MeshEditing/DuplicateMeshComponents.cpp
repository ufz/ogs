/**
 * \file   DuplicateMeshComponents.cpp
 * \author Karsten Rink
 * \date   2014-03-25
 * \brief  Implementation of Duplicate functions.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DuplicateMeshComponents.h"

#include "Mesh.h"
#include "Node.h"
#include "Elements/Element.h"
#include "Elements/Line.h"
#include "Elements/Tri.h"
#include "Elements/Quad.h"
#include "Elements/Tet.h"
#include "Elements/Hex.h"
#include "Elements/Pyramid.h"
#include "Elements/Prism.h"

namespace MeshLib
{

std::vector<MeshLib::Node*> copyNodeVector(const std::vector<MeshLib::Node*> &nodes)
{
	const std::size_t nNodes(nodes.size());
	std::vector<MeshLib::Node*> new_nodes;
	new_nodes.reserve(nNodes);
	for (std::size_t k = 0; k < nNodes; ++k)
		new_nodes.push_back(new MeshLib::Node(nodes[k]->getCoords(), new_nodes.size()));
	return new_nodes;
}

std::vector<MeshLib::Element*> copyElementVector(const std::vector<MeshLib::Element*> &elements, const std::vector<MeshLib::Node*> &nodes)
{
	const std::size_t nElements(elements.size());
	std::vector<MeshLib::Element*> new_elements;
	new_elements.reserve(nElements);
	for (std::size_t k = 0; k < nElements; ++k)
		new_elements.push_back(copyElement(elements[k], nodes));
	return new_elements;
}

MeshLib::Element* copyElement(MeshLib::Element const*const element, const std::vector<MeshLib::Node*> &nodes)
{
	if (element->getGeomType() == MeshElemType::LINE)
		return copyLine(element, nodes);
	else if (element->getGeomType() == MeshElemType::TRIANGLE)
		return copyTri(element, nodes);
	else if (element->getGeomType() == MeshElemType::QUAD)
		return copyQuad(element, nodes);
	else if (element->getGeomType() == MeshElemType::TETRAHEDRON)
		return copyTet(element, nodes);
	else if (element->getGeomType() == MeshElemType::HEXAHEDRON)
		return copyHex(element, nodes);
	else if (element->getGeomType() == MeshElemType::PYRAMID)
		return copyPyramid(element, nodes);
	else if (element->getGeomType() == MeshElemType::PRISM)
		return copyPrism(element, nodes);

	ERR ("Error: Unknown element type.");
	return nullptr;
}

MeshLib::Element* copyLine(MeshLib::Element const*const org_elem, const std::vector<MeshLib::Node*> &nodes)
{
	MeshLib::Node** new_nodes = new MeshLib::Node*[2];
	new_nodes[0] = nodes[org_elem->getNode(0)->getID()];
	new_nodes[1] = nodes[org_elem->getNode(1)->getID()];
	return new MeshLib::Line(new_nodes, org_elem->getValue());
}

MeshLib::Element* copyTri(MeshLib::Element const*const org_elem, const std::vector<MeshLib::Node*> &nodes)
{
	MeshLib::Node** new_nodes = new MeshLib::Node*[3];
	for (unsigned i=0; i<3; ++i)
		new_nodes[i] = nodes[org_elem->getNode(i)->getID()];
	return new MeshLib::Tri(new_nodes, org_elem->getValue());
}

MeshLib::Element* copyQuad(MeshLib::Element const*const org_elem, const std::vector<MeshLib::Node*> &nodes)
{
	MeshLib::Node** new_nodes = new MeshLib::Node*[4];
	for (unsigned i=0; i<4; ++i)
		new_nodes[i] = nodes[org_elem->getNode(i)->getID()];
	return new MeshLib::Quad(new_nodes, org_elem->getValue());
}

MeshLib::Element* copyTet(MeshLib::Element const*const org_elem, const std::vector<MeshLib::Node*> &nodes)
{
	MeshLib::Node** new_nodes = new MeshLib::Node*[4];
	for (unsigned i=0; i<4; ++i)
		new_nodes[i] = nodes[org_elem->getNode(i)->getID()];
	return new MeshLib::Tet(new_nodes, org_elem->getValue());
}

MeshLib::Element* copyHex(MeshLib::Element const*const org_elem, const std::vector<MeshLib::Node*> &nodes)
{
	MeshLib::Node** new_nodes = new MeshLib::Node*[8];
	for (unsigned i=0; i<8; ++i)
		new_nodes[i] = nodes[org_elem->getNode(i)->getID()];
	return new MeshLib::Hex(new_nodes, org_elem->getValue());
}

MeshLib::Element* copyPyramid(MeshLib::Element const*const org_elem, const std::vector<MeshLib::Node*> &nodes)
{
	MeshLib::Node** new_nodes = new MeshLib::Node*[5];
	for (unsigned i=0; i<5; ++i)
		new_nodes[i] = nodes[org_elem->getNode(i)->getID()];
	return new MeshLib::Pyramid(new_nodes, org_elem->getValue());
}

MeshLib::Element* copyPrism(MeshLib::Element const*const org_elem, const std::vector<MeshLib::Node*> &nodes)
{
	MeshLib::Node** new_nodes = new MeshLib::Node*[6];
	for (unsigned i=0; i<6; ++i)
		new_nodes[i] = nodes[org_elem->getNode(i)->getID()];
	return new MeshLib::Prism(new_nodes, org_elem->getValue());
}

} // namespace MeshLib
