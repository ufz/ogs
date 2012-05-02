/**
 * Mesh.cpp
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#include "Mesh.h"

#include "Node.h"
#include "Elements/Tri.h"
#include "Elements/Quad.h"
#include "Elements/Tet.h"
#include "Elements/Hex.h"
#include "Elements/Pyramid.h"
#include "Elements/Prism.h"

namespace MeshLib {

Mesh::Mesh(const std::string &name, const std::vector<Node*> &nodes, const std::vector<Element*> &elements)
	: _name(name), _nodes(nodes), _elements(elements)
{
	this->removeIdenticalNodes();
}

Mesh::Mesh(const Mesh &mesh)
	: _name(mesh.getName()), _nodes(mesh.getNodes()), _elements(mesh.getElements())
{
}

Mesh::Mesh(const std::string &file_name)
{
	// read mesh
	this->removeIdenticalNodes();
}

Mesh::~Mesh()
{
	const size_t nElements (_elements.size());
	for (size_t i=0; i<nElements; i++)
		delete _elements[i];

	const size_t nNodes (_nodes.size());
	for (size_t i=0; i<nNodes; i++)
		delete _nodes[i];
}

void Mesh::removeIdenticalNodes()
{
	//check for unique mesh nodes
}

}