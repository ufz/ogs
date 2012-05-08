/**
 * FemMesh.cpp
 *
 *      Date: 2012/05/02
 *      Author: KR
 */

#include "FemMesh.h"

namespace MeshLib {

FemMesh::FemMesh(const std::string &name, const std::vector<Node*> &nodes, const std::vector<Element*> &elements)
	: Mesh(name, nodes, elements)
{
}

FemMesh::FemMesh(const Mesh &mesh)
	: Mesh(mesh.getName(), mesh.getNodes(), mesh.getElements())
{
}

FemMesh::FemMesh(const FemMesh &mesh)
	: Mesh(mesh.getName(), mesh.getNodes(), mesh.getElements())
{
}

FemMesh::FemMesh(const std::string &file_name)
	: Mesh(file_name)
{
}

FemMesh::~FemMesh()
{

}

}

