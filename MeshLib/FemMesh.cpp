/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file FemMesh.cpp
 *
 * Created on 2012-05-02 by Karsten Rink
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

FemMesh::~FemMesh()
{

}

}

