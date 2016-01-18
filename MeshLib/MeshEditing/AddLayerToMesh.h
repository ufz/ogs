/**
 * \file   AddLayerToMesh.h
 * \author Karsten Rink
 * \date   2016-01-18
 * \brief  Definition of AddLayerToMesh class
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ADDLAYERTOMESH_H
#define ADDLAYERTOMESH_H

#include <vector>

namespace MeshLib
{

class Mesh;
class Node;
class Element;

MeshLib::Mesh* addTopLayerToMesh(MeshLib::Mesh const& mesh, double thickness);


} // end namespace MeshLib

#endif //ADDLAYERTOMESH_H
