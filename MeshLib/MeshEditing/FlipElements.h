/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef FLIPELEMENTS_H_
#define FLIPELEMENTS_H_

#include <vector>
#include "MeshLib/Elements/Element.h"

namespace MeshLib
{
class Mesh;
class Node;

/**
 * Reverses the node order of a 1D / 2D element such that the direction 
 * of a line changes and normals of 2D elements changes its sign.
 */
MeshLib::Element* flipElement(MeshLib::Element const& elem, std::vector<MeshLib::Node*> const& nodes);

/**
 * Reverses the node order of all elements in a 1D / 2D mesh such that 
 * the direction of lines changes and normals of 2D elements changes 
 * their sign.
 */
MeshLib::Mesh* flipMeshElements(MeshLib::Mesh const& mesh);

}

#endif //FLIPELEMENTS_H_
