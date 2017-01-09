/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>
#include "MeshLib/Elements/Element.h"

namespace MeshLib
{
class Mesh;
class Node;

/**
 * Creates a copy of an 1d / 2d element where the node order is reversed,
 * such that the direction of a line changes and normals of 2D elements
 * changes its sign.
 * @param elem original element
 * @param nodes node vector used for the copy of the element
 * @return a flipped copy of the element
 */
std::unique_ptr<MeshLib::Element> createFlippedElement(MeshLib::Element const& elem, std::vector<MeshLib::Node*> const& nodes);

/**
 * Creates a copy of a 1d / 2d mesh where the node order of all elements
 * is reversed such that the direction of lines changes and normals of
 * 2D elements changes their sign.
 * @param mesh input mesh
 * @return a flipped copy of the input mesh
 */
std::unique_ptr<MeshLib::Mesh> createFlippedMesh(MeshLib::Mesh const& mesh);

}
