/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

namespace MeshLib
{
class Element;
}  // namespace MeshLib

namespace MeshLib
{
/// Find neighbor elements in given radius of the element.
/// \returns ids of elements for which at least one node lies inside the radius
/// of any node of the given element.
/// \note There are corner-cases not correctly handled by this algorithm:
/// Elements where an edge, but not the edge' nodes, fall inside the search
/// radius are not included.  The algorithm is only considering node to node
/// distances and does not take into account if there is a part of an element
/// (an edge or a face) but not the nodes falls inside the search radius.
///
/// \note For radius 0 only the given element's id is returned.
std::vector<std::size_t> findElementsWithinRadius(Element const& e,
                                                  double const radius);

}  // namespace MeshLib
