/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "Element.h"
#include "FaceRule.h"
#include "MeshLib/MeshEnums.h"

namespace MeshLib
{
/**
 * This class represents a 2d triangle element in 3d space.
 */
class TriRule : public FaceRule
{
public:
    /// Constant: The number of base nodes for this element
    static const unsigned n_base_nodes = 3u;

    /// Constant: The geometric type of the element
    static const MeshElemType mesh_elem_type = MeshElemType::TRIANGLE;

    /// Constant: The number of edges
    static const unsigned n_edges = 3;

    /// Constant: The number of neighbors
    static const unsigned n_neighbors = 3;

    /**
     * \copydoc MeshLib::Element::isPntInElement()
     * \param nodes the nodes of the element.
     */
    static bool isPntInElement(Node const* const* nodes,
                               MathLib::Point3d const& pnt, double eps);

    /**
     * Tests if the element is geometrically valid.
     */
    static ElementErrorCode validate(const Element* e);

    /// Calculates the area of the triangle with straight edges.
    static double computeVolume(Node const* const* element_nodes);
};
}  // namespace MeshLib
