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

#include "EdgeRule.h"
#include "Element.h"
#include "MeshLib/MeshEnums.h"

namespace MeshLib
{
/**
 * A 1d Edge or Line Element.
 */
class LineRule : public EdgeRule
{
public:
    /// Constant: The number of base nodes for this element
    static const unsigned n_base_nodes = 2u;

    /// Constant: The geometric type of the element
    static const MeshElemType mesh_elem_type = MeshElemType::LINE;

    /// Constant: The number of neighbors
    static const unsigned n_neighbors = 2;

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

    /// Returns the ID of a face given an array of nodes.
    static unsigned identifyFace(Node const* const* /*element_nodes*/,
                                 Node const* nodes[1]);

    /// Calculates the length of a line.
    static double computeVolume(Node const* const* element_nodes);
};
}  // namespace MeshLib
