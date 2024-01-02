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
 * This class represents a 2d quad element in 3d space.
 */
class QuadRule : public FaceRule
{
public:
    /// Constant: The number of base nodes for this element
    static const unsigned n_base_nodes = 4u;

    /// Constant: The geometric type of the element
    static const MeshElemType mesh_elem_type = MeshElemType::QUAD;

    /// Constant: The number of edges
    static const unsigned n_edges = 4;

    /// Constant: The number of neighbors
    static const unsigned n_neighbors = 4;

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

    /// Calculates the area of a quad with straight edges.
    static double computeVolume(Node const* const* element_nodes);
};
}  // namespace MeshLib
