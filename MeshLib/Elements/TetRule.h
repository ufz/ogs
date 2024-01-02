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

#include "CellRule.h"
#include "Element.h"
#include "MeshLib/MeshEnums.h"

namespace MeshLib
{

class TetRule : public CellRule
{
public:
    /// Constant: The number of base nodes for this element
    static const unsigned n_base_nodes = 4u;

    /// Constant: The geometric type of the element
    static const MeshElemType mesh_elem_type = MeshElemType::TETRAHEDRON;

    /// Constant: The number of faces
    static const unsigned n_faces = 4;

    /// Constant: The number of edges
    static const unsigned n_edges = 6;

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

    /// Calculates the volume of a tetrahedron.
    static double computeVolume(Node const* const* element_nodes);
};
}  // namespace MeshLib
