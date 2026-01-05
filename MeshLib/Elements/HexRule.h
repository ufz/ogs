// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "CellRule.h"
#include "Element.h"
#include "MeshLib/MeshEnums.h"

namespace MeshLib
{
class HexRule : public CellRule
{
public:
    /// Constant: The number of base nodes for this element
    static const unsigned n_base_nodes = 8u;

    /// Constant: The geometric type of the element
    static const MeshElemType mesh_elem_type = MeshElemType::HEXAHEDRON;

    /// Constant: The number of faces
    static const unsigned n_faces = 6;

    /// Constant: The number of edges
    static const unsigned n_edges = 12;

    /// Constant: The number of neighbors
    static const unsigned n_neighbors = 6;

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

    /// Calculates the volume of a convex hexahedron by partitioning it into six
    /// tetrahedra.
    static double computeVolume(Node const* const* element_nodes);
};
}  // namespace MeshLib
