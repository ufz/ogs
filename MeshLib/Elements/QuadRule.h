/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
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

    /// Calculates the volume of a convex hexahedron by partitioning it into six
    /// tetrahedra.
    static double computeVolume(Node const* const* _nodes);

protected:
    /// Returns the ID of a face given an array of nodes.
    template <std::size_t N>
    static unsigned identifyFace(Node const* const* _nodes,
                                 Node const* nodes[3],
                                 const unsigned edge_nodes[4][N])
    {
        for (unsigned i = 0; i < 4; i++)
        {
            unsigned flag(0);
            for (unsigned j = 0; j < 2; j++)
            {
                for (unsigned k = 0; k < 2; k++)
                {
                    if (_nodes[edge_nodes[i][j]] == nodes[k])
                    {
                        flag++;
                    }
                }
            }
            if (flag == 2)
            {
                return i;
            }
        }
        return std::numeric_limits<unsigned>::max();
    }
};
}  // namespace MeshLib
