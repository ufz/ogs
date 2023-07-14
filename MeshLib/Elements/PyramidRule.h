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

#include "CellRule.h"
#include "Element.h"
#include "MeshLib/MeshEnums.h"

namespace MeshLib
{

class PyramidRule : public CellRule
{
public:
    /// Constant: The number of base nodes for this element
    static const unsigned n_base_nodes = 5u;

    /// Constant: The geometric type of the element
    static const MeshElemType mesh_elem_type = MeshElemType::PYRAMID;

    /// Constant: The number of faces
    static const unsigned n_faces = 5;

    /// Constant: The number of edges
    static const unsigned n_edges = 8;

    /// Constant: The number of neighbors
    static const unsigned n_neighbors = 5;

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
    template <std::size_t N>
    static unsigned identifyFace(Node const* const* _nodes,
                                 Node const* nodes[3],
                                 unsigned const face_nodes[5][N])
    {
        for (unsigned i = 0; i < 5; i++)
        {
            unsigned flag(0);
            for (unsigned j = 0; j < 4; j++)
            {
                for (unsigned k = 0; k < 3; k++)
                {
                    if (face_nodes[i][j] != 99 &&
                        _nodes[face_nodes[i][j]] == nodes[k])
                    {
                        flag++;
                    }
                }
            }
            if (flag == 3)
            {
                return i;
            }
        }
        return std::numeric_limits<unsigned>::max();
    }
};
}  // namespace MeshLib
