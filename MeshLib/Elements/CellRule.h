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

#include "MeshLib/Node.h"

namespace MeshLib
{

class Element;

class CellRule
{
public:
    /// Constant: Dimension of this mesh element
    static const unsigned dimension = 3u;

    /**
     * Checks if the node order of an element is correct by testing surface normals.
     * For 3D elements true is returned if the normals of all faces points away from the centre of
     * the element.
     * Note: This method might give wrong results if something else is wrong with the element
     * (non-planar faces, non-convex geometry, possibly zero volume) which causes the calculated
     * center of gravity to lie outside of the actual element
     */
    static bool testElementNodeOrder(Element const& e);

protected:
    /// Returns the ID of a face given an array of nodes.
    template <typename ElementRule>
    static unsigned identifyFace(Node const* const* element_nodes,
                                 Node const* nodes[ElementRule::dimension])
    {
        for (unsigned i = 0; i < ElementRule::n_faces; i++)
        {
            unsigned flag(0);
            constexpr std::size_t n = sizeof(ElementRule::face_nodes[0]) /
                                      sizeof(ElementRule::face_nodes[0][0]);
            for (unsigned j = 0; j < n; j++)
            {
                for (unsigned k = 0; k < ElementRule::dimension; k++)
                {
                    if (ElementRule::face_nodes[i][j] != 99 &&
                        element_nodes[ElementRule::face_nodes[i][j]] ==
                            nodes[k])
                    {
                        flag++;
                    }
                }
            }
            if (flag == ElementRule::dimension)
            {
                return i;
            }
        }
        return std::numeric_limits<unsigned>::max();
    }
}; /* class */

}  // namespace MeshLib
