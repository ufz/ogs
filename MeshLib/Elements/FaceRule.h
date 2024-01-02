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

namespace MeshLib
{

class FaceRule
{
public:
    /// Constant: Dimension of this mesh element
    static const unsigned dimension = 2u;

    /// Returns the face i of the element.
    static const Element* getFace(const Element* e, unsigned i) { return e->getEdge(i); }

    /// Constant: The number of faces
    static const unsigned n_faces = 0;

    /**
     * Checks if the node order of an element is correct by testing surface normals.
     * For 2D elements true is returned if the normal points (roughly) upwards.
     */
    static bool testElementNodeOrder(Element const& e);

    /// \returns the first vector forming the surface' plane
    static Eigen::Vector3d getFirstSurfaceVector(Element const& e);

    /// \returns the second vector forming the surface' plane
    static Eigen::Vector3d getSecondSurfaceVector(Element const& e);

    /// Returns the surface normal of a 2D element.
    static Eigen::Vector3d getSurfaceNormal(Element const& e);

protected:
    /// Returns the ID of an edge given an array of nodes.
    template <typename ElementRule>
    static unsigned identifyFace(Node const* const* element_nodes,
                                 Node const* nodes[ElementRule::dimension])
    {
        for (unsigned i = 0; i < ElementRule::n_edges; i++)
        {
            unsigned flag(0);
            constexpr std::size_t n = sizeof(ElementRule::edge_nodes[0]) /
                                      sizeof(ElementRule::edge_nodes[0][0]);
            for (unsigned j = 0; j < n; j++)
            {
                for (unsigned k = 0; k < ElementRule::dimension; k++)
                {
                    if (element_nodes[ElementRule::edge_nodes[i][j]] ==
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
