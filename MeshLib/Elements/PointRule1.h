/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHLIB_POINTRULE1_H_
#define MESHLIB_POINTRULE1_H_

#include "MeshLib/MeshEnums.h"
#include "Element.h"
#include "VertexRule.h"
#include "EdgeReturn.h"

namespace MeshLib
{

/// A 0d point element.
class PointRule1 : public VertexRule
{
public:
    /// Constant: The number of base nodes for this element
    static const unsigned n_base_nodes = 1u;

    /// Constant: The number of all nodes for this element
    static const unsigned n_all_nodes = 1u;

    /// Constant: The geometric type of the element
    static const MeshElemType mesh_elem_type = MeshElemType::POINT;

    /// Constant: The FEM type of the element
    static const CellType cell_type = CellType::POINT1;

    /// Constant: The number of neighbors
    static const unsigned n_neighbors = 2;

    /// Constant: Local node index table for edge
    static const unsigned edge_nodes[1][1];

    /// Edge rule
    typedef NoEdgeReturn EdgeReturn;

    /// Checks if a point is inside the element.
    ///
    /// Specifically this function tests if the squared euclidean distance
    /// between the points \c _nodes and \c pnt is less then \c eps.
    ///
    /// \param pnt a 3D MathLib::Point3d object
    /// \param eps tolerance for numerical algorithm used for computing the
    /// property
    /// \return true if the point is not outside the element, false otherwise
    static bool isPntInElement(Node const* const* _nodes,
                               MathLib::Point3d const& pnt, double eps);

    /// Tests if the element is geometrically valid.
    static ElementErrorCode validate(const Element* e);

    /// Returns the ID of a face given an array of nodes.
    static unsigned identifyFace(Node const* const*, Node* nodes[1]);

    /// Calculates the length of a line
    static double computeVolume(Node const* const* _nodes);
};
}

#endif  // MESHLIB_POINTRULE1_H_

