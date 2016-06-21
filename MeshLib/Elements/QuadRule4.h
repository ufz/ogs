/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef QUADRULE4_H_
#define QUADRULE4_H_

#include "MeshLib/MeshEnums.h"
#include "Element.h"
#include "FaceRule.h"
#include "EdgeReturn.h"
#include "Line.h"

#ifdef _MSC_VER
#include "meshlib_export.h"
#else
#define MESHLIB_EXPORT
#endif

namespace MeshLib
{

/**
 * This class represents a 2d quadrilateral element. The following sketch shows the node and edge numbering.
 * @anchor QuadNodeAndEdgeNumbering
 * @code
 *              2
 *        3-----------2
 *        |           |
 *        |           |
 *       3|           |1
 *        |           |
 *        |           |
 *        0-----------1
 *              0
 * @endcode
 */
class QuadRule4 : public FaceRule
{
public:
    /// Constant: The number of base nodes for this element
    static const unsigned n_base_nodes = 4u;

    /// Constant: The number of all nodes for this element
    static const unsigned n_all_nodes = 4u;

    /// Constant: The geometric type of the element
    static const MeshElemType mesh_elem_type = MeshElemType::QUAD;

    /// Constant: The FEM type of the element
    static const CellType cell_type = CellType::QUAD4;

    /// Constant: The number of edges
    static const unsigned n_edges = 4;

    /// Constant: The number of neighbors
    static const unsigned n_neighbors = 4;

    /// Constant: Local node index table for edge
    static MESHLIB_EXPORT const unsigned edge_nodes[4][2];

    /// Returns the i-th edge of the element.
    typedef LinearEdgeReturn EdgeReturn;

    /**
     * Checks if a point is inside the element.
     * @param pnt a 3D MathLib::Point3d object
     * @param eps tolerance for numerical algorithm used or computing the property
     * @return true if the point is not outside the element, false otherwise
     */
    static bool isPntInElement(Node const* const* _nodes, MathLib::Point3d const& pnt, double eps);

    /**
     * Tests if the element is geometrically valid.
     */
    static ElementErrorCode validate(const Element* e);

    /// Returns the ID of a face given an array of nodes.
    static unsigned identifyFace(Node const* const*, Node* nodes[3]);

    /// Calculates the volume of a convex hexahedron by partitioning it into six tetrahedra.
    static double computeVolume(Node const* const* _nodes);

}; /* class */

} /* namespace */

#endif /* QUADRULE4_H_ */

