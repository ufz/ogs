/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PRISMRULE6_H_
#define PRISMRULE6_H_

#include "MeshLib/MeshEnums.h"
#include "Element.h"
#include "EdgeReturn.h"
#include "CellRule.h"

namespace MeshLib
{

/**
 * This class represents a 3d prism element. The following sketch shows the node and edge numbering.
 * @anchor Prism6NodeAndEdgeNumbering
 * @code
 *            5
 *           / \
 *          / : \
 *        8/  :  \7
 *        /   :5  \
 *       /    :  6 \
 *      3-----------4
 *      |     :     |
 *      |     2     |
 *      |    . .    |
 *     3|   .   .   |4
 *      | 2.     .1 |
 *      | .       . |
 *      |.         .|
 *      0-----------1
 *            0
 *
 * @endcode
 */
class PrismRule6 : public CellRule
{
public:
    /// Constant: The number of base nodes for this element
    static const unsigned n_base_nodes = 6u;

    /// Constant: The number of all nodes for this element
    static const unsigned n_all_nodes = 6u;

    /// Constant: The geometric type of the element
    static const MeshElemType mesh_elem_type = MeshElemType::PRISM;

    /// Constant: The FEM type of the element
    static const CellType cell_type = CellType::PRISM6;

    /// Constant: The number of faces
    static const unsigned n_faces = 5;

    /// Constant: The number of edges
    static const unsigned n_edges = 9;

    /// Constant: The number of neighbors
    static const unsigned n_neighbors = 5;

    /// Constant: Local node index table for faces
    static const unsigned face_nodes[5][4];

    /// Constant: Local node index table for edge
    static const unsigned edge_nodes[9][2];

    /// Constant: Table for the number of nodes for each face
    static const unsigned n_face_nodes[5];

    /// Returns the i-th edge of the element.
    typedef LinearEdgeReturn EdgeReturn;

    /// Returns the i-th face of the element.
    static const Element* getFace(const Element* e, unsigned i);

    /**
     * \copydoc MeshLib::Element::isPntInElement()
     * @param nodes the nodes of the element.
     */
    static bool isPntInElement(Node const* const* nodes,
                               MathLib::Point3d const& pnt, double eps);

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

#endif /* PRISMRULE6_H_ */

