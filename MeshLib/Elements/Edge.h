/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Definition of the Edge class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef EDGE_H_
#define EDGE_H_

#include <limits>

#include "Element.h"
#include "Node.h"

#include "MathTools.h"


namespace MeshLib
{

/**
 * Virtual base class for 1d mesh elements.
 */
class Edge : public Element
{
public:
    /// Get the length of this 1d element.
    double getLength() const { return _length; };

    /// Returns the length, area or volume of a 1D, 2D or 3D element
    double getContent() const { return _length; };

    /// Get dimension of the mesh element.
    unsigned getDimension() const { return 1; };

    /// Returns the edge i of the element.
    const Element* getEdge(unsigned /*i*/) const { return nullptr; };

    /// Returns the face i of the element.
    const Element* getFace(unsigned /*i*/) const { return nullptr; };

    /// 1D elements have no edges
    unsigned getNEdges() const { return 0; };

    /// Get the number of nodes for face i.
    unsigned getNFaceNodes(unsigned /*i*/) const { return 0; };

    /// Get the number of faces for this element.
    unsigned getNFaces() const { return 0; };

    /// Get the number of neighbors for this element.
    unsigned getNNeighbors() const { return 0; };

    /// Returns true if these two indices form an edge and false otherwise
    bool isEdge(unsigned idx1, unsigned idx2) const
    {
        if (0==idx1 && 1==idx2) return true;
        if (1==idx1 && 0==idx2) return true;
        return false;
    }

    /// Destructor
    virtual ~Edge() {};

    /**
     * This method is pure virtual and is inherited from class @sa Element.
     * It has to be implemented in the derived classes of class Face!
     * @return a copy of the object
     */
    virtual Element* clone() const = 0;

protected:
    /// 1D elements have no edges.
    Node* getEdgeNode(unsigned /*edge_id*/, unsigned /*node_id*/) const { return nullptr; };

    /// Returns the ID of a face given an array of nodes (but is not applicable for edges!).
    unsigned identifyFace(Node* [3]/*nodes[3]*/) const { return std::numeric_limits<unsigned>::max(); };

    /// Constructor for a generic mesh element without an array of mesh nodes.
    Edge(unsigned value = 0);

    double _length;

}; /* class */

} /* namespace */

#endif /* EDGE_H_ */

