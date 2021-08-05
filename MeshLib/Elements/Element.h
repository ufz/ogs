/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Definition of the Element class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <limits>
#include <optional>

#include "MathLib/Point3d.h"
#include "MeshLib/Elements/ElementErrorCode.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEnums.h"

namespace MeshLib {

class Node;

/**
 * Virtual base class for mesh elements.
 */
class Element
{
    friend class Mesh;

public:
    /**
     * \brief Tries to add an element e as neighbour to this element.
     * If the elements really are neighbours, the element is added to the
     * neighbour-list and the face id of the neighbour connected to this element
     * is returned. Otherwise the maximum value of the value type is returned.
     */
    std::optional<unsigned> addNeighbor(Element* e);

    /// Returns the length, area or volume of a 1D, 2D or 3D element
    double getContent() const { return _content; }

    /**
     * Get node with local index i where i should be at most the number
     * of nodes of the element
     * @param i local index of node, at most the number of nodes of the
     * element that you can obtain with Element::getNumberOfBaseNodes()
     * @return a pointer to the appropriate (and constant, i.e. not
     * modifiable by the user) instance of class Node or a nullptr
     * @sa Element::getNodeIndex()
     */
    const Node* getNode(unsigned i) const;

    /**
     * (Re)Sets the node of the element.
     * @param idx the index of the pointer to a node within the element
     * @param node a pointer to a node
     */
    void setNode(unsigned idx, Node* node);

    /// Get array of element nodes.
    Node* const* getNodes() const { return _nodes; }

    /// Get dimension of the mesh element.
    virtual unsigned getDimension() const = 0;

    /// Returns the i-th edge of the element.
    virtual const Element* getEdge(unsigned i) const = 0;

    /// Returns the i-th face of the element.
    virtual const Element* getFace(unsigned i) const = 0;

    virtual const Element* getBoundary(unsigned i) const = 0;

    /// Returns the ID of the element.
    virtual std::size_t getID() const final { return _id; }

    virtual unsigned getNumberOfBoundaries() const = 0;

    /// Get the number of edges for this element.
    virtual unsigned getNumberOfEdges() const = 0;

    /// Get the number of faces for this element.
    virtual unsigned getNumberOfFaces() const = 0;

    /// Get the specified neighbor.
    const Element* getNeighbor(unsigned i) const;

    /// Get the number of neighbors for this element.
    virtual unsigned getNumberOfNeighbors() const = 0;

    /**
     * Returns the number of linear nodes.
     */
    virtual unsigned getNumberOfBaseNodes() const = 0;

    /// Returns the number of all nodes including both linear and nonlinear nodes
    virtual unsigned getNumberOfNodes() const = 0;

    /**
     * Get the global index for the Node with local index i.
     * The index i should be at most the number of nodes of the element.
     * @param i local index of Node, at most the number of nodes of the
     * element that you can obtain with Element::getNumberOfBaseNodes()
     * @return the global index or std::numeric_limits<unsigned>::max()
     * @sa Element::getNode()
     */
    std::size_t getNodeIndex(unsigned i) const;

    /**
     * Get the type of the mesh element in geometric context (as a MeshElemType-enum).
     */
    virtual MeshElemType getGeomType() const = 0;

    /**
     * Get the type of the element in context of the finite element method.
     * @return a value of the enum FEMElemType::type
     */
    virtual CellType getCellType() const = 0;

    /// Returns true if the element is located at a boundary (i.e. has at least one face without neighbour)
    virtual bool isBoundaryElement() const;

    /// Returns true if these two indices form an edge and false otherwise
    virtual bool isEdge(unsigned i, unsigned j) const = 0;

    /**
     * Checks if a point is inside the element.
     * @param pnt a 3D MathLib::Point3d object
     * @param eps tolerance for numerical algorithm used or computing the property
     * @return true if the point is not outside the element, false otherwise
     */
    virtual bool isPntInElement(MathLib::Point3d const& pnt, double eps = std::numeric_limits<double>::epsilon()) const = 0;

    /**
     * Tests if the element is geometrically valid.
     */
    virtual ElementErrorCode validate() const = 0;

    /// Destructor
    virtual ~Element();

    /**
     * Method clone is a pure virtual method in the abstract base class Element.
     * It has to be implemented in the derived classes (for instance in class Hex).
     * @return an exact copy of the object
     */
    virtual Element* clone() const = 0;

    /**
     * Constructs a new object polymorphically. This is similar to clone, but
     * accepts new nodes and id.
     * \pre The length of the \c nodes vector is equal to the derived element's
     * total number of nodes.
     */
    virtual Element* clone(Node** nodes, std::size_t id) const = 0;

    /**
     * Computes the length / area / volumen of this element. This is automatically
     * done at initialisation time but can be repeated by calling this function at any time.
     */
    virtual double computeVolume() = 0;

    /// Returns the ID of a face given an array of nodes.
    virtual unsigned identifyFace(Node* nodes[3]) const = 0;

    /**
     * Checks if the node order of an element is correct by testing surface normals.
     */
    virtual bool testElementNodeOrder() const = 0;

    /// Return a specific edge node.
    virtual Node* getEdgeNode(unsigned edge_id, unsigned node_id) const = 0;

#ifndef NDEBUG
    friend std::ostream& operator<<(std::ostream& os, Element const& e);
#endif  // NDEBUG

    /// Dimension of the space, where the element exists.
    unsigned space_dimension_ = 3u;

protected:
    /// Constructor for a generic mesh element without an array of mesh nodes.
    /// @param id     element id
    explicit Element(std::size_t id);

    /// Sets the element ID.
    virtual void setID(std::size_t id) final { _id = id; }

    Node** _nodes;
    std::size_t _id;
    /// Content corresponds to length for 1D, area for 2D, and volume for 3D elements
    double _content;

    Element** _neighbors;
    /// Sets the neighbor over the face with \c face_id to the given \c
    /// neighbor.
    void setNeighbor(Element* neighbor, unsigned const face_id);
}; /* class */

/// Returns true if elem is a neighbour of this element and false otherwise.
bool areNeighbors(Element const* const element, Element const* const other);

/// Returns true if the element has zero length/area/volume.
bool hasZeroVolume(MeshLib::Element const& element);

/// Calculates the center of gravity for the mesh element
MeshLib::Node getCenterOfGravity(MeshLib::Element const& element);

/// Compute the minimum and maximum node distances for this element.
std::pair<double, double> computeSqrNodeDistanceRange(
    MeshLib::Element const& element, bool const check_allnodes = true);

/// Compute the minimum and maximum squared edge length for this element
std::pair<double, double> computeSqrEdgeLengthRange(Element const& element);

/// Let \f$p'\f$ the orthogonal projection to the \f$x\f$-\f$y\f$ plane of the
/// point \c p and \f$e'\f$ the orthogonal projection to the \f$x\f$-\f$y\f$
/// plane of the element \c e.
/// The method checks if \f$p'\f$ is located in \f$e'\f$.  \todo At the moment
/// the test works only for triangle and quad elements.
/// @param p \c MathLib::Point3d is the test point
/// @param e the element that is used for the request
/// @return true if the \f$p' \in e'\f$ and false if \f$p' \notin e'\f$
bool isPointInElementXY(MathLib::Point3d const& p, Element const& e);

/// Returns the position of the given node in the node array of this element.
unsigned getNodeIDinElement(Element const& element, const Node* node);

}  // namespace MeshLib
