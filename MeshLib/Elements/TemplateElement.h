/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <array>
#include <limits>

#include "MathLib/Point3d.h"

#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/FaceRule.h"
#include "MeshLib/Elements/CellRule.h"
#include "MeshLib/Elements/ElementErrorCode.h"

namespace MeshLib
{

/**
 * Template for implementing mesh element classes
 *
 * \tparam T_BASE         Base element class, e.g. Face, Cell
 * \tparam ELEMENT_RULE   Geometrical and topological rules of the element
 */
template <class ELEMENT_RULE>
class TemplateElement : public Element
{
public:
    /// Constant: The number of all nodes for this element
    static const unsigned n_all_nodes = ELEMENT_RULE::n_all_nodes;

    /// Constant: The number of base nodes for this element
    static const unsigned n_base_nodes = ELEMENT_RULE::n_base_nodes;

    /// Constant: The dimension of this element
    static const unsigned dimension = ELEMENT_RULE::dimension;

    /**
     * Constructor with an array of mesh nodes.
     *
     * @param nodes  an array of pointers of mesh nodes which form this element
     * @param id     element id
     */
    explicit TemplateElement(
        Node* nodes[n_all_nodes],
        std::size_t id = std::numeric_limits<std::size_t>::max());

    /**
     * Constructor with an array of mesh nodes
     *
     * @param nodes  an array of pointers of mesh nodes which form this element
     * @param id     element id
     */
    explicit TemplateElement(
        std::array<Node*, n_all_nodes> const& nodes,
        std::size_t id = std::numeric_limits<std::size_t>::max());

    /// Copy constructor
    TemplateElement(const TemplateElement &e);

    /// Returns a copy of this object.
    Element* clone() const override { return new TemplateElement(*this); }
    Element* clone(Node** nodes, std::size_t id) const override
    {
        return new TemplateElement(nodes, id);
    }

    /// Get dimension of the mesh element.
    unsigned getDimension() const override { return dimension; }
    /// Returns the edge i of the element.
    const Element* getEdge(unsigned i) const override
    {
        return ELEMENT_RULE::EdgeReturn::getEdge(this, i);
    }

    /// Returns the face i of the element.
    const Element* getFace(unsigned i) const override
    {
        return ELEMENT_RULE::getFace(this, i);
    }

    /// Returns the boundary i of the element.
    const Element* getBoundary(unsigned i) const override
    {
        if constexpr (std::is_convertible<ELEMENT_RULE, FaceRule>::value)
        {
            return ELEMENT_RULE::EdgeReturn::getEdge(this, i);
        }
        if constexpr (std::is_convertible<ELEMENT_RULE, CellRule>::value)
        {
            return ELEMENT_RULE::getFace(this, i);
        }
        OGS_FATAL("TemplateElement::getBoundary for boundary {:d} failed.", i);
    }

    /// Returns the number of boundaries of the element.
    unsigned getNumberOfBoundaries() const override
    {
        if constexpr (std::is_convertible<ELEMENT_RULE, FaceRule>::value)
        {
            return ELEMENT_RULE::n_edges;
        }
        else
        {
            return ELEMENT_RULE::n_faces;
        }
    }

    /// Get the number of edges for this element.
    unsigned getNumberOfEdges() const override { return ELEMENT_RULE::n_edges; }
    /// Get the number of faces for this element.
    unsigned getNumberOfFaces() const override { return ELEMENT_RULE::n_faces; }
    /// Get the number of neighbors for this element.
    unsigned getNumberOfNeighbors() const override
    {
        return ELEMENT_RULE::n_neighbors;
    }

    const Element* getNeighbor(unsigned i) const override
    {
#ifndef NDEBUG
        if (i < ELEMENT_RULE::n_neighbors)
#endif
        {
            return _neighbors[i];
        }
#ifndef NDEBUG
        ERR("Error in MeshLib::TemplateElement::getNeighbor() - Index {:d} "
            "does not exist.",
            i);
        return nullptr;
#endif
    }

    /// Get the number of linear nodes for this element.
    unsigned getNumberOfBaseNodes() const override { return n_base_nodes; }
    /// Get the number of all nodes for this element.
    unsigned getNumberOfNodes() const override { return n_all_nodes; }
    /// Get the type of this element.
    MeshElemType getGeomType() const override
    {
        return ELEMENT_RULE::mesh_elem_type;
    }

    /// Get the FEM type of this element.
    CellType getCellType() const override { return ELEMENT_RULE::cell_type; }
    /// Returns true if these two indices form an edge and false otherwise
    bool isEdge(unsigned idx1, unsigned idx2) const override;

    /**
     * \copydoc MeshLib::Element::isPntInElement()
     *
     * This is actually calling the correct implementation of this function
     * passing the element's nodes.
     */
    bool isPntInElement(
        MathLib::Point3d const& pnt,
        double eps = std::numeric_limits<double>::epsilon()) const override
    {
        return ELEMENT_RULE::isPntInElement(this->_nodes, pnt, eps);
    }

    /**
     * Tests if the element is geometrically valid.
     */
    ElementErrorCode validate() const override
    {
        return ELEMENT_RULE::validate(this);
    }

    /// Returns the ID of a face given an array of nodes.
    unsigned identifyFace(Node* nodes[3]) const override
    {
        return ELEMENT_RULE::identifyFace(this->_nodes, nodes);
    }

    /// Calculates the volume of a convex hexahedron by partitioning it into six tetrahedra.
    double computeVolume() override
    {
        return ELEMENT_RULE::computeVolume(this->_nodes);
    }

    const Node* getNode(unsigned i) const override;
    void setNode(unsigned idx, Node* node) override;
    Node* const* getNodes() const override { return _nodes; }

    /// Return a specific edge node.
    inline Node* getEdgeNode(unsigned edge_id, unsigned node_id) const override
    {
        if (getNumberOfEdges() > 0)
        {
            return const_cast<Node*>(
                this->_nodes[ELEMENT_RULE::edge_nodes[edge_id][node_id]]);
        }

        return nullptr;
    }

    /**
    * Checks if the node order of an element is correct by testing surface normals.
    * For 1D elements this always returns true.
    */
    bool testElementNodeOrder() const override
    {
        return ELEMENT_RULE::testElementNodeOrder(*this);
    }

    double getContent() const override final;
};

}  // namespace MeshLib

#include "TemplateElement-impl.h"
