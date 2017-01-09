/**
 * \file   MeshRevision.h
 * \author Karsten Rink
 * \date   2014-02-14
 * \brief  Definition of the MeshRevision class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHREVISION_H_
#define MESHREVISION_H_

#include <array>
#include <limits>
#include <string>
#include <vector>

// forward declaration
namespace MeshLib {
    class Mesh;
    class Node;
    class Element;
}

namespace MeshLib {

/**
 * Collapses nodes with a distance smaller min_distance and
 * reduces elements accordingly.
 */
class MeshRevision {

public:
    /**
     * Constructor
     * @param mesh The mesh which is being revised. Note that node IDs in mesh are changed during computation but are resetted after the algorithms implemented here are finished
     */
    MeshRevision(MeshLib::Mesh &mesh);

    virtual ~MeshRevision() {}

    /**
     * Collapsed all nodes with distance < eps but ignores elements
     * (i.e. elements with collapsed nodes may result)
     * This is implicitely called when calling simplifyMesh(), so it does not need to be
     * called separately when using simplifyMesh().
     */
    MeshLib::Mesh* collapseNodes(const std::string &new_mesh_name, double eps);

    /// Returns the number of potentially collapsable nodes
    unsigned getNumberOfCollapsableNodes(double eps = std::numeric_limits<double>::epsilon()) const;

    /// Designates nodes to be collapsed by setting their ID to the index of the node they will get merged with.
    std::vector<std::size_t> collapseNodeIndices(double eps) const;

    /**
     * Create a new mesh where all nodes with a distance < eps from each other
     * are collapsed.
     * Elements are adjusted accordingly and elements with nonplanar faces are
     * subdivided into geometrically correct elements.
     * @param new_mesh_name New name.
     * @param eps           Minimum distance for nodes not to be collapsed
     * @param min_elem_dim  Minimum dimension of elements to be inserted into
     *                      new mesh (i.e. min_elem_dim=3 will prevent the new
     *                      mesh to contain 2D elements)
     */
    MeshLib::Mesh* simplifyMesh(const std::string& new_mesh_name, double eps,
                                unsigned min_elem_dim = 1);

    /**
     * Create a new mesh where all elements with nonplanar faces are subdivided into simpler
     * element types. This method does not collapse or remove any nodes.
     */
    MeshLib::Mesh* subdivideMesh(const std::string &new_mesh_name) const;

private:
    /// Constructs a new node vector for the resulting mesh by removing all nodes whose ID indicates they need to be merged/removed.
    std::vector<MeshLib::Node*> constructNewNodesArray(
        const std::vector<std::size_t> &id_map) const;

    /// Calculates the number of unique nodes in an element (i.e. uncollapsed nodes)
    unsigned getNumberOfUniqueNodes(MeshLib::Element const*const element) const;

    /// Resets the node IDs of the source mesh (needs to be called after everything is done).
    void resetNodeIDs();

    /// Subdivides an element if it has a face that is not coplanar
    /// @param element the element that will be subdivided
    /// @param nodes vector containing the nodes the elements originated by the
    /// subdivision are based on
    /// @param elements vector of MeshLib::Elements; the elements originated by
    /// the subdivision will be inserted into elements
    /// @return the number of elements originated by the subdivision
    std::size_t subdivideElement(MeshLib::Element const*const element,
        std::vector<MeshLib::Node*> const& nodes,
        std::vector<MeshLib::Element*> & elements) const;

    // Revises an element by removing collapsed nodes, using the nodes vector from the result mesh.
    std::size_t reduceElement(MeshLib::Element const*const element,
        unsigned n_unique_nodes,
        const std::vector<MeshLib::Node*> &nodes,
        std::vector<MeshLib::Element*> &elements,
        unsigned min_elem_dim) const;

    /// Cleans up all nodes and elements if something went wrong
    void cleanUp(std::vector<MeshLib::Node*> &nodes,
        std::vector<MeshLib::Element*> &new_elements) const;

    /// Subdivides a nonplanar quad into two triangles
    unsigned subdivideQuad(MeshLib::Element const*const quad,
        const std::vector<MeshLib::Node*> &nodes,
        std::vector<MeshLib::Element*> &new_elements) const;

    /// Subdivides a Hex with nonplanar faces into tets
    unsigned subdivideHex(MeshLib::Element const*const hex,
        const std::vector<MeshLib::Node*> &nodes,
        std::vector<MeshLib::Element*> &new_elements) const;

    /// Subdivides a pyramid with a nonplanar base into two tets
    unsigned subdividePyramid(MeshLib::Element const*const pyramid,
        const std::vector<MeshLib::Node*> &nodes,
        std::vector<MeshLib::Element*> &new_elements) const;

    /// Subdivides a prism with nonplanar quad faces into two tets
    unsigned subdividePrism(MeshLib::Element const*const prism,
        const std::vector<MeshLib::Node*> &nodes,
        std::vector<MeshLib::Element*> &new_elements) const;

    /// Creates a line element from the first two unique nodes found in the element (element *should* have exactly two unique nodes!)
    MeshLib::Element* constructLine(MeshLib::Element const*const element,
        const std::vector<MeshLib::Node*> &nodes) const;
    /// Creates a triangle element from the first three unique nodes found in the element (element *should* have exactly three unique nodes!)
    MeshLib::Element* constructTri(MeshLib::Element const*const element,
        const std::vector<MeshLib::Node*> &nodes) const;
    /// Creates a quad or a tet, depending if the four nodes being coplanar or not (element *should* have exactly four unique nodes!)
    MeshLib::Element* constructFourNodeElement(
        MeshLib::Element const*const element,
        const std::vector<MeshLib::Node*> &nodes,
        unsigned min_elem_dim = 1) const;

    /**
     * Reduces a hexahedron element by removing collapsed nodes and constructing one or more new elements from the remaining nodes.
     * @return The number of newly created elements
     */
    unsigned reduceHex(MeshLib::Element const*const hex,
        unsigned n_unique_nodes,
        const std::vector<MeshLib::Node*> &nodes,
        std::vector<MeshLib::Element*> &new_elements,
        unsigned min_elem_dim) const;
    /// Reduces a pyramid element by removing collapsed nodes and constructing a new elements from the remaining nodes.
    void reducePyramid(MeshLib::Element const*const pyramid,
        unsigned n_unique_nodes,
        const std::vector<MeshLib::Node*> &nodes,
        std::vector<MeshLib::Element*> &new_elements,
        unsigned min_elem_dim) const;
    /**
     * Reduces a prism element by removing collapsed nodes and constructing one or two new elements from the remaining nodes.
     * @return The number of newly created elements
     */
    unsigned reducePrism(MeshLib::Element const*const prism,
        unsigned n_unique_nodes,
        std::vector<MeshLib::Node*> const& nodes,
        std::vector<MeshLib::Element*> & new_elements,
        unsigned min_elem_dim) const;

    // In an element with 5 unique nodes, return the node that will be the top of the resulting pyramid
    unsigned findPyramidTopNode(MeshLib::Element const& element,
        std::array<std::size_t,4> const& base_node_ids) const;

    /// Lookup-table for returning the diametral node id of the given node id in a Hex
    unsigned lutHexDiametralNode(unsigned id) const;

    /// Lookup-table for returning four nodes connected to the two nodes (id1, id2) forming an edge in a Hex
    const std::array<unsigned,4> lutHexCuttingQuadNodes(unsigned id1, unsigned id2) const;

    /// When a hex is subdivided into two prisms, this returns the nodes of the hex edge that will serve as the back of one of the prisms.
    const std::pair<unsigned, unsigned> lutHexBackNodes(
        unsigned i, unsigned j, unsigned k, unsigned l) const;

    /// Lookup-table for returning the third node of bottom or top triangle given the other two
    unsigned lutPrismThirdNode(unsigned id1, unsigned id2) const;

    /// The original mesh used for constructing the class
    Mesh& _mesh;

    static const std::array<unsigned,8> _hex_diametral_nodes;
};

}

#endif /* MESHREVISION_H_ */
