/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Definition of the Mesh class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cstdlib>
#include <string>
#include <vector>

#include "BaseLib/Counter.h"
#include "BaseLib/Error.h"

#include "MeshEnums.h"
#include "Properties.h"

namespace ApplicationUtils
{
    class NodeWiseMeshPartitioner;
}

namespace MeshLib
{
    class Node;
    class Element;

/**
 * A basic mesh.
 */
class Mesh : BaseLib::Counter<Mesh>
{
    /* friend functions: */
    friend void removeMeshNodes(MeshLib::Mesh &mesh, const std::vector<std::size_t> &nodes);

    friend class ApplicationUtils::NodeWiseMeshPartitioner;

public:
    /// Constructor using a mesh name and an array of nodes and elements
    /// @param name          Mesh name.
    /// @param nodes         A vector of mesh nodes. In case nonlinear nodes are
    ///                      involved, one should put them after line ones in
    ///                      the vector and set "n_base_nodes" argument.
    /// @param elements      An array of mesh elements.
    /// @param properties    Mesh properties.
    /// @param n_base_nodes  The number of base nodes. This is an optional
    ///                      parameter for nonlinear case.  If the parameter is
    ///                      set to zero, we consider there are no nonlinear
    ///                      nodes.
    Mesh(std::string name,
         std::vector<Node*>
             nodes,
         std::vector<Element*>
             elements,
         Properties const& properties = Properties(),
         const std::size_t n_base_nodes = 0);

    /// Copy constructor
    Mesh(const Mesh &mesh);

    /// Destructor
    virtual ~Mesh();

    /// Add a node to the mesh.
    void addNode(Node* node);

    /// Add an element to the mesh.
    void addElement(Element* elem);

    /// Returns the dimension of the mesh (determined by the maximum dimension over all elements).
    unsigned getDimension() const { return _mesh_dimension; }

    /// Get the node with the given index.
    const Node* getNode(std::size_t idx) const { return _nodes[idx]; }

    /// Get the element with the given index.
    const Element* getElement(std::size_t idx) const { return _elements[idx]; }

    /// Get the minimum edge length over all elements of the mesh.
    double getMinEdgeLength() const { return _edge_length.first; }

    /// Get the maximum edge length over all elements of the mesh.
    double getMaxEdgeLength() const { return _edge_length.second; }

    /// Get the number of elements
    std::size_t getNumberOfElements() const { return _elements.size(); }

    /// Get the number of nodes
    std::size_t getNumberOfNodes() const { return _nodes.size(); }

    /// Get name of the mesh.
    const std::string getName() const { return _name; }

    /// Get the nodes-vector for the mesh.
    std::vector<Node*> const& getNodes() const { return _nodes; }

    /// Get the element-vector for the mesh.
    std::vector<Element*> const& getElements() const { return _elements; }

    /// Resets the IDs of all mesh-elements to their position in the element vector
    void resetElementIDs();

    /// Resets the IDs of all mesh-nodes to their position in the node vector
    void resetNodeIDs();

    /// Changes the name of the mesh.
    void setName(const std::string &name) { this->_name = name; }

    /// Get id of the mesh
    std::size_t getID() const {return _id; }

    /// Get the number of base nodes
    std::size_t getNumberOfBaseNodes() const { return _n_base_nodes; }

    /// Return true if the given node is a basic one (i.e. linear order node)
    bool isBaseNode(std::size_t node_idx) const {return node_idx < _n_base_nodes; }

    /// Return true if the mesh has any nonlinear nodes
    bool isNonlinear() const { return (getNumberOfNodes() != getNumberOfBaseNodes()); }

    MeshLib::Properties & getProperties() { return _properties; }
    MeshLib::Properties const& getProperties() const { return _properties; }

    bool isAxiallySymmetric() const { return _is_axially_symmetric; }
    void setAxiallySymmetric(bool is_axial_symmetric) {
        _is_axially_symmetric = is_axial_symmetric;
    }

protected:
    /// Set the minimum and maximum length over the edges of the mesh.
    void calcEdgeLengthRange();

    /**
     * Resets the connected elements for the node vector, i.e. removes the old information and
     * calls setElementsConnectedToNodes to set the new information.
     * \attention This needs to be called if node neighbourhoods are reset.
     */
    void resetElementsConnectedToNodes();

    /// Sets the dimension of the mesh.
    void setDimension();

    /// Fills in the neighbor-information for nodes (i.e. which element each node belongs to).
    void setElementsConnectedToNodes();

    /// Fills in the neighbor-information for elements.
    /// Note: Using this implementation, an element e can only have neighbors that have the same dimensionality as e.
    void setElementNeighbors();

    void setNodesConnectedByEdges();

    /// Computes the element-connectivity of nodes. Two nodes i and j are
    /// connected if they are shared by an element.
    void setNodesConnectedByElements();

    /// Check if all the nonlinear nodes are stored at the end of the node vector
    void checkNonlinearNodeIDs() const;

    /// Check if the mesh contains any nonlinear element
    bool hasNonlinearElement() const;

    std::size_t const _id;
    unsigned _mesh_dimension;
    /// The minimal and maximal edge length over all elements in the mesh
    std::pair<double, double> _edge_length;
    /// The minimal and maximal distance of nodes within an element over all elements in the mesh
    std::pair<double, double> _node_distance;
    std::string _name;
    std::vector<Node*> _nodes;
    std::vector<Element*> _elements;
    std::size_t _n_base_nodes;
    Properties _properties;

    bool _is_axially_symmetric = false;
}; /* class */



/// Scales the mesh property with name \c property_name by given \c factor.
/// \note The property must be a "double" property.
void scaleMeshPropertyVector(MeshLib::Mesh& mesh,
                             std::string const& property_name,
                             double factor);

/// Creates a new \c PropertyVector in the given mesh and initializes it with
/// the given data. A \c PropertyVector with the same name must not exist.
/// \param mesh A \c Mesh the new \c ProperyVector will be created in.
/// \param name A string that contains the name of the new \c PropertyVector.
/// \param item_type One of the values \c MeshLib::MeshItemType::Cell or \c
/// \c MeshLib::MeshItemType::Node that shows the association of the property
/// values either to \c Element's / cells or \c Node's
/// \param number_of_components the number of components of a property
/// \param values A vector containing the values that are used for
/// initialization.
template <typename T>
void addPropertyToMesh(MeshLib::Mesh& mesh, std::string const& name,
                       MeshLib::MeshItemType item_type,
                       std::size_t number_of_components,
                       std::vector<T> const& values)
{
    if (item_type == MeshLib::MeshItemType::Node)
        if (mesh.getNumberOfNodes() != values.size() / number_of_components)
            OGS_FATAL(
                "Error number of nodes (%u) does not match the number of "
                "tuples (%u).",
                mesh.getNumberOfNodes(), values.size() / number_of_components);
    if (item_type == MeshLib::MeshItemType::Cell)
        if (mesh.getNumberOfElements() != values.size() / number_of_components)
            OGS_FATAL(
                "Error number of elements (%u) does not match the number of "
                "tuples (%u).",
                mesh.getNumberOfElements(),
                values.size() / number_of_components);

    auto* const property = mesh.getProperties().createNewPropertyVector<T>(
        name, item_type, number_of_components);
    if (!property)
    {
        OGS_FATAL("Error while creating PropertyVector \"%s\".", name.c_str());
    }
    property->reserve(values.size());
    std::copy(values.cbegin(), values.cend(), std::back_inserter(*property));
}

/// \returns a PropertyVector of the corresponding type, name on nodes, or
/// cells, or integration points if such exists, or creates a new one.
/// \attention For the integration points the result's size is zero.
/// \see MeshLib::addPropertyToMesh()
template <typename T>
PropertyVector<T>* getOrCreateMeshProperty(Mesh& mesh,
                                           std::string const& property_name,
                                           MeshLib::MeshItemType const item_type,
                                           int const number_of_components)
{
    if (property_name.empty())
        OGS_FATAL(
            "Trying to get or to create a mesh property with empty name.");

    auto numberOfMeshItems = [&mesh, &item_type]() -> std::size_t {
        switch (item_type)
        {
            case MeshLib::MeshItemType::Cell:
                return mesh.getNumberOfElements();
            case MeshLib::MeshItemType::Node:
                return mesh.getNumberOfNodes();
            case MeshLib::MeshItemType::IntegrationPoint:
                return 0;  // For the integration point data the size is
                           // variable
            default:
                OGS_FATAL(
                    "MeshLib::getOrCreateMeshProperty cannot handle other "
                    "types than Node, Cell, or IntegrationPoint.");
        }
        return 0;
    };

    if (mesh.getProperties().existsPropertyVector<T>(property_name))
    {
        auto result =
            mesh.getProperties().template getPropertyVector<T>(property_name);
        assert(result);
        assert(result->size() == numberOfMeshItems() * number_of_components);
        return result;
    }

    auto result = mesh.getProperties().template createNewPropertyVector<T>(
        property_name, item_type, number_of_components);
    assert(result);
    result->resize(numberOfMeshItems() * number_of_components);
    return result;
}

} /* namespace */
