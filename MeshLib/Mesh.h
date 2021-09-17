/**
 * \file
 * \author Karsten Rink
 * \date   2012-05-02
 * \brief  Definition of the Mesh class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cstdlib>
#include <memory>
#include <string>
#include <vector>

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
class Mesh
{
    /* friend functions: */
    friend void removeMeshNodes(Mesh& mesh,
                                const std::vector<std::size_t>& nodes);

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

    /// Finds the maximum id among all of the base nodes.
    void recalculateMaxBaseNodeId();

    /// Changes the name of the mesh.
    void setName(const std::string &name) { this->_name = name; }

    /// Get id of the mesh
    std::size_t getID() const {return _id; }

    /// Get the number of base nodes
    std::size_t getNumberOfBaseNodes() const { return _n_base_nodes; }

    /// Check if the mesh contains any nonlinear element.
    bool hasNonlinearElement() const;

    Properties& getProperties() { return _properties; }
    Properties const& getProperties() const { return _properties; }

    bool isAxiallySymmetric() const { return _is_axially_symmetric; }
    void setAxiallySymmetric(bool is_axial_symmetric) {
        _is_axially_symmetric = is_axial_symmetric;
    }

protected:
    /// Set the minimum and maximum length over the edges of the mesh.
    void calcEdgeLengthRange();

    /// Sets the dimension of the mesh.
    void setDimension();

    /// Fills in the neighbor-information for nodes (i.e. which element each node belongs to).
    void setElementsConnectedToNodes();

    /// Fills in the neighbor-information for elements.
    /// Note: Using this implementation, an element e can only have neighbors that have the same dimensionality as e.
    void setElementNeighbors();

    /// Check if all the nonlinear nodes are stored at the end of the node vector
    void checkNonlinearNodeIDs() const;

    std::size_t const _id;
    unsigned _mesh_dimension;
    /// The minimal and maximal edge length over all elements in the mesh
    std::pair<double, double> _edge_length = {
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()};
    /// The minimal and maximal distance of nodes within an element over all elements in the mesh
    std::pair<double, double> _node_distance;
    std::string _name;
    std::vector<Node*> _nodes;
    std::vector<Element*> _elements;
    std::size_t _n_base_nodes;
    Properties _properties;

    bool _is_axially_symmetric = false;
}; /* class */

/// Computes the element-connectivity of nodes. Two nodes i and j are
/// connected if they are shared by an element.
std::vector<std::vector<Node*>> calculateNodesConnectedByElements(
    Mesh const& mesh);

/// Meshes are equal if their id's are equal.
inline bool operator==(Mesh const& a, Mesh const& b)
{
    return a.getID() == b.getID();
}

inline bool operator!=(Mesh const& a, Mesh const& b)
{
    return !(a == b);
}

/// Scales the mesh property with name \c property_name by given \c factor.
/// \note The property must be a "double" property.
void scaleMeshPropertyVector(Mesh& mesh,
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
void addPropertyToMesh(Mesh& mesh, std::string const& name,
                       MeshItemType item_type, std::size_t number_of_components,
                       std::vector<T> const& values)
{
    if (item_type == MeshItemType::Node)
    {
        if (mesh.getNumberOfNodes() != values.size() / number_of_components)
        {
            OGS_FATAL(
                "Error number of nodes ({:d}) does not match the number of "
                "tuples ({:d}).",
                mesh.getNumberOfNodes(), values.size() / number_of_components);
        }
    }
    if (item_type == MeshItemType::Cell)
    {
        if (mesh.getNumberOfElements() != values.size() / number_of_components)
        {
            OGS_FATAL(
                "Error number of elements ({:d}) does not match the number of "
                "tuples ({:d}).",
                mesh.getNumberOfElements(),
                values.size() / number_of_components);
        }
    }

    auto* const property = mesh.getProperties().createNewPropertyVector<T>(
        name, item_type, number_of_components);
    if (!property)
    {
        OGS_FATAL("Error while creating PropertyVector '{:s}'.", name);
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
                                           MeshItemType const item_type,
                                           int const number_of_components)
{
    if (property_name.empty())
    {
        OGS_FATAL(
            "Trying to get or to create a mesh property with empty name.");
    }

    auto numberOfMeshItems = [&mesh, &item_type]() -> std::size_t {
        switch (item_type)
        {
            case MeshItemType::Cell:
                return mesh.getNumberOfElements();
            case MeshItemType::Node:
                return mesh.getNumberOfNodes();
            case MeshItemType::IntegrationPoint:
                return 0;  // For the integration point data the size is
                           // variable
            default:
                OGS_FATAL(
                    "getOrCreateMeshProperty cannot handle other "
                    "types than Node, Cell, or IntegrationPoint.");
        }
        return 0;
    };

    if (mesh.getProperties().existsPropertyVector<T>(property_name))
    {
        auto result =
            mesh.getProperties().template getPropertyVector<T>(property_name);
        assert(result);
        if (item_type != MeshItemType::IntegrationPoint)
        {
            // Test the size if number of mesh items is known, which is not the
            // case for the integration point data.
            assert(result->size() ==
                   numberOfMeshItems() * number_of_components);
        }
        return result;
    }

    auto result = mesh.getProperties().template createNewPropertyVector<T>(
        property_name, item_type, number_of_components);
    assert(result);
    result->resize(numberOfMeshItems() * number_of_components);
    return result;
}

/// Returns the material ids property vector defined on the mesh.
///
/// The material ids are always an \c int property named "MaterialIDs".
/// If the property does not exists (or is of different type), a nullptr is
/// returned.
PropertyVector<int> const* materialIDs(Mesh const& mesh);

/// Creates a new mesh from a vector of elements.
///
/// \note The elements are owned by the returned mesh object as well as the
/// nodes and will be destructed together with the mesh.
std::unique_ptr<Mesh> createMeshFromElementSelection(
    std::string mesh_name, std::vector<Element*> const& elements);

}  // namespace MeshLib
