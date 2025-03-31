/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 26, 2025, 4:28 PM
 */

#include "MergeMeshToBulkMesh.h"

#include <Eigen/Dense>
#include <algorithm>
#include <boost/functional/hash.hpp>
#include <boost/range/combine.hpp>
#include <cmath>
#include <memory>
#include <numeric>
#include <range/v3/view/filter.hpp>
#include <range/v3/view/zip.hpp>
#include <ranges>
#include <tuple>
#include <unordered_map>
#include <vector>

#include "BaseLib/Logging.h"
#include "GeoLib/AABB.h"
#include "GeoLib/OctTree.h"
#include "MeshLib/Elements/Elements.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Properties.h"
#include "MeshLib/Utils/getOrCreateMeshProperty.h"

namespace MeshToolsLib
{
template <typename T>
bool createMergedPropertyVector(
    MeshLib::Mesh& merged_mesh,
    std::unordered_map<std::string, double>& initial_value_dict,
    MeshLib::PropertyVector<T> const* const pv_bulk_mesh,
    MeshLib::Properties const& properties_bulk_mesh)
{
    if (pv_bulk_mesh == nullptr)
    {
        return false;
    }

    if (pv_bulk_mesh->getPropertyName() == "vtkGhostType")
    {
        // Do nothing
        return true;
    }

    auto const item_type = pv_bulk_mesh->getMeshItemType();

    auto const pv_name = pv_bulk_mesh->getPropertyName();

    auto const pv_num_components = pv_bulk_mesh->getNumberOfGlobalComponents();

    if (pv_name == "OGS_VERSION" || pv_name == "IntegrationPointMetaData")
    {
        auto new_pv = MeshLib::getOrCreateMeshProperty<T>(
            merged_mesh, pv_name, item_type, pv_num_components);
        new_pv->resize(pv_bulk_mesh->size());

        std::copy(pv_bulk_mesh->begin(), pv_bulk_mesh->end(), new_pv->begin());
        return true;
    }

    if (item_type == MeshLib::MeshItemType::Node)
    {
        auto new_pv = MeshLib::getOrCreateMeshProperty<T>(
            merged_mesh, pv_name, item_type, pv_num_components);
        new_pv->resize(merged_mesh.getNumberOfNodes() * pv_num_components);

        std::copy(pv_bulk_mesh->begin(), pv_bulk_mesh->end(), new_pv->begin());
        if (pv_num_components > 1)
        {
            if (pv_name.find("sigma") != std::string::npos)
            {
                std::vector<double> sigma0(pv_num_components, 0.0);
                sigma0[0] = initial_value_dict["sxx"];
                sigma0[1] = initial_value_dict["syy"];
                sigma0[2] = initial_value_dict["szz"];

                // Number of the vector elements of the mesh to be merged.
                auto const num_elements =
                    (new_pv->size() - pv_bulk_mesh->size()) / pv_num_components;
                for (std::size_t i = 0; i < num_elements; ++i)
                {
                    new_pv->insert(new_pv->begin() + pv_bulk_mesh->size(),
                                   sigma0.begin(), sigma0.end());
                }
                return true;
            }
            std::fill(new_pv->begin() + pv_bulk_mesh->size(), new_pv->end(), 0);
            return true;
        }

        // Map possible pv_name values to their corresponding dictionary
        // keys
        const std::unordered_map<std::string, std::string> pv_to_dict_key = {
            {"pressure", "p"},
            {"p", "p"},
            {"gas_pressure", "pg"},
            {"pg", "pg"},
            {"capillary_pressure", "pc"},
            {"pc", "pc"},
            {"temperature", "T"},
            {"T", "T"}};

        double value = 0.0;
        auto it = pv_to_dict_key.find(pv_name);
        if (it != pv_to_dict_key.end() &&
            initial_value_dict.contains(it->second))
        {
            value = initial_value_dict[it->second];
        }
        std::fill(new_pv->begin() + pv_bulk_mesh->size(), new_pv->end(), value);
        return true;
    }

    if (item_type == MeshLib::MeshItemType::Cell)
    {
        auto new_pv = MeshLib::getOrCreateMeshProperty<T>(
            merged_mesh, pv_name, item_type, pv_num_components);
        new_pv->resize(merged_mesh.getNumberOfElements() * pv_num_components);

        std::copy(pv_bulk_mesh->begin(), pv_bulk_mesh->end(), new_pv->begin());

        if (pv_name == "MaterialIDs")
        {
            std::fill(new_pv->begin() + pv_bulk_mesh->size(), new_pv->end(),
                      initial_value_dict["mat_id"]);
            return true;
        }
        std::fill(new_pv->begin() + pv_bulk_mesh->size(), new_pv->end(), 0.0);
        return true;
    }

    return true;
}

std::vector<MeshLib::Node*> findNodesInBoundedDomain(
    std::vector<MeshLib::Node*> const& nodes, GeoLib::AABB const& aabb)
{
    return nodes |
           ranges::views::filter(
               [&](MeshLib::Node* n)
               { return aabb.containsPoint(n->asEigenVector3d(), 1e-16); }) |
           ranges::to<std::vector<MeshLib::Node*>>();
}

struct NodesPartitionResult
{
    std::vector<MeshLib::Node*> paired_nodes;
    std::optional<std::vector<MeshLib::Node*>> non_paired_nodes;
};

// Assumes at least one element in \c node_vector can be paired with an element
// in \c tool_node_vector based on identical coordinates.
// Partitions \c node_vector into two vectors: one containing nodes that can be
// paired with elements in \c tool_node_vector, and another may containing the
// remaining nodes or empty depending on the argument \c
// return_non_paired_nodes.
NodesPartitionResult partitionNodesByCoordinateMatch(
    std::vector<MeshLib::Node*> const& node_vector,
    std::vector<MeshLib::Node*> const& tool_node_vector,
    GeoLib::AABB const& aabb,
    bool const return_non_paired_nodes = true)
{
    auto oct_tree = std::unique_ptr<GeoLib::OctTree<MeshLib::Node, 16>>(
        GeoLib::OctTree<MeshLib::Node, 16>::createOctTree(
            aabb.getMinPoint(), aabb.getMaxPoint(), 1e-16));

    // Push all tool nodes into oct_tree:
    for (auto const node : tool_node_vector)
    {
        MeshLib::Node* node_ptr = nullptr;
        oct_tree->addPoint(node, node_ptr);
    }

    // Find the paired nodes in the node_vector
    std::vector<MeshLib::Node*> paired_nodes;
    std::vector<MeshLib::Node*> other_nodes;
    for (auto node : node_vector)
    {
        MeshLib::Node* node_ptr = nullptr;
        if (oct_tree->addPoint(node, node_ptr))
        {
            if (return_non_paired_nodes)
            {
                other_nodes.push_back(node);
            }
            continue;
        }
        paired_nodes.push_back(node);
    }

    if (return_non_paired_nodes)
    {
        return {paired_nodes, other_nodes};
    }

    return {paired_nodes, std::nullopt};
}
// Custom hash for Eigen::Vector3d
struct Vector3dHash
{
    size_t operator()(const Eigen::Vector3d& v) const
    {
        return boost::hash_range(v.data(), v.data() + 3);
    }
};

struct VectorEqual
{
    bool operator()(const Eigen::Vector3d& a, const Eigen::Vector3d& b) const
    {
        return a.isApprox(b, 1e-16);
    }
};

std::vector<std::size_t> findNodeIdMapping(
    const std::vector<MeshLib::Node*>& nodes_a,
    const std::vector<MeshLib::Node*>& nodes_b)
{
    // Preallocate hash map for performance (avoids rehashing)
    std::unordered_map<Eigen::Vector3d, std::size_t, Vector3dHash, VectorEqual>
        coord_to_map(nodes_b.size());

    for (MeshLib::Node* const node : nodes_b)
    {
        coord_to_map[node->asEigenVector3d()] = node->getID();
    }

    std::vector<std::size_t> id_mapping;
    id_mapping.reserve(nodes_a.size());

    for (MeshLib::Node* const node : nodes_a)
    {
        auto it = coord_to_map.find(node->asEigenVector3d());
        if (it != coord_to_map.end())
        {
            id_mapping.push_back(it->second);  // Found match
        }
        else
        {
            OGS_FATAL("No matching node found for ID {:d} ", node->getID());
        }
    }

    return id_mapping;
}

std::unique_ptr<MeshLib::Mesh> mergeMeshToBulkMesh(
    MeshLib::Mesh const& bulk_mesh, MeshLib::Mesh const& other_mesh,
    std::unordered_map<std::string, double>& initial_value_dict)
{
    auto const& other_mesh_nodes = other_mesh.getNodes();
    GeoLib::AABB aabb(other_mesh_nodes.begin(), other_mesh_nodes.end());

    auto const& bulk_mesh_nodes = bulk_mesh.getNodes();
    auto const bulk_nodes_in_aabb =
        findNodesInBoundedDomain(bulk_mesh_nodes, aabb);

    // Find the interface nodes in the bulk nodes
    auto const pn_bulk_mesh = partitionNodesByCoordinateMatch(
        bulk_nodes_in_aabb, other_mesh_nodes, aabb,
        false /*return_non_paired_nodes*/);
    auto interface_nodes_of_bulk_mesh = pn_bulk_mesh.paired_nodes;

    // Partitioned node vector of the other mesh
    auto const pn_other_mesh = partitionNodesByCoordinateMatch(
        other_mesh_nodes, bulk_nodes_in_aabb, aabb);
    auto const interface_nodes_of_other_mesh = pn_other_mesh.paired_nodes;
    auto const internal_nodes_of_other_mesh = *(pn_other_mesh.non_paired_nodes);

    // Common node id mapping: from the merge mesh to the bulk mesh
    auto cn_id_mapping_m2b = findNodeIdMapping(interface_nodes_of_other_mesh,
                                               interface_nodes_of_bulk_mesh);

    std::unordered_map<std::size_t, std::size_t> other_mesh_node_id_dict;

    for (auto const&& [id_mapping, node] :
         ranges::views::zip(cn_id_mapping_m2b, interface_nodes_of_other_mesh))
    {
        other_mesh_node_id_dict.insert({node->getID(), id_mapping});
    }

    // Create new mesh node vector
    std::vector<MeshLib::Node*> new_node_vector;
    new_node_vector.reserve(bulk_mesh_nodes.size() +
                            internal_nodes_of_other_mesh.size());
    // First, copy the nodes from the bulk mesh
    new_node_vector.insert(new_node_vector.end(), bulk_mesh_nodes.begin(),
                           bulk_mesh_nodes.end());
    // Second, append the inside nodes of the merge mesh
    std::size_t new_node_id = bulk_mesh_nodes.size();
    for (auto node : internal_nodes_of_other_mesh)
    {
        other_mesh_node_id_dict.insert({node->getID(), new_node_id});
        new_node_vector.push_back(node);
        new_node_id++;
    }

    // Create new mesh element vector
    auto const elements_bulk = bulk_mesh.getElements();
    auto const elements_other_mesh = other_mesh.getElements();
    std::vector<MeshLib::Element*> new_element_vector;
    new_element_vector.reserve(elements_bulk.size() +
                               elements_other_mesh.size());
    // First, copy the elements from the bulk mesh
    new_element_vector.insert(new_element_vector.end(), elements_bulk.begin(),
                              elements_bulk.end());
    // Append the elements of the mesh to be merged to the element vector of the
    // bulk mesh:
    for (auto element : elements_other_mesh)
    {
        auto const nn = element->getNumberOfNodes();
        for (unsigned i = 0; i < nn; ++i)
        {
            // `node` cannot be declared as `const` because `element` is not
            // `const`.
            MeshLib::Node* node = element->getNode(i);
            auto const new_id = other_mesh_node_id_dict[node->getID()];
            element->setNode(i, new_node_vector[new_id]);
        }
        new_element_vector.push_back(element);
    }

    // Node IDs are reset in the constructor of MeshLib::Mesh according to their
    // order of appearance (see MeshLib::Mesh::resetNodeIDs()). Since the order
    // of nodes in the bulk mesh remains unchanged, the corresponding segments
    // of property vectors of the original bulk mesh are also preserved.
    auto merged_mesh = std::make_unique<MeshLib::Mesh>(
        "merged_mesh", new_node_vector, new_element_vector);

    MeshLib::Properties const& properties_bulk_mesh = bulk_mesh.getProperties();

    MeshLib::applyToPropertyVectors(
        properties_bulk_mesh,
        [&](auto type, auto const& property)
        {
            return createMergedPropertyVector<decltype(type)>(
                *merged_mesh, initial_value_dict,
                dynamic_cast<MeshLib::PropertyVector<decltype(type)> const*>(
                    property),
                properties_bulk_mesh);
        });

    for (auto node : interface_nodes_of_other_mesh)
    {
        delete node;
    }

    return merged_mesh;
}

}  // namespace MeshToolsLib
