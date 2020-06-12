/*!
  \file
  \author Wenqing Wang
  \date   2014.06
  \brief  Definition of mesh class for partitioned mesh (by node) for parallel
  computing within the framework of domain decomposition (DDC).

  \copyright
  Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/

#pragma once

#include <algorithm>
#include <string>
#include <vector>

#include "Mesh.h"
#include "Node.h"

namespace MeshLib
{
/// A subdomain mesh.
class NodePartitionedMesh : public Mesh
{
public:
    // Copy a global mesh for the case of the thread number is one,
    // i.e the gobal mesh is not partitioned.
    // \param mesh The gobal mesh
    explicit NodePartitionedMesh(const Mesh& mesh)
        : Mesh(mesh),
          global_node_ids_(mesh.getNumberOfNodes()),
          n_global_base_nodes_(mesh.getNumberOfBaseNodes()),
          n_global_nodes_(mesh.getNumberOfNodes()),
          n_active_base_nodes_(mesh.getNumberOfBaseNodes()),
          n_active_nodes_(mesh.getNumberOfNodes()),
          is_single_thread_(true)
    {
        const auto& mesh_nodes = mesh.getNodes();
        for (std::size_t i = 0; i < nodes_.size(); i++)
        {
            global_node_ids_[i] = nodes_[i]->getID();

            // TODO To add copying of the connected nodes (and elements)
            //      in the copy constructor of class Node in order to
            //      drop the following lines.
            auto node = nodes_[i];
            // Copy constructor of Mesh does not copy the connected
            // nodes to node.
            if (node->connected_nodes_.size() == 0)
            {
                std::copy(mesh_nodes[i]->connected_nodes_.begin(),
                          mesh_nodes[i]->connected_nodes_.end(),
                          std::back_inserter(node->connected_nodes_));
            }
        }
    }

    /*!
        \brief Constructor
        \param name          Name assigned to the mesh.
        \param nodes         Vector for nodes, which storage looks like:
                             ||--active base nodes--|--ghost base nodes--|
                              --active extra nodes--|--ghost extra nodes--||
                             (extra nodes: nodes for high order interpolations)
        \param glb_node_ids  Global IDs of nodes of a partition.
        \param elements      Vector for elements. Ghost elements are stored
                             after regular (non-ghost) elements.
        \param properties    Mesh property.
        \param n_global_base_nodes Number of the base nodes of the global mesh.
        \param n_global_nodes      Number of all nodes of the global mesh.
        \param n_base_nodes        Number of the base nodes.
        \param n_active_base_nodes Number of the active base nodes.
        \param n_active_nodes      Number of all active nodes.
    */
    NodePartitionedMesh(const std::string& name,
                        const std::vector<Node*>& nodes,
                        const std::vector<std::size_t>& glb_node_ids,
                        const std::vector<Element*>& elements,
                        Properties properties,
                        const std::size_t n_global_base_nodes,
                        const std::size_t n_global_nodes,
                        const std::size_t n_base_nodes,
                        const std::size_t n_active_base_nodes,
                        const std::size_t n_active_nodes)
        : Mesh(name, nodes, elements, properties, n_base_nodes),
          global_node_ids_(glb_node_ids),
          n_global_base_nodes_(n_global_base_nodes),
          n_global_nodes_(n_global_nodes),
          n_active_base_nodes_(n_active_base_nodes),
          n_active_nodes_(n_active_nodes),
          is_single_thread_(false)
    {
    }

    /// Get the number of nodes of the global mesh for linear elements.
    std::size_t getNumberOfGlobalBaseNodes() const
    {
        return n_global_base_nodes_;
    }

    /// Get the number of all nodes of the global mesh.
    std::size_t getNumberOfGlobalNodes() const { return n_global_nodes_; }
    /// Get the global node ID of a node with its local ID.
    std::size_t getGlobalNodeID(const std::size_t node_id) const
    {
        return global_node_ids_[node_id];
    }

    /// Get the number of the active nodes of the partition for linear elements.
    std::size_t getNumberOfActiveBaseNodes() const
    {
        return n_active_base_nodes_;
    }

    /// Get the number of all active nodes of the partition.
    std::size_t getNumberOfActiveNodes() const { return n_active_nodes_; }
    /// Check whether a node with ID of node_id is a ghost node
    bool isGhostNode(const std::size_t node_id) const
    {
        if (node_id < n_active_base_nodes_)
            return false;
        else if (node_id >= n_base_nodes_ && node_id < getLargestActiveNodeID())
            return false;
        else
            return true;
    }

    /// Get the largest ID of active nodes for higher order elements in a
    /// partition.
    std::size_t getLargestActiveNodeID() const
    {
        return n_base_nodes_ + n_active_nodes_ - n_active_base_nodes_;
    }

    // TODO I guess that is a simplified version of computeSparsityPattern()
    /// Get the maximum number of connected nodes to node.
    std::size_t getMaximumNConnectedNodesToNode() const
    {
        std::vector<Node*>::const_iterator it_max_ncn = std::max_element(
            nodes_.cbegin(), nodes_.cend(),
            [](Node const* const node_a, Node const* const node_b) {
                return (node_a->getConnectedNodes().size() <
                        node_b->getConnectedNodes().size());
            });
        // Return the number of connected nodes +1 for the node itself.
        return (*it_max_ncn)->getConnectedNodes().size() + 1;
    }

    bool isForSingleThread() const { return is_single_thread_; }

private:
    /// Global IDs of nodes of a partition
    std::vector<std::size_t> global_node_ids_;

    /// Number of the nodes of the global mesh linear interpolations.
    std::size_t n_global_base_nodes_;

    /// Number of all nodes of the global mesh.
    std::size_t n_global_nodes_;

    /// Number of the active nodes for linear interpolations
    std::size_t n_active_base_nodes_;

    /// Number of the all active nodes.
    std::size_t n_active_nodes_;

    const bool is_single_thread_;
};

}  // namespace MeshLib
