/*!
  \file NodePartitionedMesh.h
  \author Wenqing Wang
  \date   2014.06
  \brief  Definition of mesh class for partitioned mesh (by node) for parallel computing within the
          framework of domain decomposition (DDC).

  \copyright
  Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/

#ifndef NODE_PARTITIONED_MESH_H_
#define NODE_PARTITIONED_MESH_H_

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
            \param n_nghost_elem Number of non-ghost elements, or the start ID of
                                 the entry of ghost element in the element vector.
            \param n_global_base_nodes Number of the base nodes of the global mesh.
            \param n_global_nodes      Number of all nodes of the global mesh.
            \param n_base_nodes        Number of the base nodes.
            \param n_active_base_nodes Number of the active base nodes.
            \param n_active_nodes      Number of all active nodes.
        */
        NodePartitionedMesh(const std::string &name,
                            const std::vector<Node*> &nodes,
                            const std::vector<std::size_t> &glb_node_ids,
                            const std::vector<Element*> &elements,
                            Properties properties,
                            const std::size_t n_nghost_elem,
                            const std::size_t n_global_base_nodes,
                            const std::size_t n_global_nodes,
                            const std::size_t n_base_nodes,
                            const std::size_t n_active_base_nodes,
                            const std::size_t n_active_nodes)
            : Mesh(name, nodes, elements, properties, n_base_nodes),
              _global_node_ids(glb_node_ids), _n_nghost_elem(n_nghost_elem),
              _n_global_base_nodes(n_global_base_nodes),
              _n_global_nodes(n_global_nodes),
              _n_active_base_nodes(n_active_base_nodes),
              _n_active_nodes(n_active_nodes)
        {
        }

        /// Get the number of nodes of the global mesh for linear elements.
        std::size_t getNumberOfGlobalBaseNodes() const
        {
            return _n_global_base_nodes;
        }

        /// Get the number of all nodes of the global mesh.
        std::size_t getNumberOfGlobalNodes() const
        {
            return _n_global_nodes;
        }

        /// Get the global node ID of a node with its local ID.
        std::size_t getGlobalNodeID(const std::size_t node_id) const
        {
            return _global_node_ids[node_id];
        }

        /// Get the number of the active nodes of the partition for linear elements.
        std::size_t getNumberOfActiveBaseNodes() const
        {
            return _n_active_base_nodes;
        }

        /// Get the number of all active nodes of the partition.
        std::size_t getNumberOfActiveNodes() const
        {
            return _n_active_nodes;
        }

        /// Check whether a node with ID of node_id is a ghost node
        bool isGhostNode(const std::size_t node_id) const
        {
            if(node_id < _n_active_base_nodes)
                return false;
            else if(node_id >= _n_base_nodes && node_id < getLargestActiveNodeID() )
                return false;
            else
                return true;
        }

        /// Get the largest ID of active nodes for higher order elements in a partition.
        std::size_t getLargestActiveNodeID() const
        {
            return _n_base_nodes + _n_active_nodes - _n_active_base_nodes;
        }

        /// Get the number of non-ghost elements, or the start entry ID of ghost elements in element vector.
        std::size_t getNumberOfNonGhostElements() const
        {
            return _n_nghost_elem;
        }

        // TODO I guess that is a simplified version of computeSparsityPattern()
        /// Get the maximum number of connected nodes to node.
        std::size_t getMaximumNConnectedNodesToNode() const
        {
            std::vector<Node *>::const_iterator it_max_ncn = std::max_element(
                _nodes.cbegin(), _nodes.cend(),
                [](Node const *const node_a, Node const *const node_b)
                {
                    return (node_a->getConnectedNodes().size() <
                            node_b->getConnectedNodes().size());
                });
            // Return the number of connected nodes +1 for the node itself.
            return (*it_max_ncn)->getConnectedNodes().size() + 1;
        }

    private:
        /// Global IDs of nodes of a partition
        std::vector<std::size_t> _global_node_ids;

        /// Number of non-ghost elements, or the ID of the start entry of ghost elements in _elements vector.
        std::size_t _n_nghost_elem;

        /// Number of the nodes of the global mesh linear interpolations.
        std::size_t _n_global_base_nodes;

        /// Number of all nodes of the global mesh.
        std::size_t _n_global_nodes;

        /// Number of the active nodes for linear interpolations
        std::size_t _n_active_base_nodes;

        /// Number of the all active nodes.
        std::size_t _n_active_nodes;
};

}   // namespace MeshLib

#endif // NODE_PARTITIONED_MESH_H_
