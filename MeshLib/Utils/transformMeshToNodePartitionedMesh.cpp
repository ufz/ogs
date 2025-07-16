/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "transformMeshToNodePartitionedMesh.h"

#include <mpi.h>

#include <numeric>
#include <range/v3/algorithm/count_if.hpp>
#include <range/v3/numeric.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/indirect.hpp>
#include <range/v3/view/map.hpp>
#include <unordered_map>
#include <vector>

#include "BaseLib/Logging.h"
#include "BaseLib/MPI.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/NodePartitionedMesh.h"
#include "MeshLib/Utils/DuplicateMeshComponents.h"

namespace
{
bool isRegularNode(MeshLib::NodePartitionedMesh const& bulk_mesh,
                   std::span<std::size_t const>
                       local_bulk_node_ids_for_subdomain,
                   std::size_t const subdomain_node_id)
{
    return (!bulk_mesh.isGhostNode(
        local_bulk_node_ids_for_subdomain[subdomain_node_id]));
};
}  // namespace

namespace MeshLib
{
std::pair<std::vector<Node*>, std::vector<Element*>> copyNodesAndElements(
    std::vector<Element*> const& input_elements)
{
    auto elements = cloneElements(input_elements);

    // original node ids to newly created nodes.
    std::unordered_map<std::size_t, Node*> id_node_hash_map;
    id_node_hash_map.reserve(
        elements.size());  // There will be at least one node per element.

    for (auto& e : elements)
    {
        // For each node find a cloned node in map or create if there is none.
        unsigned const n_nodes = e->getNumberOfNodes();
        for (unsigned i = 0; i < n_nodes; ++i)
        {
            Node const* n = e->getNode(i);
            auto const it = id_node_hash_map.find(n->getID());
            if (it == id_node_hash_map.end())
            {
                auto new_node_in_map = id_node_hash_map[n->getID()] =
                    new Node(*n);
                e->setNode(i, new_node_in_map);
            }
            else
            {
                e->setNode(i, it->second);
            }
        }
    }

    std::map<std::size_t, Node*> nodes_map;
    for (auto const& n : id_node_hash_map)
    {
        nodes_map[n.first] = n.second;
    }

    // Copy the unique nodes pointers.
    auto element_nodes =
        nodes_map | ranges::views::values | ranges::to<std::vector>;

    return std::make_pair(element_nodes, elements);
}

unsigned long computeNumberOfRegularNodes(NodePartitionedMesh const* bulk_mesh,
                                          Mesh const* subdomain_mesh)
{
    auto const& subdomain_nodes = subdomain_mesh->getNodes();
    auto const& local_bulk_node_ids_for_subdomain =
        *bulkNodeIDs(*subdomain_mesh);
    unsigned long const number_of_regular_nodes = ranges::count_if(
        subdomain_nodes | MeshLib::views::ids,
        [&](std::size_t const id) {
            return isRegularNode(*bulk_mesh,
                                 {local_bulk_node_ids_for_subdomain}, id);
        });

    DBUG("[{}] number of regular nodes: {}", subdomain_mesh->getName(),
         number_of_regular_nodes);
    return number_of_regular_nodes;
}

std::vector<std::size_t>
computeRegularBaseNodeGlobalNodeIDsOfSubDomainPartition(
    NodePartitionedMesh const* bulk_mesh, Mesh const* subdomain_mesh)
{
    // create/fill global nodes information of the sub-domain partition
    auto const& subdomain_nodes = subdomain_mesh->getNodes();
    auto const& local_bulk_node_ids_for_subdomain =
        *bulkNodeIDs(*subdomain_mesh);

    // global node ids for the sub-domain:
    // | regular base node ids | ghost base node ids | regular higher
    // order node ids | ghost higher order node ids |
    // | partition offset + local ids | regular node ids from other
    // partitions |

    // In order to compute the global node ids for the subdomain mesh nodes the
    // sum of the number of regular nodes of the previous partitions is required
    // first.
    // Count the local regular base nodes to compute offset for other subsequent
    // partitions
    std::size_t const number_of_regular_nodes =
        computeNumberOfRegularNodes(bulk_mesh, subdomain_mesh);

    BaseLib::MPI::Mpi mpi;

    // send own number of regular nodes to all others
    std::vector<std::size_t> const gathered_number_of_regular_nodes =
        BaseLib::MPI::allgather(number_of_regular_nodes, mpi);

    // compute the 'offset' in the global_node_ids
    std::vector<std::size_t> const numbers_of_regular_nodes_at_rank =
        BaseLib::sizesToOffsets(gathered_number_of_regular_nodes);

    // add the offset to the partitioned-owned subdomain
    std::vector<std::size_t> subdomain_global_node_ids;
    subdomain_global_node_ids.reserve(subdomain_nodes.size());
    auto const partition_offset = numbers_of_regular_nodes_at_rank[mpi.rank];
    DBUG("[{}] partition offset: {}", subdomain_mesh->getName(),
         partition_offset);
    // set the global id for the regular base nodes
    for (auto const id : subdomain_nodes | MeshLib::views::ids)
    {
        if (isRegularNode(*bulk_mesh, {local_bulk_node_ids_for_subdomain}, id))
        {
            subdomain_global_node_ids.emplace_back(partition_offset + id);
        }
    }
    return subdomain_global_node_ids;
}

std::vector<std::size_t> computeGhostBaseNodeGlobalNodeIDsOfSubDomainPartition(
    NodePartitionedMesh const* bulk_mesh,
    Mesh const* subdomain_mesh,
    std::vector<std::size_t> const& global_regular_base_node_ids)
{
    BaseLib::MPI::Mpi mpi;

    // count regular nodes that is the offset in the mapping
    auto const local_bulk_node_ids_for_subdomain =
        *bulkNodeIDs(*subdomain_mesh);

    // compute mapping subdomain ghost node id <-> bulk ghost node id
    std::vector<std::size_t> subdomain_node_id_to_bulk_node_id;
    for (auto const id : subdomain_mesh->getNodes() | MeshLib::views::ids)
    {
        if (!isRegularNode(*bulk_mesh, {local_bulk_node_ids_for_subdomain}, id))
        {
            auto const bulk_node_id = local_bulk_node_ids_for_subdomain[id];
            // this puts the global bulk node id at pos:
            // subdomain_node->getID() - number_of_regular_nodes
            subdomain_node_id_to_bulk_node_id.push_back(
                bulk_mesh->getGlobalNodeID(bulk_node_id));
        }
    }

    // Send ids of bulk ghost nodes belonging to the subdomain mesh to all
    // other ranks and at the same time receive from all other ranks.
    std::vector<std::size_t> ghost_node_ids_of_all_ranks;
    std::vector<int> const offsets =
        BaseLib::MPI::allgatherv(std::span{subdomain_node_id_to_bulk_node_id},
                                 ghost_node_ids_of_all_ranks, mpi);

    std::size_t const global_number_of_subdomain_node_id_to_bulk_node_id =
        offsets.back();
    DBUG("[{}] global_number_of_subdomain_node_id_to_bulk_node_id: '{}' ",
         subdomain_mesh->getName(),
         global_number_of_subdomain_node_id_to_bulk_node_id);

    // construct a map for fast search of local bulk node ids
    std::map<std::size_t, std::size_t> global_to_local_bulk_node_ids;
    for (auto const id : bulk_mesh->getNodes() | MeshLib::views::ids)
    {
        if (!bulk_mesh->isGhostNode(id))
        {
            global_to_local_bulk_node_ids[bulk_mesh->getGlobalNodeID(id)] = id;
        }
    }

    // construct a map for fast search of local sub-domain node ids
    std::map<std::size_t, std::size_t> local_subdomain_node_ids;
    for (auto const id : subdomain_mesh->getNodes() | MeshLib::views::ids)
    {
        local_subdomain_node_ids[local_bulk_node_ids_for_subdomain[id]] = id;
    }

    std::vector<std::size_t> subdomain_node_ids_of_all_ranks(
        global_number_of_subdomain_node_id_to_bulk_node_id,
        std::numeric_limits<std::size_t>::max());
    // search in all ranks within the bulk ids for the corresponding id
    for (int rank = 0; rank < mpi.size; ++rank)
    {
        if (rank == mpi.rank)
        {
            continue;
        }
        for (int i = offsets[rank]; i < offsets[rank + 1]; ++i)
        {
            auto it = global_to_local_bulk_node_ids.find(
                ghost_node_ids_of_all_ranks[i]);
            if (it == global_to_local_bulk_node_ids.end())
            {
                continue;
            }
            auto const local_bulk_node_id = it->second;
            subdomain_node_ids_of_all_ranks[i] = global_regular_base_node_ids
                [local_subdomain_node_ids.find(local_bulk_node_id)->second];
            DBUG(
                "[{}] found global subdomain node id: '{}' for global bulk "
                "node id {} ",
                subdomain_mesh->getName(),
                subdomain_node_ids_of_all_ranks[i],
                ghost_node_ids_of_all_ranks[i]);
        }
    }

    // find maximum over all ranks
    BaseLib::MPI::allreduceInplace(subdomain_node_ids_of_all_ranks, MPI_MAX,
                                   mpi);

    std::vector<std::size_t> global_ids_for_subdomain_ghost_nodes(
        subdomain_node_ids_of_all_ranks.begin() + offsets[mpi.rank],
        subdomain_node_ids_of_all_ranks.begin() + offsets[mpi.rank + 1]);
    return global_ids_for_subdomain_ghost_nodes;
}

std::vector<std::size_t> computeNumberOfRegularBaseNodesAtRank(
    Mesh const* subdomain_mesh)
{
    auto const number_of_regular_base_nodes =
        subdomain_mesh->computeNumberOfBaseNodes();

    BaseLib::MPI::Mpi mpi;

    std::vector<std::size_t> const gathered_number_of_regular_base_nodes =
        BaseLib::MPI::allgather(number_of_regular_base_nodes, mpi);

    return BaseLib::sizesToOffsets(gathered_number_of_regular_base_nodes);
}

// similar to the above only with regular higher order nodes
std::vector<std::size_t> computeNumberOfRegularHigherOrderNodesAtRank(
    Mesh const* subdomain_mesh)
{
    BaseLib::MPI::Mpi mpi;

    auto const number_of_regular_base_nodes =
        subdomain_mesh->computeNumberOfBaseNodes();

    auto const number_of_regular_higher_order_nodes =
        subdomain_mesh->getNumberOfNodes() - number_of_regular_base_nodes;
    std::vector<std::size_t> gathered_number_of_regular_higher_order_nodes =
        BaseLib::MPI::allgather(number_of_regular_higher_order_nodes, mpi);

    return BaseLib::sizesToOffsets(
        gathered_number_of_regular_higher_order_nodes);
}

std::vector<std::size_t> computeGlobalNodeIDsOfSubDomainPartition(
    NodePartitionedMesh const* bulk_mesh, Mesh const* subdomain_mesh)
{
    auto global_regular_base_node_ids =
        computeRegularBaseNodeGlobalNodeIDsOfSubDomainPartition(bulk_mesh,
                                                                subdomain_mesh);
    auto const local_ids_for_subdomain_ghost_nodes =
        computeGhostBaseNodeGlobalNodeIDsOfSubDomainPartition(
            bulk_mesh, subdomain_mesh, global_regular_base_node_ids);

    // At the moment higher order bulk meshes aren't handled correctly.
    // If it should become necessary in the future, the necessary things must be
    // implemented at this point.
    /*
    auto const regular_higher_order_node_global_node_ids =
        computeRegularHigherOrderNodeGlobalNodeIDsofSubDomainPartition(
            bulk_mesh, subdomain_mesh);
    auto const ghost_higher_order_node_global_node_ids =
        computeGhostHigherOrderNodeGlobalNodeIDsofSubDomainPartition(
            bulk_mesh, subdomain_mesh);
    */

    std::vector<std::size_t> global_node_ids(
        std::move(global_regular_base_node_ids));
    global_node_ids.insert(global_node_ids.end(),
                           local_ids_for_subdomain_ghost_nodes.begin(),
                           local_ids_for_subdomain_ghost_nodes.end());

    return global_node_ids;
}

unsigned long getNumberOfGlobalNodes(Mesh const* subdomain_mesh)
{
    // sum all nodes over all partitions in number_of_global_nodes
    unsigned long const number_of_global_nodes = BaseLib::MPI::allreduce(
        subdomain_mesh->getNodes().size(), MPI_SUM, BaseLib::MPI::Mpi{});
    DBUG("[{}] number_of_global_nodes: {}'", subdomain_mesh->getName(),
         number_of_global_nodes);
    return number_of_global_nodes;
}

std::unique_ptr<NodePartitionedMesh> transformMeshToNodePartitionedMesh(
    NodePartitionedMesh const* const bulk_mesh,
    Mesh const* const subdomain_mesh)
{
    DBUG("Creating NodePartitionedMesh from '{}'", subdomain_mesh->getName());

    auto const subdomain_global_node_ids =
        computeGlobalNodeIDsOfSubDomainPartition(bulk_mesh, subdomain_mesh);

    // according to comment in struct PartitionedMeshInfo this value
    // is unused
    unsigned long number_of_global_base_nodes = 0;

    unsigned long const number_of_global_nodes =
        getNumberOfGlobalNodes(subdomain_mesh);
    auto numbers_of_regular_base_nodes_at_rank =
        computeNumberOfRegularBaseNodesAtRank(subdomain_mesh);
    auto numbers_of_regular_higher_order_nodes_at_rank =
        computeNumberOfRegularHigherOrderNodesAtRank(subdomain_mesh);

    auto const number_of_regular_nodes =
        computeNumberOfRegularNodes(bulk_mesh, subdomain_mesh);
    DBUG("[{}] number_of_regular_nodes: {}", subdomain_mesh->getName(),
         number_of_regular_nodes);
    auto const [nodes, elements] =
        copyNodesAndElements(subdomain_mesh->getElements());

    return std::make_unique<NodePartitionedMesh>(
        subdomain_mesh->getName(), nodes, subdomain_global_node_ids, elements,
        subdomain_mesh->getProperties(), number_of_global_base_nodes,
        number_of_global_nodes, number_of_regular_nodes,
        std::move(numbers_of_regular_base_nodes_at_rank),
        std::move(numbers_of_regular_higher_order_nodes_at_rank));
}

}  // namespace MeshLib
