// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "MeshComponentMap.h"

#include <cstddef>
#include <limits>
#include <numeric>
#include <range/v3/algorithm/count_if.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/join.hpp>
#include <range/v3/view/map.hpp>
#include <range/v3/view/subrange.hpp>
#include <range/v3/view/transform.hpp>
#include <unordered_map>

#include "BaseLib/Error.h"
#include "BaseLib/MPI.h"
#include "MeshLib/MeshSubset.h"
#include "MeshLib/Node.h"

#ifdef USE_PETSC
#include "MeshLib/NodePartitionedMesh.h"
#endif

#ifdef USE_PETSC
namespace
{
/// Partitioned mesh of a component. Checked downcast from MeshLib::Mesh.
MeshLib::NodePartitionedMesh const& partitionedMesh(
    MeshLib::MeshSubset const& c)
{
    auto const* const mesh =
        dynamic_cast<MeshLib::NodePartitionedMesh const*>(&c.getMesh());
    if (mesh == nullptr)
    {
        OGS_FATAL("The mesh '{:s}' of the component is not partitioned.",
                  c.getMesh().getName());
    }
    return *mesh;
}

/// Partitioned meshes of all components in component order.
std::vector<MeshLib::NodePartitionedMesh const*> partitionedMeshes(
    std::vector<MeshLib::MeshSubset> const& components)
{
    std::vector<MeshLib::NodePartitionedMesh const*> meshes;
    meshes.reserve(components.size());
    for (auto const& c : components)
    {
        meshes.push_back(&partitionedMesh(c));
    }
    return meshes;
}
}  // namespace
#endif

namespace
{
/// Returns a range of lines for the given location.
auto linesByLocation(NumLib::detail::ComponentGlobalIndexDict const& dict,
                     MeshLib::Location const& l)
{
    auto const& m = dict.get<NumLib::detail::ByLocation>();
    auto const [first, last] = m.equal_range(NumLib::detail::Line{l});
    return ranges::make_subrange(first, last);
}

/// Returns a range of lines for all given locations in the order of the
/// locations.
auto linesByLocations(NumLib::detail::ComponentGlobalIndexDict const& dict,
                      std::vector<MeshLib::Location> const& ls)
{
    return ls |
           ranges::views::transform([&dict](MeshLib::Location const& l)
                                    { return linesByLocation(dict, l); }) |
           ranges::views::join;
}

}  // namespace

namespace NumLib
{
using namespace detail;

#ifdef USE_PETSC
MeshComponentMap::MeshComponentMap(
    std::vector<MeshLib::MeshSubset> const& components, ComponentOrder order)
{
    BaseLib::MPI::Mpi const mpi;
    createParallelMeshComponentMap(components, order, mpi);
}
#else
MeshComponentMap::MeshComponentMap(
    std::vector<MeshLib::MeshSubset> const& components, ComponentOrder order)
{
    createSerialMeshComponentMap(components, order);
    _num_local_dof = _dict.size();
}
#endif  // end of USE_PETSC

MeshComponentMap MeshComponentMap::getSubset(
    std::vector<MeshLib::MeshSubset> const& bulk_mesh_subsets,
    MeshLib::MeshSubset const& new_mesh_subset) const
{
    {  // Testing first an assumption met later in the code that the meshes for
       // all the bulk_mesh_subsets are equal.
        auto const first_mismatch =
            std::adjacent_find(begin(bulk_mesh_subsets), end(bulk_mesh_subsets),
                               [](auto const& a, auto const& b)
                               { return a.getMeshID() != b.getMeshID(); });
        if (first_mismatch != end(bulk_mesh_subsets))
        {
            OGS_FATAL(
                "Assumption in the MeshComponentMap violated. Expecting all of "
                "mesh ids to be the same, but it is not true for "
                "the mesh '{:s}' with id {:d}.",
                first_mismatch->getMesh().getName(),
                first_mismatch->getMeshID());
        }
    }

    // Mapping of the nodes in the new_mesh_subset to the bulk mesh nodes
    auto bulk_node_ids = [](auto const& mesh)
    {
        auto const* bulk_node_ids_ptr = MeshLib::bulkNodeIDs(mesh);
        if (bulk_node_ids_ptr == nullptr)
        {
            OGS_FATAL(
                "Bulk node ids map expected in the construction of the mesh "
                "subset.");
        }
        return *bulk_node_ids_ptr;
    };
    auto const& bulk_node_ids_map = bulk_node_ids(new_mesh_subset.getMesh());

    // New dictionary for the subset.
    ComponentGlobalIndexDict subset_dict;

    std::size_t const new_mesh_id = new_mesh_subset.getMeshID();
    // Lookup the locations in the current mesh component map and
    // insert the full lines into the new subset dictionary.
    for (auto* const node : new_mesh_subset.getNodes())
    {
        auto const node_id = node->getID();
        bool const is_base_node = MeshLib::isBaseNode(
            *node,
            new_mesh_subset.getMesh().getElementsConnectedToNode(node_id));

        MeshLib::Location const new_location{
            new_mesh_id, MeshLib::MeshItemType::Node, node_id};

        // Assuming the meshes for all the bulk_mesh_subsets are equal.
        MeshLib::Location const bulk_location{
            bulk_mesh_subsets.front().getMeshID(), MeshLib::MeshItemType::Node,
            bulk_node_ids_map[node_id]};

        for (auto component_id : getComponentIDs(bulk_location))
        {
            auto const global_index =
                getGlobalIndex(bulk_location, component_id);
            if (global_index == nop)
            {
                if (is_base_node)
                {
                    OGS_FATAL(
                        "Could not find a global index for global component "
                        "{:d} for the mesh '{:s}', node {:d}, in the "
                        "corresponding bulk mesh '{:s}' and node {:d}. This "
                        "happens because the boundary mesh is larger than the "
                        "definition region of the bulk component, usually "
                        "because the geometry for the boundary condition is "
                        "too large.",
                        component_id,
                        new_mesh_subset.getMesh().getName(),
                        node_id,
                        bulk_mesh_subsets.front().getMesh().getName(),
                        bulk_node_ids_map[node_id]);
                }
                continue;
            }
            subset_dict.insert({new_location, component_id, global_index});
        }
    }

    return MeshComponentMap(subset_dict);
}

void MeshComponentMap::renumberByLocation(GlobalIndexType offset)
{
    GlobalIndexType global_index = offset;

    auto& m = _dict.get<ByLocation>();  // view as sorted by mesh item
    for (auto itr_mesh_item = m.begin(); itr_mesh_item != m.end();
         ++itr_mesh_item)
    {
        Line pos = *itr_mesh_item;
        pos.global_index = global_index++;
        m.replace(itr_mesh_item, pos);
    }
}

std::vector<int> MeshComponentMap::getComponentIDs(const Location& l) const
{
    return linesByLocation(_dict, l) |
           ranges::views::transform([](auto const& e) { return e.comp_id; }) |
           ranges::to<std::vector>;
}

static GlobalIndexType getGlobalIndex(
    NumLib::detail::ComponentGlobalIndexDict const& dict,
    MeshLib::Location const& l,
    int const comp_id)
{
    auto const& m = dict.get<ByLocationAndComponent>();
    auto const itr = m.find(Line(l, comp_id));
    return itr != m.end() ? itr->global_index : MeshComponentMap::nop;
}

GlobalIndexType MeshComponentMap::getGlobalIndex(Location const& l,
                                                 int const comp_id) const
{
    return NumLib::getGlobalIndex(_dict, l, comp_id);
}

std::vector<GlobalIndexType> MeshComponentMap::getGlobalIndices(
    const Location& l) const
{
    return linesByLocation(_dict, l) |
           ranges::views::transform([](auto const& e)
                                    { return e.global_index; }) |
           ranges::to<std::vector>;
}

std::vector<GlobalIndexType> MeshComponentMap::getGlobalIndicesByLocation(
    std::vector<Location> const& ls) const
{
    return linesByLocations(_dict, ls) |
           ranges::views::transform([](auto const& e)
                                    { return e.global_index; }) |
           ranges::to<std::vector>;
}

std::vector<GlobalIndexType> MeshComponentMap::getGlobalIndicesByComponent(
    std::vector<Location> const& ls) const
{
    // vector of (Component, global Index) pairs.
    using CIPair = std::pair<int, GlobalIndexType>;

    std::vector<CIPair> pairs =
        linesByLocations(_dict, ls) |
        ranges::views::transform(
            [](auto const& e) { return CIPair{e.comp_id, e.global_index}; }) |
        ranges::to<std::vector>();

    auto CIPairLess = [](CIPair const& a, CIPair const& b)
    { return a.first < b.first; };

    // Sort by component; stable to keep the location order within each.
    if (!std::is_sorted(pairs.begin(), pairs.end(), CIPairLess))
    {
        std::stable_sort(pairs.begin(), pairs.end(), CIPairLess);
    }

    return pairs | ranges::views::values | ranges::to<std::vector>;
}

#ifdef USE_PETSC
/// Ghost entries store their global index negated. Global index 0 cannot be
/// negated and is stored as -num_global_dof instead.
static GlobalIndexType encodeGhostGlobalIndex(
    GlobalIndexType const global_index, std::size_t const num_global_dof)
{
    return global_index == 0 ? -static_cast<GlobalIndexType>(num_global_dof)
                             : -global_index;
}

/// Inverse of encodeGhostGlobalIndex().
static GlobalIndexType decodeGhostGlobalIndex(
    GlobalIndexType const encoded_global_index,
    std::size_t const num_global_dof)
{
    return -encoded_global_index == static_cast<GlobalIndexType>(num_global_dof)
               ? 0
               : -encoded_global_index;
}
#endif

GlobalIndexType MeshComponentMap::getLocalIndex(
    Location const& l,
    int const comp_id,
    std::size_t const range_begin,
    std::size_t const range_end) const
{
    GlobalIndexType const global_index = getGlobalIndex(l, comp_id);
    // request for index of linear quantities at higher order nodes
    // results in returning nop
    // That index shall not be modified like a usual global index.
    if (global_index == nop)
    {
        return nop;
    }
#ifndef USE_PETSC
    (void)range_begin;
    (void)range_end;
    return global_index;
#else
    if (global_index >= 0)  // non-ghost location.
        return global_index - range_begin;

    //
    // For a ghost location look up the global index in ghost indices.
    //

    GlobalIndexType const real_global_index =
        decodeGhostGlobalIndex(global_index, _num_global_dof);

    // TODO Find in ghost indices is O(n^2/2) for n being the length of
    // _ghosts_indices. Providing an inverted table would be faster.
    auto const ghost_index_it = std::find(
        _ghosts_indices.begin(), _ghosts_indices.end(), real_global_index);
    if (ghost_index_it == _ghosts_indices.end())
    {
        OGS_FATAL("index {:d} not found in ghost_indices", real_global_index);
    }

    // Using std::distance on a std::vector is O(1). As long as _ghost_indices
    // remains of std::vector type, this shall be fast.
    return range_end - range_begin +
           std::distance(_ghosts_indices.begin(), ghost_index_it);

#endif
}

/// Global index of an entry before it is numbered. Never a valid global index
/// and never produced by encodeGhostGlobalIndex() below, so an unfilled ghost
/// entry cannot be mistaken for a filled one. Decoding it yields nop.
static GlobalIndexType constexpr unassigned_ghost_index =
    -std::numeric_limits<GlobalIndexType>::max();

void MeshComponentMap::createSerialMeshComponentMap(
    std::vector<MeshLib::MeshSubset> const& components, ComponentOrder order,
    bool const number_entries)
{
    // construct dict (and here we number global_index by component type)
    GlobalIndexType global_index = 0;
    int comp_id = 0;
    for (auto const& c : components)
    {
        std::size_t const mesh_id = c.getMeshID();
        // mesh items are ordered first by node, cell, ....
        for (auto const node_id : c.getNodes() | MeshLib::views::ids)
        {
            _dict.insert(
                Line(Location(mesh_id, MeshLib::MeshItemType::Node, node_id),
                     comp_id,
                     number_entries ? global_index++ : unassigned_ghost_index));
        }
        comp_id++;
    }

    if (!number_entries)
    {
        return;  // the parallel path numbers all entries in a single pass
    }

    if (order == ComponentOrder::BY_LOCATION)
    {
        renumberByLocation();
    }
}

#ifdef USE_PETSC
/// Numbers the regular (non-ghost) entries [offset, offset + count) in
/// ByLocation order. Ghost entries keep `unassigned_ghost_index` and get their
/// real index from the exchange below.
///
/// The single pass (instead of the serial numbering plus compacting) is
/// needed because the BY_LOCATION renumbering counts ghost entries, too, and
/// would run past this rank's reserved block into the next rank's.
static void numberRegularEntriesByLocation(
    NumLib::detail::ComponentGlobalIndexDict& dict,
    std::vector<MeshLib::NodePartitionedMesh const*> const& component_meshes,
    std::size_t const offset)
{
    GlobalIndexType global_index = offset;
    auto& by_location = dict.get<ByLocation>();
    for (auto itr_mesh_item = by_location.begin();
         itr_mesh_item != by_location.end();
         ++itr_mesh_item)
    {
        if (itr_mesh_item->location.item_type != MeshLib::MeshItemType::Node)
        {
            OGS_FATAL(
                "numberRegularEntriesByLocation: expected a Node location, "
                "got item_type {}.",
                static_cast<int>(itr_mesh_item->location.item_type));
        }
        if (component_meshes[itr_mesh_item->comp_id]->isGhostNode(
                itr_mesh_item->location.item_id))
        {
            continue;  // ghost entry, filled in by the exchange below
        }
        Line l = *itr_mesh_item;  // a copy for later modification
        l.global_index = global_index++;
        by_location.replace(itr_mesh_item, l);
    }
}

static std::vector<std::size_t> numberOfRegularNodesPerComponent(
    std::vector<MeshLib::MeshSubset> const& components)
{
    auto count_regular_nodes = [](MeshLib::MeshSubset const& c) -> std::size_t
    {
        auto const& mesh = partitionedMesh(c);
        // The component is defined on all nodes of the mesh.
        if (mesh.getNumberOfNodes() == c.getNodes().size())
        {
            return mesh.getNumberOfRegularNodes();
        }
        // The component is defined on a subset of the mesh nodes only, e.g. a
        // linear variable on the base nodes of a quadratic mesh. Count the
        // regular (non-ghost) nodes of that subset.
        return ranges::count_if(c.getNodes() | MeshLib::views::ids,
                                [&mesh](auto const n_id)
                                { return !mesh.isGhostNode(n_id); });
    };

    return components | ranges::views::transform(count_regular_nodes) |
           ranges::to<std::vector>;
}

/// Owns a committed MPI datatype and frees it at the end of its lifetime.
class CommittedMpiDatatype
{
public:
    explicit CommittedMpiDatatype(MPI_Datatype const type) : type_{type}
    {
        MPI_Type_commit(&type_);
    }

    CommittedMpiDatatype(CommittedMpiDatatype const&) = delete;
    CommittedMpiDatatype& operator=(CommittedMpiDatatype const&) = delete;
    CommittedMpiDatatype(CommittedMpiDatatype&&) = delete;
    CommittedMpiDatatype& operator=(CommittedMpiDatatype&&) = delete;

    ~CommittedMpiDatatype() { MPI_Type_free(&type_); }

    MPI_Datatype get() const { return type_; }

private:
    MPI_Datatype type_;
};

struct GhostNodeInfo
{
    std::size_t node_id;
    std::size_t global_node_id;
    int comp_id;
};

/// (GhostNodeInfo, global index) pairs grouped by the rank they are for.
using GhostNodeInfosByRank =
    std::vector<std::vector<std::pair<GhostNodeInfo, GlobalIndexType>>>;

static std::vector<GhostNodeInfo> getGhostNodeInfo(
    NumLib::detail::ComponentGlobalIndexDict const& dict,
    std::vector<MeshLib::NodePartitionedMesh const*> const& component_meshes)
{
    std::vector<GhostNodeInfo> ids;
    for (Line const& l : dict.get<ByLocation>())
    {
        if (l.location.item_type != MeshLib::MeshItemType::Node)
        {
            OGS_FATAL(
                "getGhostNodeInfo: expected a Node location, got item_type "
                "{}.",
                static_cast<int>(l.location.item_type));
        }
        auto const& mesh = *component_meshes[l.comp_id];
        if (mesh.isGhostNode(l.location.item_id))
        {
            ids.emplace_back(l.location.item_id,
                             mesh.getGlobalNodeID(l.location.item_id),
                             l.comp_id);
        }
    }
    return ids;
}

/// Mapping of global node ids to local node ids of the regular nodes of a
/// partitioned mesh.
static std::unordered_map<std::size_t, std::size_t> globalToLocalRegularNodeIds(
    MeshLib::NodePartitionedMesh const& mesh)
{
    std::unordered_map<std::size_t, std::size_t> global_to_local;
    global_to_local.reserve(mesh.getNumberOfRegularNodes());
    for (std::size_t node_id = 0; node_id < mesh.getNumberOfRegularNodes();
         ++node_id)
    {
        global_to_local.emplace(mesh.getGlobalNodeID(node_id), node_id);
    }
    return global_to_local;
}

/// Whether on all ranks the global node ids of the regular nodes follow the
/// partition block layout, i.e. getPartitionID() names the owning rank. Not
/// the case for subdomain meshes, see transformMeshToNodePartitionedMesh().
static bool ownerFollowsFromGlobalNodeId(
    MeshLib::NodePartitionedMesh const& mesh, BaseLib::MPI::Mpi const& mpi)
{
    bool follows = true;
    for (std::size_t node_id = 0; node_id < mesh.getNumberOfRegularNodes();
         ++node_id)
    {
        if (mesh.getPartitionID(mesh.getGlobalNodeID(node_id)) !=
            static_cast<std::size_t>(mpi.rank))
        {
            follows = false;
            break;
        }
    }
    return BaseLib::MPI::allOf(follows, mpi);
}

/// Answers a ghost node query on the rank owning the node: the global index
/// of the queried entry here. Returns nop if the node is not a regular node
/// here or outside the component's mesh subset (e.g. a base-node-only
/// component and a non-base node), so the requesting rank can name the d.o.f.
/// that failed to resolve instead of silently keeping the placeholder.
///
/// The global node ids of subdomain meshes do not follow the partition block
/// layout, hence the lookup in global_to_local.
static GlobalIndexType globalIndexOfGhostNodeOnOwner(
    NumLib::detail::ComponentGlobalIndexDict const& dict,
    MeshLib::NodePartitionedMesh const& mesh,
    std::unordered_map<std::size_t, std::size_t> const& global_to_local,
    GhostNodeInfo const& gni)
{
    auto const it = global_to_local.find(gni.global_node_id);
    if (it == global_to_local.end())
    {
        return MeshComponentMap::nop;
    }
    return getGlobalIndex(
        dict,
        MeshLib::Location{mesh.getID(), MeshLib::MeshItemType::Node,
                          it->second},
        gni.comp_id);
}

/// Exchanges the numbers of entries each rank sends to each rank. Returns
/// the receive counts: entry r is the number of entries rank r sends here.
static std::vector<int> exchangeGhostNodeSizes(
    std::vector<int> const& send_counts, BaseLib::MPI::Mpi const& mpi)
{
    std::vector<int> receive_counts(mpi.size);
    MPI_Alltoall(send_counts.data(), 1, MPI_INT, receive_counts.data(), 1,
                 MPI_INT, mpi.communicator);
    return receive_counts;
}

// Custom MPI_Datatype for std::pair<GhostNodeInfo, GlobalIndexType>.
static CommittedMpiDatatype createMpiDatatypeForPair()
{
    constexpr int nitems = 4;
    std::array<int, nitems> const block_lengths = {1, 1, 1, 1};

    using GhostNodeInfoPair = std::pair<GhostNodeInfo, GlobalIndexType>;
    MPI_Aint const ghost_node_info_offset = offsetof(GhostNodeInfoPair, first);
    std::array<MPI_Aint, nitems> const offsets = {
        ghost_node_info_offset + offsetof(GhostNodeInfo, node_id),
        ghost_node_info_offset + offsetof(GhostNodeInfo, global_node_id),
        ghost_node_info_offset + offsetof(GhostNodeInfo, comp_id),
        offsetof(GhostNodeInfoPair, second)};

    // MPIU_INT matches PetscInt (also with --with-64-bit-indices), for which
    // BaseLib::MPI::mpiType() has no mapping.
    static_assert(std::is_same_v<GlobalIndexType, PetscInt>);
    std::array<MPI_Datatype, nitems> const types = {
        BaseLib::MPI::mpiType<std::size_t>(),
        BaseLib::MPI::mpiType<std::size_t>(), MPI_INT, MPIU_INT};

    MPI_Datatype struct_mpi_type;
    MPI_Type_create_struct(nitems, block_lengths.data(), offsets.data(),
                           types.data(), &struct_mpi_type);
    // The struct type's extent may lack the trailing padding of the pair (e.g.
    // for a 32-bit PetscInt), which would misplace all array elements after
    // the first one. Set the extent to the pair's size explicitly.
    MPI_Datatype pair_mpi_type;
    MPI_Type_create_resized(struct_mpi_type, 0, sizeof(GhostNodeInfoPair),
                            &pair_mpi_type);
    MPI_Type_free(&struct_mpi_type);
    return CommittedMpiDatatype{pair_mpi_type};
}

/// Sparse all-to-all (MPI_Isend/MPI_Irecv) of (GhostNodeInfo, global index)
/// pairs: sends to_each_rank[r] to rank r and returns the pairs received from
/// each rank (result[r] came from rank r).
static GhostNodeInfosByRank exchangeGhostInfoPointToPoint(
    GhostNodeInfosByRank const& to_each_rank, BaseLib::MPI::Mpi const& mpi)
{
    auto const pair_mpi_type = createMpiDatatypeForPair();

    std::vector<int> const send_counts =
        to_each_rank |
        ranges::views::transform([](auto const& v)
                                 { return static_cast<int>(v.size()); }) |
        ranges::to<std::vector>;
    std::vector<int> const receive_counts =
        exchangeGhostNodeSizes(send_counts, mpi);

    GhostNodeInfosByRank from_each_rank(mpi.size);

    // Non-blocking communication using MPI_Isend and MPI_Irecv. The grouped
    // vectors are contiguous, so no extra linearization is needed.
    std::vector<MPI_Request> requests;

    // Send data
    for (int dest_rank = 0; dest_rank < mpi.size; ++dest_rank)
    {
        if (send_counts[dest_rank] > 0)
        {  // Only send if there is data to send
            requests.push_back(MPI_Request());
            MPI_Isend(to_each_rank[dest_rank].data(), send_counts[dest_rank],
                      pair_mpi_type.get(), dest_rank, 0, mpi.communicator,
                      &requests.back());
        }
    }

    // Receive data
    for (int src_rank = 0; src_rank < mpi.size; ++src_rank)
    {
        from_each_rank[src_rank].resize(receive_counts[src_rank]);
        if (receive_counts[src_rank] > 0)
        {  // Only receive if there is data to receive
            requests.push_back(MPI_Request());
            MPI_Irecv(from_each_rank[src_rank].data(), receive_counts[src_rank],
                      pair_mpi_type.get(), src_rank, 0, mpi.communicator,
                      &requests.back());
        }
    }

    // Wait for all communications to complete
    MPI_Waitall(static_cast<int>(requests.size()), requests.data(),
                MPI_STATUSES_IGNORE);

    return from_each_rank;
}

void MeshComponentMap::createParallelMeshComponentMap(
    std::vector<MeshLib::MeshSubset> const& components, ComponentOrder order,
    BaseLib::MPI::Mpi const& mpi)
{
    if (order != ComponentOrder::BY_LOCATION)
    {
        // Not allowed in parallel case since this is not suitable to
        // arrange non ghost entries of a partition within
        // a rank in the parallel computing.
        OGS_FATAL(
            "Global index in the system of equations can only be numbered by "
            "the order type of ComponentOrder::BY_LOCATION");
    }

    std::vector<std::size_t> const number_of_regular_nodes_per_component =
        numberOfRegularNodesPerComponent(components);
    _num_local_dof = std::reduce(number_of_regular_nodes_per_component.begin(),
                                 number_of_regular_nodes_per_component.end());

    auto const num_local_dofs_per_rank =
        BaseLib::MPI::allgather(_num_local_dof, mpi);
    _num_global_dof = std::reduce(num_local_dofs_per_rank.begin(),
                                  num_local_dofs_per_rank.end());

    std::vector<std::size_t> global_offsets(num_local_dofs_per_rank.size(), 0);
    std::partial_sum(num_local_dofs_per_rank.begin(),
                     num_local_dofs_per_rank.end() - 1,
                     global_offsets.begin() + 1);

    auto const component_meshes = partitionedMeshes(components);

    // Insert entries only, without numbering them here:
    // numberRegularEntriesByLocation() owns the numbering and would replace
    // the insert-time numbering again.
    createSerialMeshComponentMap(components, order, false);
    numberRegularEntriesByLocation(_dict, component_meshes,
                                   global_offsets[mpi.rank]);

    auto const ghost_node_info = getGhostNodeInfo(_dict, component_meshes);

    std::vector<bool> const owner_from_global_node_id =
        component_meshes |
        ranges::views::transform(
            [&](auto const* mesh)
            { return ownerFollowsFromGlobalNodeId(*mesh, mpi); }) |
        ranges::to<std::vector>;

    // Query/answer round 1: route every ghost query to the rank owning the
    // node instead of replicating all ranks' ghost lists everywhere, and
    // answer the queries routed here. If the owner is unknown, all other
    // ranks are asked.
    GhostNodeInfosByRank queries(mpi.size);
    for (auto const& gni : ghost_node_info)
    {
        if (!owner_from_global_node_id[gni.comp_id])
        {
            for (int rank = 0; rank < mpi.size; ++rank)
            {
                if (rank != mpi.rank)
                {
                    queries[rank].emplace_back(gni, 0);
                }
            }
            continue;
        }
        std::size_t const owner =
            component_meshes[gni.comp_id]->getPartitionID(gni.global_node_id);
        if (owner < static_cast<std::size_t>(mpi.size))
        {  // ids owned by no rank cannot be answered and are reported below
            queries[owner].emplace_back(gni, 0);
        }
    }
    GhostNodeInfosByRank answers = exchangeGhostInfoPointToPoint(queries, mpi);
    std::unordered_map<MeshLib::NodePartitionedMesh const*,
                       std::unordered_map<std::size_t, std::size_t>>
        global_to_local_node_ids;
    for (auto& answers_for_rank : answers)
    {
        for (auto& [gni, global_index] : answers_for_rank)
        {
            auto const& mesh = *component_meshes[gni.comp_id];
            auto [it, inserted] = global_to_local_node_ids.try_emplace(&mesh);
            if (inserted)
            {
                it->second = globalToLocalRegularNodeIds(mesh);
            }
            global_index =
                globalIndexOfGhostNodeOnOwner(_dict, mesh, it->second, gni);
        }
        // Queries sent to all ranks are answered by the owner only.
        std::erase_if(
            answers_for_rank,
            [&](auto const& answer)
            {
                return !owner_from_global_node_id[answer.first.comp_id] &&
                       answer.second == MeshComponentMap::nop;
            });
    }

    // Query/answer round 2: the answers for this rank's ghost nodes.
    auto const receive_data = exchangeGhostInfoPointToPoint(answers, mpi);

    // Every local ghost entry must get exactly one answer back from its owning
    // rank. Checked before the answers are applied, so a missing answer is
    // reported as an exchange failure instead of leaving an unfilled entry.
    std::size_t num_received = 0;
    for (auto const& from_rank : receive_data)
    {
        num_received += from_rank.size();
    }
    std::exception_ptr local_error;
    if (num_received != ghost_node_info.size())
    {
        local_error = std::make_exception_ptr(std::runtime_error(fmt::format(
            "Ghost node global index exchange incomplete on rank {}: {} local "
            "ghost entries but {} answers received from the owning ranks. At "
            "least one ghost node's owning rank could not be determined.",
            mpi.rank, ghost_node_info.size(), num_received)));
    }

    // Replace the unassigned placeholders of the ghost entries (see
    // createSerialMeshComponentMap()) with the negated global indices received
    // from the owning ranks. The checks below are rank-local, but a rank-local
    // OGS_FATAL would leave the other ranks waiting at their next collective
    // call instead of failing the run. The failure is therefore gathered into
    // local_error and re-thrown consistently on every rank by
    // allRanksThrowOrNone().
    auto& index_by_loc_and_comp = _dict.get<ByLocationAndComponent>();
    for (auto const& [gni, gi] : receive_data | ranges::views::join)
    {
        if (local_error)
        {
            break;
        }
        auto const& c = components[gni.comp_id];

        if (gi == MeshComponentMap::nop)
        {
            local_error = std::make_exception_ptr(std::runtime_error(
                fmt::format("Ghost d.o.f. resolution failed on rank {}: node "
                            "{:d} of mesh '{:s}', component {:d} is a ghost "
                            "here, but its owning rank has no global index "
                            "for it.",
                            mpi.rank,
                            gni.node_id,
                            c.getMesh().getName(),
                            gni.comp_id)));
            break;
        }

        std::size_t const mesh_id = c.getMeshID();
        Line const query_line{
            Location(mesh_id, MeshLib::MeshItemType::Node, gni.node_id),
            gni.comp_id, 0};
        auto it = index_by_loc_and_comp.find(query_line);

        if (it == index_by_loc_and_comp.end())
        {
            local_error = std::make_exception_ptr(std::runtime_error(
                fmt::format("Impossible. The entry must have been "
                            "enumerated before. Could not find {}",
                            query_line)));
            break;
        }
        index_by_loc_and_comp.modify(
            it,
            [&](Line& line)
            {
                line.global_index = encodeGhostGlobalIndex(gi, _num_global_dof);
            });
    }
    BaseLib::MPI::allRanksThrowOrNone<std::runtime_error>(local_error);

    // Rebuild _ghosts_indices in ByLocation order (ascending local node id) to
    // match the mesh property layout expected by bulk copyValues().
    for (auto const& line : _dict.get<ByLocation>())
    {
        if (line.global_index >= 0)
        {
            continue;
        }
        _ghosts_indices.push_back(
            decodeGhostGlobalIndex(line.global_index, _num_global_dof));
    }
}
#endif

}  // namespace NumLib
