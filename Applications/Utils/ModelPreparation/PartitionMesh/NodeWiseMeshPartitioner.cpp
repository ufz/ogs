/*!
  \file
  \date   2016.05

  \brief  Define the members of class NodeWiseMeshPartitioner

  \copyright
  Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/

#include "NodeWiseMeshPartitioner.h"

#include <limits>
#include <numeric>
#include <unordered_map>

#include "BaseLib/Error.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/Stream.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"

namespace ApplicationUtils
{
struct NodeStruct
{
    NodeStruct(NodeWiseMeshPartitioner::IntegerType const id_,
               double const x_,
               double const y_,
               double const z_)
        : id(id_), x(x_), y(y_), z(z_)
    {
    }

    NodeWiseMeshPartitioner::IntegerType id;
    double x;
    double y;
    double z;
};

std::size_t Partition::numberOfMeshItems(
    MeshLib::MeshItemType const item_type) const
{
    if (item_type == MeshLib::MeshItemType::Node)
    {
        return nodes.size();
    }

    if (item_type == MeshLib::MeshItemType::Cell)
    {
        return regular_elements.size() + ghost_elements.size();
    }
    OGS_FATAL("Mesh items other than nodes and cells are not supported.");
}

std::ostream& Partition::writeNodes(
    std::ostream& os, std::vector<std::size_t> const& global_node_ids) const
{
    std::vector<NodeStruct> nodes_buffer;
    nodes_buffer.reserve(nodes.size());

    for (const auto* node : nodes)
    {
        double const* coords = node->getCoords();
        nodes_buffer.emplace_back(global_node_ids[node->getID()], coords[0],
                                  coords[1], coords[2]);
    }
    return os.write(reinterpret_cast<const char*>(nodes_buffer.data()),
                    sizeof(NodeStruct) * nodes_buffer.size());
}

/// Calculate the total number of integer variables of an element vector. Each
/// element has three integer variables for element ID, element type, number of
/// nodes of the element. Therefore the total number of the integers in
/// \c elements is 3 * elements.size() + sum (number of nodes of each element).
NodeWiseMeshPartitioner::IntegerType getNumberOfIntegerVariablesOfElements(
    std::vector<const MeshLib::Element*> const& elements)
{
    return 3 * elements.size() +
           std::accumulate(begin(elements), end(elements), 0,
                           [](auto const nnodes, auto const* e) {
                               return nnodes + e->getNumberOfNodes();
                           });
}

std::ostream& Partition::writeConfig(std::ostream& os) const
{
    long const data[] = {
        static_cast<long>(nodes.size()),
        static_cast<long>(number_of_base_nodes),
        static_cast<long>(regular_elements.size()),
        static_cast<long>(ghost_elements.size()),
        static_cast<long>(number_of_regular_base_nodes),
        static_cast<long>(number_of_regular_nodes),
        static_cast<long>(number_of_mesh_base_nodes),
        static_cast<long>(number_of_mesh_all_nodes),
        static_cast<long>(
            getNumberOfIntegerVariablesOfElements(regular_elements)),
        static_cast<long>(
            getNumberOfIntegerVariablesOfElements(ghost_elements)),
    };

    return os.write(reinterpret_cast<const char*>(data), sizeof(data));
}

std::size_t nodeIdBulkMesh(
    MeshLib::Node const& node,
    std::vector<std::size_t> const* node_id_mapping = nullptr)
{
    return node_id_mapping ? (*node_id_mapping)[node.getID()] : node.getID();
}

std::size_t partitionLookup(
    MeshLib::Node const& node,
    std::vector<std::size_t> const& partition_ids,
    std::vector<std::size_t> const* node_id_mapping = nullptr)
{
    auto node_id = [&node_id_mapping](MeshLib::Node const& n) {
        return nodeIdBulkMesh(n, node_id_mapping);
    };

    return partition_ids[node_id(node)];
}

/// 1 copy pointers to nodes belonging to the partition part_id into base nodes
/// vector, and
/// 2 collect non-linear element nodes belonging to the partition part_id in
/// extra nodes vector.
/// If \c node_id_mapping is given, it will be used to map the mesh node ids to
/// other ids; used by boundary meshes, for example.
/// \return a pair of base node and extra nodes.
std::pair<std::vector<MeshLib::Node*>, std::vector<MeshLib::Node*>>
findRegularNodesInPartition(
    std::size_t const part_id,
    const bool is_mixed_high_order_linear_elems,
    std::vector<MeshLib::Node*> const& nodes,
    std::vector<std::size_t> const& partition_ids,
    std::vector<std::size_t> const* node_id_mapping = nullptr)
{
    // Find nodes belonging to a given partition id.
    std::vector<MeshLib::Node*> partition_nodes;
    copy_if(begin(nodes), end(nodes), std::back_inserter(partition_nodes),
            [&](auto const& n) {
                return partitionLookup(*n, partition_ids, node_id_mapping) ==
                       part_id;
            });

    // Space for resulting vectors.
    std::vector<MeshLib::Node*> base_nodes;
    base_nodes.reserve(partition_nodes.size() /
                       2);  // if linear mesh, then one reallocation, no realloc
                            // for higher order elements meshes.
    std::vector<MeshLib::Node*> higher_order_nodes;
    higher_order_nodes.reserve(
        partition_nodes.size() /
        2);  // if linear mesh, then wasted space, good estimate for quadratic
             // order mesh, and realloc needed for higher order element meshes.

    // Split the nodes into base nodes and extra nodes.
    std::partition_copy(
        begin(partition_nodes), end(partition_nodes),
        std::back_inserter(base_nodes), std::back_inserter(higher_order_nodes),
        [&](MeshLib::Node* const n)
        { return !is_mixed_high_order_linear_elems || isBaseNode(*n); });

    return {base_nodes, higher_order_nodes};
}

std::ptrdiff_t numberOfRegularNodes(
    MeshLib::Element const& e, std::size_t const part_id,
    std::vector<std::size_t> const& partition_ids,
    std::vector<std::size_t> const* node_id_mapping = nullptr)
{
    return std::count_if(e.getNodes(), e.getNodes() + e.getNumberOfNodes(),
                         [&](MeshLib::Node* const n) {
                             return partitionLookup(*n, partition_ids,
                                                    node_id_mapping) == part_id;
                         });
}

/// 1 find elements belonging to the partition part_id:
/// fills vector partition.regular_elements
/// 2 find ghost elements belonging to the partition part_id
/// fills vector partition.ghost_elements
std::tuple<std::vector<MeshLib::Element const*>,
           std::vector<MeshLib::Element const*>>
findElementsInPartition(
    std::size_t const part_id,
    std::vector<MeshLib::Element*> const& elements,
    std::vector<std::size_t> const& partition_ids,
    std::vector<std::size_t> const* node_id_mapping = nullptr)
{
    std::vector<MeshLib::Element const*> regular_elements;
    std::vector<MeshLib::Element const*> ghost_elements;

    for (auto elem : elements)
    {
        auto const regular_nodes = numberOfRegularNodes(
            *elem, part_id, partition_ids, node_id_mapping);

        if (regular_nodes == 0)
        {
            continue;
        }

        if (regular_nodes ==
            static_cast<std::ptrdiff_t>(elem->getNumberOfNodes()))
        {
            regular_elements.push_back(elem);
        }
        else
        {
            ghost_elements.push_back(elem);
        }
    }
    return std::tuple<std::vector<MeshLib::Element const*>,
                      std::vector<MeshLib::Element const*>>{regular_elements,
                                                            ghost_elements};
}

/// Prerequisite: the ghost elements has to be found (using
/// findElementsInPartition).
/// Finds ghost nodes and non-linear element ghost nodes by walking over
/// ghost elements.
std::tuple<std::vector<MeshLib::Node*>, std::vector<MeshLib::Node*>>
findGhostNodesInPartition(
    std::size_t const part_id,
    const bool is_mixed_high_order_linear_elems,
    std::vector<MeshLib::Node*> const& nodes,
    std::vector<MeshLib::Element const*> const& ghost_elements,
    std::vector<std::size_t> const& partition_ids,
    std::vector<std::size_t> const* node_id_mapping = nullptr)
{
    std::vector<MeshLib::Node*> base_nodes;
    std::vector<MeshLib::Node*> ghost_nodes;

    std::vector<bool> is_ghost_node(nodes.size(), false);
    for (const auto* ghost_elem : ghost_elements)
    {
        for (unsigned i = 0; i < ghost_elem->getNumberOfNodes(); i++)
        {
            auto const& n = ghost_elem->getNode(i);
            if (is_ghost_node[n->getID()])
            {
                continue;
            }

            if (partitionLookup(*n, partition_ids, node_id_mapping) != part_id)
            {
                if (!is_mixed_high_order_linear_elems || isBaseNode(*n))
                {
                    base_nodes.push_back(nodes[n->getID()]);
                }
                else
                {
                    ghost_nodes.push_back(nodes[n->getID()]);
                }
                is_ghost_node[n->getID()] = true;
            }
        }
    }
    return std::tuple<std::vector<MeshLib::Node*>, std::vector<MeshLib::Node*>>{
        base_nodes, ghost_nodes};
}

void NodeWiseMeshPartitioner::processPartition(
    std::size_t const part_id, const bool is_mixed_high_order_linear_elems)
{
    auto& partition = _partitions[part_id];
    std::vector<MeshLib::Node*> higher_order_regular_nodes;
    std::tie(partition.nodes, higher_order_regular_nodes) =
        findRegularNodesInPartition(part_id, is_mixed_high_order_linear_elems,
                                    _mesh->getNodes(), _nodes_partition_ids);

    partition.number_of_regular_base_nodes = partition.nodes.size();
    partition.number_of_regular_nodes = partition.number_of_regular_base_nodes +
                                        higher_order_regular_nodes.size();

    std::tie(partition.regular_elements, partition.ghost_elements) =
        findElementsInPartition(part_id, _mesh->getElements(),
                                _nodes_partition_ids);
    std::vector<MeshLib::Node*> base_ghost_nodes;
    std::vector<MeshLib::Node*> higher_order_ghost_nodes;
    std::tie(base_ghost_nodes, higher_order_ghost_nodes) =
        findGhostNodesInPartition(part_id, is_mixed_high_order_linear_elems,
                                  _mesh->getNodes(), partition.ghost_elements,
                                  _nodes_partition_ids);

    std::copy(begin(base_ghost_nodes), end(base_ghost_nodes),
              std::back_inserter(partition.nodes));

    partition.number_of_base_nodes = partition.nodes.size();

    if (is_mixed_high_order_linear_elems)
    {
        std::copy(begin(higher_order_regular_nodes),
                  end(higher_order_regular_nodes),
                  std::back_inserter(partition.nodes));
        std::copy(begin(higher_order_ghost_nodes),
                  end(higher_order_ghost_nodes),
                  std::back_inserter(partition.nodes));
    }

    // Set the node numbers of base and all mesh nodes.
    partition.number_of_mesh_base_nodes = _mesh->getNumberOfBaseNodes();
    partition.number_of_mesh_all_nodes = _mesh->getNumberOfNodes();
}

/// Copies the properties from global property vector \c pv to the
/// partition-local one \c partitioned_pv.
template <typename T>
std::size_t copyNodePropertyVectorValues(
    Partition const& p,
    std::size_t const offset,
    MeshLib::PropertyVector<T> const& pv,
    MeshLib::PropertyVector<T>& partitioned_pv)
{
    auto const& nodes = p.nodes;
    auto const nnodes = nodes.size();
    auto const n_components = pv.getNumberOfGlobalComponents();
    for (std::size_t i = 0; i < nnodes; ++i)
    {
        const auto global_id = nodes[i]->getID();
        std::copy_n(&pv[n_components * global_id], n_components,
                    &partitioned_pv[offset + n_components * i]);
    }
    return n_components * nnodes;
}

/// Copies the properties from global property vector \c pv to the
/// partition-local one \c partitioned_pv. Regular elements' and ghost elements'
/// values are copied.
template <typename T>
std::size_t copyCellPropertyVectorValues(
    Partition const& p,
    std::size_t const offset,
    MeshLib::PropertyVector<T> const& pv,
    MeshLib::PropertyVector<T>& partitioned_pv)
{
    std::size_t const n_regular(p.regular_elements.size());
    auto const n_components = pv.getNumberOfGlobalComponents();
    for (std::size_t i = 0; i < n_regular; ++i)
    {
        const auto id = p.regular_elements[i]->getID();
        std::copy_n(&pv[n_components * id], n_components,
                    &partitioned_pv[offset + n_components * i]);
    }

    std::size_t const n_ghost(p.ghost_elements.size());
    for (std::size_t i = 0; i < n_ghost; ++i)
    {
        const auto id = p.ghost_elements[i]->getID();
        std::copy_n(&pv[n_components * id], n_components,
                    &partitioned_pv[offset + n_components * (n_regular + i)]);
    }
    return n_components * (n_regular + n_ghost);
}

template <typename T>
bool copyPropertyVector(
    MeshLib::Properties& partitioned_properties,
    std::vector<Partition> const& partitions,
    MeshLib::PropertyVector<T> const* const pv,
    std::map<MeshLib::MeshItemType, std::size_t> const& total_number_of_tuples)
{
    if (pv == nullptr)
    {
        return false;
    }
    auto const item_type = pv->getMeshItemType();

    if (item_type == MeshLib::MeshItemType::IntegrationPoint)
    {
        return true;  // Skip integration point data. Requires parsing of json
                      // for the integration point data. Return true, because
                      // the property was "successfully" parsed.
    }

    auto partitioned_pv = partitioned_properties.createNewPropertyVector<T>(
        pv->getPropertyName(), pv->getMeshItemType(),
        pv->getNumberOfGlobalComponents());
    partitioned_pv->resize(total_number_of_tuples.at(item_type) *
                           pv->getNumberOfGlobalComponents());

    auto copy_property_vector_values = [&](Partition const& p,
                                           std::size_t offset) {
        if (item_type == MeshLib::MeshItemType::Node)
        {
            return copyNodePropertyVectorValues(p, offset, *pv,
                                                *partitioned_pv);
        }
        if (item_type == MeshLib::MeshItemType::Cell)
        {
            return copyCellPropertyVectorValues(p, offset, *pv,
                                                *partitioned_pv);
        }
        OGS_FATAL(
            "Copying of property vector values for mesh item type {:s} is not "
            "implemented.",
            item_type);
    };

    std::size_t position_offset(0);
    for (auto p : partitions)
    {
        position_offset += copy_property_vector_values(p, position_offset);
    }
    return true;
}

/// Applies a function of the form f(type, name) -> bool for each of the
/// properties names.
/// The type argument is used to call f<decltype(type)>(name).
/// At least one of the functions must return the 'true' value, but at most one
/// is executed.
template <typename Function>
void applyToPropertyVectors(MeshLib::Properties const& properties, Function f)
{
    for (auto [name, property] : properties)
    {
        // Open question, why is the 'unsigned long' case not compiling giving
        // an error "expected '(' for function-style cast or type construction"
        // with clang-7, and "error C4576: a parenthesized type followed by an
        // initializer list is a non-standard explicit type conversion syntax"
        // with MSVC-15.
        bool success = f(double{}, property) || f(float{}, property) ||
                       f(int{}, property) || f(long{}, property) ||
                       f(unsigned{}, property) || f(long{}, property) ||
                       f(static_cast<unsigned long>(0), property) ||
                       f(std::size_t{}, property) || f(char{}, property) ||
                       f(static_cast<unsigned char>(0), property);
        if (!success)
        {
            OGS_FATAL("Could not apply function to PropertyVector '{:s}'.",
                      property->getPropertyName());
        }
    }
}

void addVtkGhostTypeProperty(MeshLib::Properties& partitioned_properties,
                             std::vector<Partition> const& partitions,
                             std::size_t const total_number_of_cells)
{
    auto* vtk_ghost_type =
        partitioned_properties.createNewPropertyVector<unsigned char>(
            "vtkGhostType", MeshLib::MeshItemType::Cell);
    if (vtk_ghost_type == nullptr)
    {
        OGS_FATAL("Could not create vtkGhostType cell data array.");
    }

    vtk_ghost_type->resize(total_number_of_cells);
    std::size_t offset = 0;
    for (auto const& partition : partitions)
    {
        offset += partition.regular_elements.size();
        for (std::size_t i = 0; i < partition.ghost_elements.size(); ++i)
        {
            if (partition.duplicate_ghost_cell[i])
            {
                (*vtk_ghost_type)[offset + i] |=
                    vtkDataSetAttributes::DUPLICATECELL;
            }
        }
        offset += partition.ghost_elements.size();
    }
}

/// Partition existing properties and add vtkGhostType cell data array property.
MeshLib::Properties partitionProperties(
    MeshLib::Properties const& properties,
    std::vector<Partition> const& partitions)
{
    using namespace MeshLib;

    Properties partitioned_properties;

    auto count_tuples = [&](MeshItemType const mesh_item_type) {
        return std::accumulate(begin(partitions), end(partitions), 0,
                               [&](std::size_t const sum, Partition const& p) {
                                   return sum +
                                          p.numberOfMeshItems(mesh_item_type);
                               });
    };
    std::map<MeshItemType, std::size_t> const total_number_of_tuples = {
        {MeshItemType::Cell, count_tuples(MeshItemType::Cell)},
        {MeshItemType::Node, count_tuples(MeshItemType::Node)}};

    DBUG(
        "total number of tuples after partitioning defined for cells is {:d} "
        "and for nodes {:d}.",
        total_number_of_tuples.at(MeshItemType::Cell),
        total_number_of_tuples.at(MeshItemType::Node));

    // 1 create new PV
    // 2 resize the PV with total_number_of_tuples
    // 3 copy the values according to the partition info
    applyToPropertyVectors(properties, [&](auto type, auto const property) {
        return copyPropertyVector<decltype(type)>(
            partitioned_properties, partitions,
            dynamic_cast<PropertyVector<decltype(type)> const*>(property),
            total_number_of_tuples);
    });

    addVtkGhostTypeProperty(partitioned_properties,
                            partitions,
                            total_number_of_tuples.at(MeshItemType::Cell));

    return partitioned_properties;
}

void markDuplicateGhostCells(MeshLib::Mesh const& mesh,
                             std::vector<Partition>& partitions)
{
    std::vector<bool> cell_visited(mesh.getElements().size(), false);

    for (auto& partition : partitions)
    {
        partition.duplicate_ghost_cell.resize(partition.ghost_elements.size(),
                                              true);

        for (std::size_t i = 0; i < partition.ghost_elements.size(); i++)
        {
            const auto& ghost_element = *partition.ghost_elements[i];
            if (!cell_visited[ghost_element.getID()])
            {
                cell_visited[ghost_element.getID()] = true;
                partition.duplicate_ghost_cell[i] = false;
            }
        }
    }
}

void NodeWiseMeshPartitioner::partitionByMETIS(
    const bool is_mixed_high_order_linear_elems)
{
    for (std::size_t part_id = 0; part_id < _partitions.size(); part_id++)
    {
        INFO("Processing partition: {:d}", part_id);
        processPartition(part_id, is_mixed_high_order_linear_elems);
    }

    markDuplicateGhostCells(*_mesh, _partitions);

    renumberNodeIndices(is_mixed_high_order_linear_elems);

    _partitioned_properties =
        partitionProperties(_mesh->getProperties(), _partitions);
}

void NodeWiseMeshPartitioner::renumberBulkNodeIdsProperty(
    MeshLib::PropertyVector<std::size_t>* const bulk_node_ids_pv,
    std::vector<Partition> const& local_partitions) const
{
    if (bulk_node_ids_pv == nullptr)
    {
        return;
    }

    auto& bulk_node_ids = *bulk_node_ids_pv;

    std::size_t offset = 0;  // offset in property vector for current partition

    assert(_partitions.size() == local_partitions.size());
    int const n_partitions = static_cast<int>(_partitions.size());
    for (int partition_id = 0; partition_id < n_partitions; ++partition_id)
    {
        auto const& bulk_partition = _partitions[partition_id];
        auto const& local_partition = local_partitions[partition_id];

        // Create global-to-local node id mapping for the bulk partition.
        auto const& bulk_nodes = bulk_partition.nodes;
        auto const n_bulk_nodes = bulk_nodes.size();
        std::map<std::size_t, std::size_t> global_to_local;
        for (std::size_t local_node_id = 0; local_node_id < n_bulk_nodes;
             ++local_node_id)
        {
            global_to_local[bulk_nodes[local_node_id]->getID()] = local_node_id;
        }

        auto const& local_nodes = local_partition.nodes;
        auto const n_local_nodes = local_nodes.size();
        for (std::size_t local_node_id = 0; local_node_id < n_local_nodes;
             ++local_node_id)
        {
            bulk_node_ids[offset + local_node_id] =
                global_to_local[bulk_node_ids[offset + local_node_id]];
        }
        offset += n_local_nodes;
    }
}

void NodeWiseMeshPartitioner::renumberBulkElementIdsProperty(
    MeshLib::PropertyVector<std::size_t>* const bulk_element_ids_pv,
    std::vector<Partition> const& local_partitions) const
{
    if (bulk_element_ids_pv == nullptr)
    {
        return;
    }

    auto& bulk_element_ids = *bulk_element_ids_pv;

    std::size_t offset = 0;  // offset in property vector for current partition

    assert(_partitions.size() == local_partitions.size());
    int const n_partitions = static_cast<int>(_partitions.size());
    for (int partition_id = 0; partition_id < n_partitions; ++partition_id)
    {
        auto const& bulk_partition = _partitions[partition_id];
        auto const& local_partition = local_partitions[partition_id];

        // Create global-to-local element id mapping for the bulk partition.
        std::map<std::size_t, std::size_t> global_to_local;
        auto map_elements =
            [&global_to_local](
                std::vector<MeshLib::Element const*> const& elements,
                std::size_t const offset) {
                auto const n_elements = elements.size();
                for (std::size_t e = 0; e < n_elements; ++e)
                {
                    global_to_local[elements[e]->getID()] = offset + e;
                }
            };

        map_elements(bulk_partition.regular_elements, 0);
        map_elements(bulk_partition.ghost_elements,
                     bulk_partition.regular_elements.size());

        // Renumber the local bulk_element_ids map.
        auto renumber_elements =
            [&bulk_element_ids, &global_to_local](
                std::vector<MeshLib::Element const*> const& elements,
                std::size_t const offset) {
                auto const n_elements = elements.size();
                for (std::size_t e = 0; e < n_elements; ++e)
                {
                    bulk_element_ids[offset + e] =
                        global_to_local[bulk_element_ids[offset + e]];
                }
                return n_elements;
            };

        offset += renumber_elements(local_partition.regular_elements, offset);
        offset += renumber_elements(local_partition.ghost_elements, offset);
    }
}

std::vector<Partition> NodeWiseMeshPartitioner::partitionOtherMesh(
    MeshLib::Mesh const& mesh,
    bool const is_mixed_high_order_linear_elems) const
{
    auto const& bulk_node_ids =
        mesh.getProperties().getPropertyVector<std::size_t>(
            "bulk_node_ids", MeshLib::MeshItemType::Node, 1);

    std::vector<Partition> partitions(_partitions.size());
    for (std::size_t part_id = 0; part_id < _partitions.size(); part_id++)
    {
        auto& partition = partitions[part_id];
        INFO("Processing partition: {:d}", part_id);
        // Set the node numbers of base and all mesh nodes.
        partition.number_of_mesh_base_nodes = mesh.getNumberOfBaseNodes();
        partition.number_of_mesh_all_nodes = mesh.getNumberOfNodes();

        std::vector<MeshLib::Node*> higher_order_regular_nodes;
        std::tie(partition.nodes, higher_order_regular_nodes) =
            findRegularNodesInPartition(
                part_id, is_mixed_high_order_linear_elems, mesh.getNodes(),
                _nodes_partition_ids, bulk_node_ids);

        partition.number_of_regular_base_nodes = partition.nodes.size();
        partition.number_of_regular_nodes =
            partition.number_of_regular_base_nodes +
            higher_order_regular_nodes.size();

        std::tie(partition.regular_elements, partition.ghost_elements) =
            findElementsInPartition(part_id, mesh.getElements(),
                                    _nodes_partition_ids, bulk_node_ids);

        std::vector<MeshLib::Node*> base_ghost_nodes;
        std::vector<MeshLib::Node*> higher_order_ghost_nodes;
        std::tie(base_ghost_nodes, higher_order_ghost_nodes) =
            findGhostNodesInPartition(part_id, is_mixed_high_order_linear_elems,
                                      mesh.getNodes(), partition.ghost_elements,
                                      _nodes_partition_ids, bulk_node_ids);

        std::copy(begin(base_ghost_nodes), end(base_ghost_nodes),
                  std::back_inserter(partition.nodes));

        partition.number_of_base_nodes = partition.nodes.size();

        if (is_mixed_high_order_linear_elems)
        {
            std::copy(begin(higher_order_regular_nodes),
                      end(higher_order_regular_nodes),
                      std::back_inserter(partition.nodes));
            std::copy(begin(higher_order_ghost_nodes),
                      end(higher_order_ghost_nodes),
                      std::back_inserter(partition.nodes));
        }
    }

    markDuplicateGhostCells(mesh, partitions);
    return partitions;
}

void NodeWiseMeshPartitioner::renumberNodeIndices(
    const bool is_mixed_high_order_linear_elems)
{
    std::size_t node_global_id_offset = 0;
    // Renumber the global indices.
    // -- Base nodes
    for (auto& partition : _partitions)
    {
        for (std::size_t i = 0; i < partition.number_of_regular_base_nodes; i++)
        {
            _nodes_global_ids[partition.nodes[i]->getID()] =
                node_global_id_offset;
            node_global_id_offset++;
        }
    }

    if (!is_mixed_high_order_linear_elems)
    {
        return;
    }

    // -- Nodes for high order elements.
    for (auto& partition : _partitions)
    {
        const std::size_t end_id = partition.number_of_base_nodes +
                                   partition.number_of_regular_nodes -
                                   partition.number_of_regular_base_nodes;
        for (std::size_t i = partition.number_of_base_nodes; i < end_id; i++)
        {
            _nodes_global_ids[partition.nodes[i]->getID()] =
                node_global_id_offset;
            node_global_id_offset++;
        }
    }
}

template <typename T>
void writePropertyVectorValues(std::ostream& os,
                               MeshLib::PropertyVector<T> const& pv)
{
    os.write(reinterpret_cast<const char*>(pv.data()), pv.size() * sizeof(T));
}

template <typename T>
bool writePropertyVector(MeshLib::PropertyVector<T> const* const pv,
                         MeshLib::MeshItemType const mesh_item_type,
                         std::ostream& out_val, std::ostream& out_meta)
{
    if (pv == nullptr)
    {
        return false;
    }
    // skip property of different mesh item type. Return true, because this
    // operation was successful.
    if (pv->getMeshItemType() != mesh_item_type)
    {
        return true;
    }

    MeshLib::IO::PropertyVectorMetaData pvmd;
    pvmd.property_name = pv->getPropertyName();
    pvmd.fillPropertyVectorMetaDataTypeInfo<T>();
    pvmd.number_of_components = pv->getNumberOfGlobalComponents();
    pvmd.number_of_tuples = pv->getNumberOfTuples();
    writePropertyVectorValues(out_val, *pv);
    MeshLib::IO::writePropertyVectorMetaData(out_meta, pvmd);
    return true;
}

void writeProperties(const std::string& file_name_base,
                     MeshLib::Properties const& partitioned_properties,
                     std::vector<Partition> const& partitions,
                     MeshLib::MeshItemType const mesh_item_type)
{
    auto const number_of_properties =
        partitioned_properties.size(mesh_item_type);
    if (number_of_properties == 0)
    {
        return;
    }

    auto const file_name_infix = toString(mesh_item_type);

    auto const file_name_cfg = file_name_base + "_partitioned_" +
                               file_name_infix + "_properties_cfg" +
                               std::to_string(partitions.size()) + ".bin";
    std::ofstream out(file_name_cfg, std::ios::binary);
    if (!out)
    {
        OGS_FATAL("Could not open file '{:s}' for output.", file_name_cfg);
    }

    auto const file_name_val = file_name_base + "_partitioned_" +
                               file_name_infix + "_properties_val" +
                               std::to_string(partitions.size()) + ".bin";
    std::ofstream out_val(file_name_val, std::ios::binary);
    if (!out_val)
    {
        OGS_FATAL("Could not open file '{:s}' for output.", file_name_val);
    }

    BaseLib::writeValueBinary(out, number_of_properties);

    applyToPropertyVectors(
        partitioned_properties, [&](auto type, auto const& property) {
            return writePropertyVector<decltype(type)>(
                dynamic_cast<MeshLib::PropertyVector<decltype(type)> const*>(
                    property),
                mesh_item_type, out_val, out);
        });

    unsigned long offset = 0;
    for (const auto& partition : partitions)
    {
        MeshLib::IO::PropertyVectorPartitionMetaData pvpmd{
            offset, static_cast<unsigned long>(
                        partition.numberOfMeshItems(mesh_item_type))};
        DBUG(
            "Write meta data for node-based PropertyVector: global offset "
            "{:d}, number of tuples {:d}",
            pvpmd.offset, pvpmd.number_of_tuples);
        MeshLib::IO::writePropertyVectorPartitionMetaData(out, pvpmd);
        offset += pvpmd.number_of_tuples;
    }
}

struct ConfigOffsets
{
    long node_rank_offset;
    long element_rank_offset;
    long ghost_element_rank_offset;

    std::ostream& writeConfig(std::ostream& os) const;
};

std::ostream& ConfigOffsets::writeConfig(std::ostream& os) const
{
    os.write(reinterpret_cast<const char*>(this), sizeof(ConfigOffsets));

    static long reserved = 0;  // Value reserved in the binary format, not used
                               // in the partitioning process.
    return os.write(reinterpret_cast<const char*>(&reserved), sizeof(long));
}

struct PartitionOffsets
{
    long node;
    long regular_elements;
    long ghost_elements;
};

PartitionOffsets computePartitionOffsets(Partition const& partition)
{
    return {static_cast<long>(partition.nodes.size()),
            static_cast<long>(partition.regular_elements.size() +
                              getNumberOfIntegerVariablesOfElements(
                                  partition.regular_elements)),
            static_cast<long>(partition.ghost_elements.size() +
                              getNumberOfIntegerVariablesOfElements(
                                  partition.ghost_elements))};
}

ConfigOffsets incrementConfigOffsets(ConfigOffsets const& oldConfig,
                                     PartitionOffsets const& offsets)
{
    return {
        static_cast<long>(oldConfig.node_rank_offset +
                          offsets.node * sizeof(NodeStruct)),
        // Offset the ending entry of the element integer variables of
        // the non-ghost elements of this partition in the vector of elem_info.
        static_cast<long>(oldConfig.element_rank_offset +
                          offsets.regular_elements * sizeof(long)),

        // Offset the ending entry of the element integer variables of
        // the ghost elements of this partition in the vector of elem_info.
        static_cast<long>(oldConfig.ghost_element_rank_offset +
                          offsets.ghost_elements * sizeof(long))};
}

/// Write the configuration data of the partition data in binary files.
/// \return a pair of vectors for:
///  1. The number of all non-ghost element integer variables for each
///     partition.
///  2. The number of all ghost element integer variables for each partition.
std::tuple<std::vector<long>, std::vector<long>> writeConfigData(
    const std::string& file_name_base, std::vector<Partition> const& partitions)
{
    auto const file_name_cfg = file_name_base + "_partitioned_msh_cfg" +
                               std::to_string(partitions.size()) + ".bin";
    std::ofstream of_bin_cfg(file_name_cfg, std::ios::binary);
    if (!of_bin_cfg)
    {
        OGS_FATAL("Could not open file '{:s}' for output.", file_name_cfg);
    }

    std::vector<long> partitions_element_offsets;
    partitions_element_offsets.reserve(partitions.size());
    std::vector<long> partitions_ghost_element_offsets;
    partitions_ghost_element_offsets.reserve(partitions.size());

    ConfigOffsets config_offsets = {0, 0, 0};  // 0 for first partition.
    for (const auto& partition : partitions)
    {
        partition.writeConfig(of_bin_cfg);

        config_offsets.writeConfig(of_bin_cfg);
        auto const& partition_offsets = computePartitionOffsets(partition);
        config_offsets =
            incrementConfigOffsets(config_offsets, partition_offsets);

        partitions_element_offsets.push_back(
            partition_offsets.regular_elements);
        partitions_ghost_element_offsets.push_back(
            partition_offsets.ghost_elements);
    }

    return std::make_tuple(partitions_element_offsets,
                           partitions_ghost_element_offsets);
}

/// Get integer variables, which are used to define an element
///
/// \param elem            Element
/// \param local_node_ids  Local node indices of a partition
/// \param elem_info       A vector holds all integer variables of
///                        element definitions
/// \param counter         Recorder of the number of integer variables.
void getElementIntegerVariables(
    const MeshLib::Element& elem,
    const std::unordered_map<std::size_t, long>& local_node_ids,
    std::vector<long>& elem_info,
    long& counter)
{
    constexpr unsigned mat_id =
        0;  // TODO: Material ID to be set from the mesh data
    const long nn = elem.getNumberOfNodes();
    elem_info[counter++] = mat_id;
    elem_info[counter++] = static_cast<long>(elem.getCellType());
    elem_info[counter++] = nn;

    for (long i = 0; i < nn; i++)
    {
        auto const& n = *elem.getNode(i);
        elem_info[counter++] = local_node_ids.at(n.getID());
    }
}

/// Generates a mapping of given node ids to a new local (renumbered) node ids.
std::unordered_map<std::size_t, long> enumerateLocalNodeIds(
    std::vector<MeshLib::Node*> const& nodes)
{
    std::unordered_map<std::size_t, long> local_ids;
    local_ids.reserve(nodes.size());

    long local_node_id = 0;
    for (const auto* node : nodes)
    {
        local_ids[node->getID()] = local_node_id++;
    }
    return local_ids;
}

/// Write the element integer variables of all partitions into binary files.
/// \param file_name_base       The prefix of the file name.
/// \param partitions           Partitions vector.
/// \param regular_element_offsets The numbers of all non-ghost element
///                             integer variables of each partitions.
/// \param ghost_element_offsets  The numbers of all ghost element
void writeElements(std::string const& file_name_base,
                   std::vector<Partition> const& partitions,
                   std::vector<long> const& regular_element_offsets,
                   std::vector<long> const& ghost_element_offsets)
{
    const std::string npartitions_str = std::to_string(partitions.size());

    auto const file_name_ele =
        file_name_base + "_partitioned_msh_ele" + npartitions_str + ".bin";
    std::ofstream element_info_os(file_name_ele, std::ios::binary);
    if (!element_info_os)
    {
        OGS_FATAL("Could not open file '{:s}' for output.", file_name_ele);
    }

    auto const file_name_ele_g =
        file_name_base + "_partitioned_msh_ele_g" + npartitions_str + ".bin";
    std::ofstream ghost_element_info_os(file_name_ele_g, std::ios::binary);
    if (!ghost_element_info_os)
    {
        OGS_FATAL("Could not open file '{:s}' for output.", file_name_ele_g);
    }

    for (std::size_t i = 0; i < partitions.size(); i++)
    {
        const auto& partition = partitions[i];
        auto const local_node_ids = enumerateLocalNodeIds(partition.nodes);

        // Vector containing the offsets of the regular elements of this
        // partition
        std::vector<long> ele_info(regular_element_offsets[i]);

        auto writeElementData =
            [&local_node_ids](
                std::vector<MeshLib::Element const*> const& elements,
                long const element_offsets,
                std::ofstream& output_stream) {
                long counter = elements.size();
                std::vector<long> ele_info(element_offsets);

                for (std::size_t j = 0; j < elements.size(); j++)
                {
                    const auto* elem = elements[j];
                    ele_info[j] = counter;
                    getElementIntegerVariables(*elem, local_node_ids, ele_info,
                                               counter);
                }
                // Write vector data of regular elements
                output_stream.write(
                    reinterpret_cast<const char*>(ele_info.data()),
                    ele_info.size() * sizeof(long));
            };

        // regular elements.
        writeElementData(partition.regular_elements, regular_element_offsets[i],
                         element_info_os);
        // Ghost elements
        writeElementData(partition.ghost_elements, ghost_element_offsets[i],
                         ghost_element_info_os);
    }
}

/// Write the nodes of all partitions into a binary file.
/// \param file_name_base The prefix of the file name.
/// \param partitions the list of partitions
/// \param global_node_ids global numbering of nodes
void writeNodes(const std::string& file_name_base,
                std::vector<Partition> const& partitions,
                std::vector<std::size_t> const& global_node_ids)
{
    auto const file_name = file_name_base + "_partitioned_msh_nod" +
                           std::to_string(partitions.size()) + ".bin";
    std::ofstream os(file_name, std::ios::binary);
    if (!os)
    {
        OGS_FATAL("Could not open file '{:s}' for output.", file_name);
    }

    for (const auto& partition : partitions)
    {
        partition.writeNodes(os, global_node_ids);
    }
}

void NodeWiseMeshPartitioner::write(const std::string& file_name_base)
{
    writeProperties(file_name_base, _partitioned_properties, _partitions,
                    MeshLib::MeshItemType::Node);
    writeProperties(file_name_base, _partitioned_properties, _partitions,
                    MeshLib::MeshItemType::Cell);

    const auto elements_offsets = writeConfigData(file_name_base, _partitions);

    const std::vector<IntegerType>& regular_element_offsets =
        std::get<0>(elements_offsets);
    const std::vector<IntegerType>& ghost_element_offsets =
        std::get<1>(elements_offsets);
    writeElements(file_name_base, _partitions, regular_element_offsets,
                  ghost_element_offsets);

    writeNodes(file_name_base, _partitions, _nodes_global_ids);
}

void NodeWiseMeshPartitioner::writeOtherMesh(
    std::string const& output_filename_base,
    std::vector<Partition> const& partitions,
    MeshLib::Properties const& partitioned_properties) const
{
    writeNodes(output_filename_base, partitions, _nodes_global_ids);

    const auto elem_integers =
        writeConfigData(output_filename_base, partitions);

    const std::vector<IntegerType>& num_elem_integers =
        std::get<0>(elem_integers);
    const std::vector<IntegerType>& num_g_elem_integers =
        std::get<1>(elem_integers);
    writeElements(output_filename_base, partitions, num_elem_integers,
                  num_g_elem_integers);

    writeProperties(output_filename_base, partitioned_properties, partitions,
                    MeshLib::MeshItemType::Node);
    writeProperties(output_filename_base, partitioned_properties, partitions,
                    MeshLib::MeshItemType::Cell);
}
}  // namespace ApplicationUtils
