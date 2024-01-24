/*!
  \file
  \date   2016.05

  \brief  Define the members of class NodeWiseMeshPartitioner

  \copyright
  Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/

#include "NodeWiseMeshPartitioner.h"

#include <limits>
#include <numeric>
#include <range/v3/algorithm/transform.hpp>
#include <range/v3/range/conversion.hpp>
#include <unordered_map>

#include "BaseLib/Error.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/RunTime.h"
#include "MeshLib/Elements/Elements.h"
#include "MeshLib/IO/NodeData.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/Utils/IntegrationPointWriter.h"
#include "MeshToolsLib//IntegrationPointDataTools.h"

namespace ApplicationUtils
{
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

    if (item_type == MeshLib::MeshItemType::IntegrationPoint)
    {
        return number_of_integration_points;
    }
    OGS_FATAL("Mesh items other than nodes and cells are not supported.");
}

std::ostream& Partition::writeNodes(
    std::ostream& os, std::vector<std::size_t> const& global_node_ids) const
{
    std::vector<MeshLib::IO::NodeData> nodes_buffer;
    nodes_buffer.reserve(nodes.size());

    for (const auto* node : nodes)
    {
        double const* coords = node->data();
        nodes_buffer.emplace_back(global_node_ids[node->getID()], coords[0],
                                  coords[1], coords[2]);
    }
    return os.write(reinterpret_cast<const char*>(nodes_buffer.data()),
                    sizeof(MeshLib::IO::NodeData) * nodes_buffer.size());
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
                           [](auto const nnodes, auto const* e)
                           { return nnodes + e->getNumberOfNodes(); });
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

std::size_t partitionLookup(std::size_t const& node_id,
                            std::vector<std::size_t> const& partition_ids,
                            std::vector<std::size_t> const& node_id_mapping)
{
    return partition_ids[node_id_mapping[node_id]];
}

std::pair<std::vector<MeshLib::Node const*>, std::vector<MeshLib::Node const*>>
splitIntoBaseAndHigherOrderNodes(std::vector<MeshLib::Node const*> const& nodes,
                                 MeshLib::Mesh const& mesh)
{
    // Space for resulting vectors.
    std::vector<MeshLib::Node const*> base_nodes;
    // if linear mesh, then one reallocation, no realloc for higher order
    // elements meshes.
    base_nodes.reserve(nodes.size() / 2);
    std::vector<MeshLib::Node const*> higher_order_nodes;
    // if linear mesh, then wasted space, good estimate for quadratic
    // order mesh, and realloc needed for higher order element meshes.
    higher_order_nodes.reserve(nodes.size() / 2);

    // Split the nodes into base nodes and extra nodes.
    std::partition_copy(
        begin(nodes), end(nodes), std::back_inserter(base_nodes),
        std::back_inserter(higher_order_nodes),
        [&](MeshLib::Node const* const n)
        { return isBaseNode(*n, mesh.getElementsConnectedToNode(*n)); });

    return {base_nodes, higher_order_nodes};
}

/// Prerequisite: the ghost elements has to be found
/// Finds ghost nodes and non-linear element ghost nodes by walking over
/// ghost elements.
std::tuple<std::vector<MeshLib::Node*>, std::vector<MeshLib::Node*>>
findGhostNodesInPartition(
    std::size_t const part_id,
    std::vector<MeshLib::Node*> const& nodes,
    std::vector<MeshLib::Element const*> const& ghost_elements,
    std::vector<std::size_t> const& partition_ids,
    MeshLib::Mesh const& mesh,
    std::vector<std::size_t> const& node_id_mapping)
{
    std::vector<MeshLib::Node*> base_ghost_nodes;
    std::vector<MeshLib::Node*> higher_order_ghost_nodes;

    std::vector<bool> is_ghost_node(nodes.size(), false);
    for (const auto* ghost_elem : ghost_elements)
    {
        for (unsigned i = 0; i < ghost_elem->getNumberOfNodes(); i++)
        {
            auto const& n = ghost_elem->getNode(i);
            auto const node_id = n->getID();
            if (is_ghost_node[node_id])
            {
                continue;
            }

            if (partitionLookup(node_id, partition_ids, node_id_mapping) !=
                part_id)
            {
                if (isBaseNode(*n, mesh.getElementsConnectedToNode(*n)))
                {
                    base_ghost_nodes.push_back(nodes[node_id]);
                }
                else
                {
                    higher_order_ghost_nodes.push_back(nodes[node_id]);
                }
                is_ghost_node[node_id] = true;
            }
        }
    }
    return std::tuple<std::vector<MeshLib::Node*>, std::vector<MeshLib::Node*>>{
        base_ghost_nodes, higher_order_ghost_nodes};
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

/// Copies the data from the property vector \c pv belonging to the given
/// Partition \c p to the property vector \c partitioned_pv.
/// \c partitioned_pv is ordered by partition. Regular elements' and ghost
/// elements' values are copied.
template <typename T>
std::size_t copyFieldPropertyDataToPartitions(
    MeshLib::Properties const& properties,
    Partition const& p,
    std::size_t const id_offset_partition,
    std::vector<std::size_t> const& element_ip_data_offsets,
    MeshLib::PropertyVector<T> const& pv,
    MeshLib::PropertyVector<T>& partitioned_pv)
{
    // Special field data such as OGS_VERSION, IntegrationPointMetaData,
    // etc., which are not "real" integration points, are copied "as is"
    // (i.e. fully) for every partition.
    if (pv.getPropertyName().find("_ip") == std::string::npos)
    {
        std::copy_n(&pv[0], pv.size(), &partitioned_pv[id_offset_partition]);
        return pv.size();
    }

    auto const n_components = pv.getNumberOfGlobalComponents();

    std::size_t id_offset = 0;

    auto copyFieldData =
        [&](std::vector<const MeshLib::Element*> const& elements)
    {
        auto const ip_meta_data = MeshLib::getIntegrationPointMetaData(
            properties, pv.getPropertyName());

        for (auto const element : elements)
        {
            int const number_of_element_field_data =
                MeshToolsLib::getNumberOfElementIntegrationPoints(ip_meta_data,
                                                                  *element) *
                n_components;
            // The original element ID is not changed.
            auto const element_id = element->getID();
            int const begin_pos = element_ip_data_offsets[element_id];
            int const end_pos = element_ip_data_offsets[element_id + 1];

            std::copy(pv.begin() + begin_pos, pv.begin() + end_pos,
                      &partitioned_pv[id_offset + id_offset_partition]);
            id_offset += number_of_element_field_data;
        }
    };

    copyFieldData(p.regular_elements);
    copyFieldData(p.ghost_elements);

    return id_offset;
}

void setIntegrationPointNumberOfPartition(MeshLib::Properties const& properties,
                                          std::vector<Partition>& partitions)
{
    for (auto const& [name, property] : properties)
    {
        auto const item_type = property->getMeshItemType();

        if (item_type != MeshLib::MeshItemType::IntegrationPoint)
        {
            continue;
        }

        // For special field data such as OGS_VERSION, IntegrationPointMetaData,
        // etc., which are not "real" integration points:
        if (property->getPropertyName().find("_ip") == std::string::npos)
        {
            continue;
        }

        std::string const property_name = property->getPropertyName();
        auto countIntegrationPoints =
            [&](std::vector<const MeshLib::Element*> const& elements)
        {
            auto const ip_meta_data =
                MeshLib::getIntegrationPointMetaData(properties, property_name);
            std::size_t counter = 0;
            for (auto const element : elements)
            {
                int const number_of_integration_points =
                    MeshToolsLib::getNumberOfElementIntegrationPoints(
                        ip_meta_data, *element);
                counter += number_of_integration_points;
            }
            return counter;
        };

        for (auto& p : partitions)
        {
            p.number_of_integration_points =
                countIntegrationPoints(p.regular_elements) +
                countIntegrationPoints(p.ghost_elements);
        }
        return;
    }
}

template <typename T>
bool copyPropertyVector(
    std::vector<MeshLib::Element*> const& global_mesh_elements,
    MeshLib::Properties& partitioned_properties,
    MeshLib::Properties const& properties,
    std::vector<Partition> const& partitions,
    MeshLib::PropertyVector<T> const* const pv,
    std::map<MeshLib::MeshItemType, std::size_t> const& total_number_of_tuples)
{
    if (pv == nullptr)
    {
        return false;
    }
    auto const item_type = pv->getMeshItemType();

    std::size_t partitioned_pv_size = total_number_of_tuples.at(item_type) *
                                      pv->getNumberOfGlobalComponents();

    std::vector<std::size_t> element_ip_data_offsets;
    if (item_type == MeshLib::MeshItemType::IntegrationPoint)
    {
        // Special field data such as OGS_VERSION, IntegrationPointMetaData,
        // etc., which are not "real" integration points, are copied "as is"
        // (i.e. fully) for every partition.
        if (pv->getPropertyName().find("_ip") == std::string::npos)
        {
            partitioned_pv_size = pv->size() * partitions.size();
        }

        element_ip_data_offsets =
            MeshToolsLib::getIntegrationPointDataOffsetsOfMeshElements(
                global_mesh_elements, *pv, properties);
    }

    auto partitioned_pv = partitioned_properties.createNewPropertyVector<T>(
        pv->getPropertyName(), pv->getMeshItemType(),
        pv->getNumberOfGlobalComponents());
    partitioned_pv->resize(partitioned_pv_size);

    auto copy_property_vector_values =
        [&](Partition const& p, std::size_t offset)
    {
        if (item_type == MeshLib::MeshItemType::IntegrationPoint)
        {
            return copyFieldPropertyDataToPartitions(properties, p, offset,
                                                     element_ip_data_offsets,
                                                     *pv, *partitioned_pv);
        }

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
            toString(item_type));
    };

    std::size_t position_offset(0);
    for (auto p : partitions)
    {
        position_offset += copy_property_vector_values(p, position_offset);
    }
    return true;
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
    std::unique_ptr<MeshLib::Mesh> const& mesh,
    std::vector<Partition>& partitions)
{
    using namespace MeshLib;

    MeshLib::Properties const& properties = mesh->getProperties();

    // Count the number of integration point data of all partitions:
    setIntegrationPointNumberOfPartition(properties, partitions);

    Properties partitioned_properties;
    auto count_tuples = [&](MeshItemType const mesh_item_type)
    {
        return std::accumulate(
            begin(partitions), end(partitions), 0,
            [&](std::size_t const sum, Partition const& p)
            { return sum + p.numberOfMeshItems(mesh_item_type); });
    };

    std::map<MeshItemType, std::size_t> const total_number_of_tuples = {
        {MeshItemType::Cell, count_tuples(MeshItemType::Cell)},
        {MeshItemType::Node, count_tuples(MeshItemType::Node)},
        {MeshItemType::IntegrationPoint,
         count_tuples(MeshItemType::IntegrationPoint)}};

    DBUG(
        "total number of tuples after partitioning defined for cells is {:d} "
        "and for nodes {:d} and for integration points {:d}.",
        total_number_of_tuples.at(MeshItemType::Cell),
        total_number_of_tuples.at(MeshItemType::Node),
        total_number_of_tuples.at(MeshItemType::IntegrationPoint));

    //  1 create new PV
    //  2 resize the PV with total_number_of_tuples
    //  3 copy the values according to the partition info
    applyToPropertyVectors(
        properties,
        [&](auto type, auto const property)
        {
            return copyPropertyVector<decltype(type)>(
                mesh->getElements(), partitioned_properties, properties,
                partitions,
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

void checkFieldPropertyVectorSize(
    std::vector<MeshLib::Element*> const& global_mesh_elements,
    MeshLib::Properties const& properties)
{
    for (auto const& [name, property] : properties)
    {
        auto const item_type = property->getMeshItemType();

        if (item_type != MeshLib::MeshItemType::IntegrationPoint)
        {
            continue;
        }

        // For special field data such as OGS_VERSION, IntegrationPointMetaData,
        // etc., which are not "real" integration points:
        if (property->getPropertyName().find("_ip") == std::string::npos)
        {
            continue;
        }

        std::size_t number_of_total_integration_points = 0;
        auto const ip_meta_data = MeshLib::getIntegrationPointMetaData(
            properties, property->getPropertyName());
        for (auto const element : global_mesh_elements)
        {
            int const number_of_integration_points =
                MeshToolsLib::getNumberOfElementIntegrationPoints(ip_meta_data,
                                                                  *element);
            number_of_total_integration_points += number_of_integration_points;
        }

        const auto pv =
            dynamic_cast<MeshLib::PropertyVector<double> const*>(property);
        std::size_t const component_number = pv->getNumberOfGlobalComponents();
        if (pv->size() != number_of_total_integration_points * component_number)
        {
            OGS_FATAL(
                "The property vector's size {:d} for integration point data "
                "{:s} does not match its actual size {:d}. The field data in "
                "the vtu file are wrong.",
                pv->size(), name,
                number_of_total_integration_points * component_number);
        }
    }
}

std::vector<std::vector<std::size_t>> computePartitionIDPerElement(
    std::vector<std::size_t> const& node_partition_map,
    std::vector<MeshLib::Element*> const& elements,
    std::vector<std::size_t> const& bulk_node_ids)
{
    auto node_partition_ids = ranges::views::transform(
        [&](MeshLib::Element const* const element)
        {
            auto node_lookup = ranges::views::transform(
                [&](std::size_t const i)
                { return node_partition_map[bulk_node_ids[i]]; });

            return element->nodes() | MeshLib::views::ids | node_lookup |
                   ranges::to<std::vector>;
        });

    return elements | node_partition_ids | ranges::to<std::vector>;
}

void distributeNodesToPartitions(
    std::vector<Partition>& partitions,
    std::vector<std::size_t> const& nodes_partition_ids,
    std::vector<MeshLib::Node*> const& nodes,
    std::vector<std::size_t> const& bulk_node_ids)
{
    for (auto const* const node : nodes)
    {
        partitions[nodes_partition_ids[bulk_node_ids[node->getID()]]]
            .nodes.push_back(node);
    }
}

void reorderNodesIntoBaseAndHigherOrderNodes(Partition& partition,
                                             MeshLib::Mesh const& mesh)
{
    std::vector<MeshLib::Node const*> higher_order_nodes;
    // after splitIntoBaseAndHigherOrderNodes() partition.nodes contains only
    // base nodes
    std::tie(partition.nodes, higher_order_nodes) =
        splitIntoBaseAndHigherOrderNodes(partition.nodes, mesh);
    partition.number_of_regular_base_nodes = partition.nodes.size();
    std::copy(begin(higher_order_nodes), end(higher_order_nodes),
              std::back_inserter(partition.nodes));
    partition.number_of_regular_nodes = partition.nodes.size();
}

void reorderNodesIntoBaseAndHigherOrderNodesPerPartition(
    std::vector<Partition>& partitions, MeshLib::Mesh const& mesh)
{
    for (auto& partition : partitions)
    {
        reorderNodesIntoBaseAndHigherOrderNodes(partition, mesh);
    }
}

void setNumberOfNodesInPartitions(std::vector<Partition>& partitions,
                                  MeshLib::Mesh const& mesh)
{
    auto const number_of_mesh_base_nodes = mesh.computeNumberOfBaseNodes();
    auto const number_of_mesh_all_nodes = mesh.getNumberOfNodes();
    for (auto& partition : partitions)
    {
        partition.number_of_regular_nodes = partition.nodes.size();
        partition.number_of_mesh_base_nodes = number_of_mesh_base_nodes;
        partition.number_of_mesh_all_nodes = number_of_mesh_all_nodes;
    }
}

void distributeElementsIntoPartitions(
    std::vector<Partition>& partitions,
    MeshLib::Mesh const& mesh,
    std::vector<std::vector<std::size_t>> const& partition_ids_per_element)
{
    for (auto const& element : mesh.getElements())
    {
        auto const element_id = element->getID();
        auto node_partition_ids = partition_ids_per_element[element_id];
        // make partition ids unique
        std::sort(node_partition_ids.begin(), node_partition_ids.end());
        auto last =
            std::unique(node_partition_ids.begin(), node_partition_ids.end());
        node_partition_ids.erase(last, node_partition_ids.end());

        // all element nodes belong to the same partition => regular element
        if (node_partition_ids.size() == 1)
        {
            partitions[node_partition_ids[0]].regular_elements.push_back(
                element);
        }
        else
        {
            for (auto const partition_id : node_partition_ids)
            {
                partitions[partition_id].ghost_elements.push_back(element);
            }
        }
    }
}

// determine and append ghost nodes to partition.nodes in the following order
// [base nodes, higher order nodes, base ghost nodes, higher order ghost
// nodes]
void determineAndAppendGhostNodesToPartitions(
    std::vector<Partition>& partitions, MeshLib::Mesh const& mesh,
    std::vector<std::size_t> const& nodes_partition_ids,
    std::vector<std::size_t> const& node_id_mapping)
{
    for (std::size_t part_id = 0; part_id < partitions.size(); part_id++)
    {
        auto& partition = partitions[part_id];
        std::vector<MeshLib::Node*> base_ghost_nodes;
        std::vector<MeshLib::Node*> higher_order_ghost_nodes;
        std::tie(base_ghost_nodes, higher_order_ghost_nodes) =
            findGhostNodesInPartition(
                part_id, mesh.getNodes(), partition.ghost_elements,
                nodes_partition_ids, mesh, node_id_mapping);

        std::copy(begin(base_ghost_nodes), end(base_ghost_nodes),
                  std::back_inserter(partition.nodes));

        partition.number_of_base_nodes =
            partition.number_of_regular_base_nodes + base_ghost_nodes.size();

        std::copy(begin(higher_order_ghost_nodes),
                  end(higher_order_ghost_nodes),
                  std::back_inserter(partition.nodes));
    }
}

void partitionMesh(std::vector<Partition>& partitions,
                   MeshLib::Mesh const& mesh,
                   std::vector<std::size_t> const& nodes_partition_ids,
                   std::vector<std::size_t> const& bulk_node_ids)
{
    BaseLib::RunTime run_timer;
    run_timer.start();
    auto const partition_ids_per_element = computePartitionIDPerElement(
        nodes_partition_ids, mesh.getElements(), bulk_node_ids);
    INFO("partitionMesh(): Partition IDs per element computed in {:g} s",
         run_timer.elapsed());

    run_timer.start();
    distributeNodesToPartitions(partitions, nodes_partition_ids,
                                mesh.getNodes(), bulk_node_ids);
    INFO("partitionMesh(): distribute nodes to partitions took {:g} s",
         run_timer.elapsed());

    run_timer.start();
    reorderNodesIntoBaseAndHigherOrderNodesPerPartition(partitions, mesh);
    INFO(
        "partitionMesh(): sorting [base nodes | higher order nodes] took {:g} "
        "s",
        run_timer.elapsed());

    run_timer.start();
    setNumberOfNodesInPartitions(partitions, mesh);
    INFO(
        "partitionMesh(): setting number of nodes and of all mesh base nodes "
        "took {:g} s",
        run_timer.elapsed());

    run_timer.start();
    distributeElementsIntoPartitions(partitions, mesh,
                                     partition_ids_per_element);
    INFO("partitionMesh(): distribute elements into partitions took {:g} s",
         run_timer.elapsed());

    run_timer.start();
    determineAndAppendGhostNodesToPartitions(
        partitions, mesh, nodes_partition_ids, bulk_node_ids);
    INFO("partitionMesh(): determine / append ghost nodes took {:g} s",
         run_timer.elapsed());

    run_timer.start();
    markDuplicateGhostCells(mesh, partitions);
    INFO("partitionMesh(): markDuplicateGhostCells took {:g} s",
         run_timer.elapsed());
}

void NodeWiseMeshPartitioner::partitionByMETIS()
{
    std::vector<std::size_t> bulk_node_ids(_mesh->getNumberOfNodes());
    std::iota(bulk_node_ids.begin(), bulk_node_ids.end(), 0);

    partitionMesh(_partitions, *_mesh, _nodes_partition_ids, bulk_node_ids);

    renumberNodeIndices();

    // In case the field data in the vtu file are manually added, e.g. by using
    // some tools, the size of the field property vector has to be checked.
    checkFieldPropertyVectorSize(_mesh->getElements(), _mesh->getProperties());

    _partitioned_properties = partitionProperties(_mesh, _partitions);

    renumberBulkIdsProperty(_partitions, _partitioned_properties);
}

void NodeWiseMeshPartitioner::renumberBulkIdsProperty(
    std::vector<Partition> const& partitions,
    MeshLib::Properties& partitioned_properties)
{
    auto const bulk_node_ids_string =
        MeshLib::getBulkIDString(MeshLib::MeshItemType::Node);
    if (partitioned_properties.hasPropertyVector(bulk_node_ids_string))
    {
        renumberBulkNodeIdsProperty(
            partitioned_properties.getPropertyVector<std::size_t>(
                bulk_node_ids_string, MeshLib::MeshItemType::Node, 1),
            partitions);
    }
    auto const bulk_element_ids_string =
        MeshLib::getBulkIDString(MeshLib::MeshItemType::Cell);
    if (partitioned_properties.hasPropertyVector<std::size_t>(
            static_cast<std::string>(bulk_element_ids_string),
            MeshLib::MeshItemType::Cell))
    {
        renumberBulkElementIdsProperty(
            partitioned_properties.getPropertyVector<std::size_t>(
                bulk_element_ids_string, MeshLib::MeshItemType::Cell, 1),
            partitions);
    }
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
                std::size_t const offset)
        {
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
                std::size_t const offset)
        {
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
    MeshLib::Mesh const& mesh) const
{
    auto const bulk_node_ids_string =
        MeshLib::getBulkIDString(MeshLib::MeshItemType::Node);
    auto const& bulk_node_ids =
        mesh.getProperties().getPropertyVector<std::size_t>(
            bulk_node_ids_string, MeshLib::MeshItemType::Node, 1);

    std::vector<Partition> partitions(_partitions.size());

    partitionMesh(partitions, mesh, _nodes_partition_ids, *bulk_node_ids);

    return partitions;
}

void NodeWiseMeshPartitioner::renumberNodeIndices()
{
    std::size_t node_global_id_offset = 0;
    // Renumber the global indices.
    for (auto& partition : _partitions)
    {
        for (std::size_t i = 0; i < partition.number_of_regular_nodes; i++)
        {
            _nodes_global_ids[partition.nodes[i]->getID()] =
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
        partitioned_properties,
        [&](auto type, auto const& property)
        {
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
                          offsets.node * sizeof(MeshLib::IO::NodeData)),
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
    std::vector<MeshLib::Node const*> const& nodes)
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
                std::ofstream& output_stream)
        {
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
            output_stream.write(reinterpret_cast<const char*>(ele_info.data()),
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
    writeProperties(file_name_base, _partitioned_properties, _partitions,
                    MeshLib::MeshItemType::IntegrationPoint);

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
