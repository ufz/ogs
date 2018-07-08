/*!
  \file NodeWiseMeshPartitioner.cpp
  \date   2016.05

  \brief  Define the members of class NodeWiseMeshPartitioner

  \copyright
  Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/

#include "NodeWiseMeshPartitioner.h"

#include <limits>
#include <numeric>
#include <unordered_map>

#include <logog/include/logog.hpp>

#include "BaseLib/Error.h"

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

std::ostream& Partition::writeNodesBinary(
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

std::ostream& Partition::writeConfigBinary(std::ostream& os) const
{
    long const data[] = {
        static_cast<long>(nodes.size()),
        static_cast<long>(number_of_base_nodes),
        static_cast<long>(regular_elements.size()),
        static_cast<long>(ghost_elements.size()),
        static_cast<long>(number_of_non_ghost_base_nodes),
        static_cast<long>(number_of_non_ghost_nodes),
        static_cast<long>(number_of_mesh_base_nodes),
        static_cast<long>(number_of_mesh_all_nodes),
        static_cast<long>(
            getNumberOfIntegerVariablesOfElements(regular_elements)),
        static_cast<long>(
            getNumberOfIntegerVariablesOfElements(ghost_elements)),
    };

    return os.write(reinterpret_cast<const char*>(data),
                    sizeof(data));
}

void NodeWiseMeshPartitioner::findNonGhostNodesInPartition(
    std::size_t const part_id,
    const bool is_mixed_high_order_linear_elems,
    std::vector<MeshLib::Node*>& extra_nodes)
{
    std::vector<MeshLib::Node*> const& nodes = _mesh->getNodes();
    auto& partition = _partitions[part_id];
    // -- Extra nodes for high order elements
    for (std::size_t i = 0; i < _mesh->getNumberOfNodes(); i++)
    {
        if (_nodes_partition_ids[i] == part_id)
        {
            splitOffHigherOrderNode(nodes, is_mixed_high_order_linear_elems, i,
                                    partition.nodes, extra_nodes);
        }
    }
    partition.number_of_non_ghost_base_nodes = partition.nodes.size();
    partition.number_of_non_ghost_nodes =
        partition.number_of_non_ghost_base_nodes + extra_nodes.size();
}

void NodeWiseMeshPartitioner::findElementsInPartition(std::size_t const part_id)
{
    auto& partition = _partitions[part_id];
    std::vector<MeshLib::Element*> const& elements = _mesh->getElements();
    std::vector<bool> _is_regular_element(elements.size(), false);

    for (std::size_t elem_id = 0; elem_id < elements.size(); elem_id++)
    {
        const auto* elem = elements[elem_id];
        if (_is_regular_element[elem_id])
        {
            continue;
        }

        std::size_t non_ghost_node_number = 0;
        for (unsigned i = 0; i < elem->getNumberOfNodes(); i++)
        {
            if (_nodes_partition_ids[elem->getNodeIndex(i)] == part_id)
            {
                non_ghost_node_number++;
            }
        }

        if (non_ghost_node_number == 0)
        {
            continue;
        }

        if (non_ghost_node_number == elem->getNumberOfNodes())
        {
            partition.regular_elements.push_back(elem);
            _is_regular_element[elem_id] = true;
        }
        else
        {
            partition.ghost_elements.push_back(elem);
        }
    }
}

void NodeWiseMeshPartitioner::findGhostNodesInPartition(
    std::size_t const part_id,
    const bool is_mixed_high_order_linear_elems,
    std::vector<MeshLib::Node*>& extra_nodes)
{
    auto& partition = _partitions[part_id];
    std::vector<MeshLib::Node*> const& nodes = _mesh->getNodes();
    std::vector<bool> nodes_reserved(_mesh->getNumberOfNodes(), false);
    for (const auto* ghost_elem : partition.ghost_elements)
    {
        for (unsigned i = 0; i < ghost_elem->getNumberOfNodes(); i++)
        {
            const unsigned node_id = ghost_elem->getNodeIndex(i);
            if (nodes_reserved[node_id])
            {
                continue;
            }

            if (_nodes_partition_ids[node_id] != part_id)
            {
                splitOffHigherOrderNode(nodes, is_mixed_high_order_linear_elems,
                                        node_id, partition.nodes, extra_nodes);
                nodes_reserved[node_id] = true;
            }
        }
    }
}

void NodeWiseMeshPartitioner::splitOffHigherOrderNode(
    std::vector<MeshLib::Node*> const& nodes,
    bool const is_mixed_high_order_linear_elems,
    unsigned const node_id,
    std::vector<MeshLib::Node*>& base_nodes,
    std::vector<MeshLib::Node*>& extra_nodes)
{
    auto const n_base_nodes = _mesh->getNumberOfBaseNodes();
    if (!is_mixed_high_order_linear_elems || node_id > n_base_nodes)
    {
        base_nodes.push_back(nodes[node_id]);
    }
    else
    {
        extra_nodes.push_back(nodes[node_id]);
    }
}

void NodeWiseMeshPartitioner::processPartition(
    std::size_t const part_id, const bool is_mixed_high_order_linear_elems)
{
    std::vector<MeshLib::Node*> extra_nodes;
    findNonGhostNodesInPartition(part_id, is_mixed_high_order_linear_elems,
                                 extra_nodes);

    findElementsInPartition(part_id);
    findGhostNodesInPartition(part_id, is_mixed_high_order_linear_elems,
                              extra_nodes);
    auto& partition = _partitions[part_id];
    partition.number_of_base_nodes = partition.nodes.size();

    if (is_mixed_high_order_linear_elems)
    {
        partition.nodes.insert(partition.nodes.end(), extra_nodes.begin(),
                               extra_nodes.end());
    }

    // Set the node numbers of base and all mesh nodes.
    partition.number_of_mesh_base_nodes = _mesh->getNumberOfBaseNodes();
    partition.number_of_mesh_all_nodes = _mesh->getNumberOfNodes();
}

template <typename T>
bool copyNodePropertyVector(MeshLib::Properties const& original_properties,
                            MeshLib::Properties& partitioned_properties,
                            std::vector<Partition> const& partitions,
                            std::string const& name,
                            std::size_t const total_number_of_tuples)
{
    if (!original_properties.existsPropertyVector<T>(name))
        return false;

    auto const& pv = original_properties.getPropertyVector<T>(name);
    auto partitioned_pv = partitioned_properties.createNewPropertyVector<T>(
        name, pv->getMeshItemType(), pv->getNumberOfComponents());
    partitioned_pv->resize(total_number_of_tuples *
                           pv->getNumberOfComponents());
    std::size_t position_offset(0);
    for (auto p : partitions)
    {
        for (std::size_t i = 0; i < p.nodes.size(); ++i)
        {
            const auto global_id = p.nodes[i]->getID();
            (*partitioned_pv)[position_offset + i] = (*pv)[global_id];
        }
        position_offset += p.nodes.size();
    }
    return true;
}

template <typename T>
bool copyCellPropertyVector(MeshLib::Properties const& original_properties,
                            MeshLib::Properties& partitioned_properties,
                            std::vector<Partition> const& partitions,
                            std::string const& name,
                            std::size_t const total_number_of_tuples)
{
    if (!original_properties.existsPropertyVector<T>(name))
        return false;

    auto const& pv = original_properties.getPropertyVector<T>(name);
    auto partitioned_pv = partitioned_properties.createNewPropertyVector<T>(
        name, pv->getMeshItemType(), pv->getNumberOfComponents());
    partitioned_pv->resize(total_number_of_tuples *
                           pv->getNumberOfComponents());
    std::size_t position_offset(0);
    for (auto const& p : partitions)
    {
        std::size_t const n_regular(p.regular_elements.size());
        for (std::size_t i = 0; i < n_regular; ++i)
        {
            const auto id = p.regular_elements[i]->getID();
            (*partitioned_pv)[position_offset + i] = (*pv)[id];
        }
        position_offset += n_regular;
        std::size_t const n_ghost(p.ghost_elements.size());
        for (std::size_t i = 0; i < n_ghost; ++i)
        {
            const auto id = p.ghost_elements[i]->getID();
            (*partitioned_pv)[position_offset + i] = (*pv)[id];
        }
        position_offset += n_ghost;
    }
    return true;
}

void NodeWiseMeshPartitioner::processNodeProperties()
{
    std::size_t const total_number_of_tuples =
        std::accumulate(std::begin(_partitions), std::end(_partitions), 0,
                        [](std::size_t const sum, Partition const& p) {
                            return sum + p.nodes.size();
                        });

    DBUG("total number of node-based tuples after partitioning: %d ",
         total_number_of_tuples);
    // 1 create new PV
    // 2 resize the PV with total_number_of_tuples
    // 3 copy the values according to the partition info
    auto const& original_properties(_mesh->getProperties());
    auto const property_names =
        original_properties.getPropertyVectorNames(MeshLib::MeshItemType::Node);
    for (auto const& name : property_names)
    {
        bool success =
            copyNodePropertyVector<double>(original_properties,
                                           _partitioned_properties, _partitions,
                                           name, total_number_of_tuples) ||
            copyNodePropertyVector<float>(original_properties,
                                          _partitioned_properties, _partitions,
                                          name, total_number_of_tuples) ||
            copyNodePropertyVector<int>(original_properties,
                                        _partitioned_properties, _partitions,
                                        name, total_number_of_tuples) ||
            copyNodePropertyVector<long>(original_properties,
                                         _partitioned_properties, _partitions,
                                         name, total_number_of_tuples) ||
            copyNodePropertyVector<unsigned>(
                original_properties, _partitioned_properties, _partitions, name,
                total_number_of_tuples) ||
            copyNodePropertyVector<unsigned long>(
                original_properties, _partitioned_properties, _partitions, name,
                total_number_of_tuples) ||
            copyNodePropertyVector<std::size_t>(
                original_properties, _partitioned_properties, _partitions, name,
                total_number_of_tuples);
        if (!success)
            WARN(
                "processNodeProperties: Could not create partitioned "
                "PropertyVector '%s'.",
                name.c_str());
    }
}

void NodeWiseMeshPartitioner::processCellProperties()
{
    std::size_t const total_number_of_tuples = std::accumulate(
        std::begin(_partitions), std::end(_partitions), 0,
        [](std::size_t const sum, Partition const& p) {
            return sum + p.regular_elements.size() + p.ghost_elements.size();
        });

    DBUG("total number of cell-based tuples after partitioning: %d ",
         total_number_of_tuples);
    // 1 create new PV
    // 2 resize the PV with total_number_of_tuples
    // 3 copy the values according to the partition info
    auto const& original_properties(_mesh->getProperties());
    auto const property_names =
        original_properties.getPropertyVectorNames(MeshLib::MeshItemType::Cell);
    for (auto const& name : property_names)
    {
        bool success =
            copyCellPropertyVector<double>(original_properties,
                                           _partitioned_properties, _partitions,
                                           name, total_number_of_tuples) ||
            copyCellPropertyVector<float>(original_properties,
                                          _partitioned_properties, _partitions,
                                          name, total_number_of_tuples) ||
            copyCellPropertyVector<int>(original_properties,
                                        _partitioned_properties, _partitions,
                                        name, total_number_of_tuples) ||
            copyCellPropertyVector<long>(original_properties,
                                         _partitioned_properties, _partitions,
                                         name, total_number_of_tuples) ||
            copyCellPropertyVector<unsigned>(
                original_properties, _partitioned_properties, _partitions, name,
                total_number_of_tuples) ||
            copyCellPropertyVector<unsigned long>(
                original_properties, _partitioned_properties, _partitions, name,
                total_number_of_tuples) ||
            copyCellPropertyVector<std::size_t>(
                original_properties, _partitioned_properties, _partitions, name,
                total_number_of_tuples);
        if (!success)
            WARN(
                "processCellProperties: Could not create partitioned "
                "PropertyVector '%s'.",
                name.c_str());
    }
}

void NodeWiseMeshPartitioner::partitionByMETIS(
    const bool is_mixed_high_order_linear_elems)
{
    for (std::size_t part_id = 0; part_id < _partitions.size(); part_id++)
    {
        INFO("Processing partition: %d", part_id);
        processPartition(part_id, is_mixed_high_order_linear_elems);
    }

    renumberNodeIndices(is_mixed_high_order_linear_elems);

    processNodeProperties();
    processCellProperties();
}

void NodeWiseMeshPartitioner::renumberNodeIndices(
    const bool is_mixed_high_order_linear_elems)
{
    std::size_t node_global_id_offset = 0;
    // Renumber the global indices.
    // -- Base nodes
    for (auto& partition : _partitions)
    {
        for (std::size_t i = 0; i < partition.number_of_non_ghost_base_nodes;
             i++)
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
                                   partition.number_of_non_ghost_nodes -
                                   partition.number_of_non_ghost_base_nodes;
        for (std::size_t i = partition.number_of_base_nodes; i < end_id; i++)
        {
            _nodes_global_ids[partition.nodes[i]->getID()] =
                node_global_id_offset;
            node_global_id_offset++;
        }
    }
}

template <typename T>
void writePropertyVectorValuesBinary(std::ostream& os,
                                     MeshLib::PropertyVector<T> const& pv)
{
    os.write(reinterpret_cast<const char*>(pv.data()),
            pv.size() * sizeof(T));
}

template <typename T>
bool writePropertyVectorBinary(
    MeshLib::Properties const& partitioned_properties, std::string const& name,
    std::ostream& out_val, std::ostream& out_meta)
{
    if (!partitioned_properties.existsPropertyVector<T>(name))
        return false;

    MeshLib::IO::PropertyVectorMetaData pvmd;
    pvmd.property_name = name;
    auto* pv = partitioned_properties.getPropertyVector<T>(name);
    pvmd.fillPropertyVectorMetaDataTypeInfo<T>();
    pvmd.number_of_components = pv->getNumberOfComponents();
    pvmd.number_of_tuples = pv->getNumberOfTuples();
    writePropertyVectorValuesBinary(out_val, *pv);
    MeshLib::IO::writePropertyVectorMetaDataBinary(out_meta, pvmd);
    return true;
}

void writeNodePropertiesBinary(
    const std::string& file_name_base,
    MeshLib::Properties const& partitioned_properties,
    std::vector<Partition> const& partitions)
{
    auto const& property_names = partitioned_properties.getPropertyVectorNames(
        MeshLib::MeshItemType::Node);
    if (property_names.empty())
    {
        return;
    }

    std::size_t const number_of_properties(property_names.size());

    std::ofstream out = BaseLib::createBinaryFile(
        file_name_base + "_partitioned_node_properties_cfg" +
        std::to_string(partitions.size()) + ".bin");

    std::ofstream out_val = BaseLib::createBinaryFile(
        file_name_base + "_partitioned_node_properties_val" +
        std::to_string(partitions.size()) + ".bin");

    out.write(reinterpret_cast<const char*>(&number_of_properties),
              sizeof(number_of_properties));
    for (auto const& name : property_names)
    {
        bool success = writePropertyVectorBinary<double>(partitioned_properties,
                                                         name, out_val, out) ||
                       writePropertyVectorBinary<float>(partitioned_properties,
                                                        name, out_val, out) ||
                       writePropertyVectorBinary<int>(partitioned_properties,
                                                      name, out_val, out) ||
                       writePropertyVectorBinary<long>(partitioned_properties,
                                                       name, out_val, out) ||
                       writePropertyVectorBinary<unsigned>(
                           partitioned_properties, name, out_val, out) ||
                       writePropertyVectorBinary<unsigned long>(
                           partitioned_properties, name, out_val, out) ||
                       writePropertyVectorBinary<std::size_t>(
                           partitioned_properties, name, out_val, out);
        if (!success)
        {
            OGS_FATAL(
                "writeNodePropertiesBinary: Could not write PropertyVector "
                "'%s'.",
                name.c_str());
        }
    }

    unsigned long offset = 0;
    for (const auto& partition : partitions)
    {
        MeshLib::IO::PropertyVectorPartitionMetaData pvpmd{};
        pvpmd.offset = offset;
        pvpmd.number_of_tuples = partition.nodes.size();
        DBUG(
            "Write meta data for node-based PropertyVector: global offset %d, "
            "number of tuples %d",
            pvpmd.offset, pvpmd.number_of_tuples);
        MeshLib::IO::writePropertyVectorPartitionMetaData(out, pvpmd);
        offset += pvpmd.number_of_tuples;
    }
}

void writeCellPropertiesBinary(
    const std::string& file_name_base,
    MeshLib::Properties const& partitioned_properties,
    std::vector<Partition> const& partitions)
{
    auto const& property_names = partitioned_properties.getPropertyVectorNames(
        MeshLib::MeshItemType::Cell);
    if (property_names.empty())
    {
        return;
    }

    std::size_t const number_of_properties(property_names.size());

    std::ofstream out = BaseLib::createBinaryFile(
        file_name_base + "_partitioned_cell_properties_cfg" +
        std::to_string(partitions.size()) + ".bin");

    std::ofstream out_val = BaseLib::createBinaryFile(
        file_name_base + "_partitioned_cell_properties_val" +
        std::to_string(partitions.size()) + ".bin");

    out.write(reinterpret_cast<const char*>(&number_of_properties),
              sizeof(number_of_properties));
    for (auto const& name : property_names)
    {
        bool success = writePropertyVectorBinary<double>(partitioned_properties,
                                                         name, out_val, out) ||
                       writePropertyVectorBinary<float>(partitioned_properties,
                                                        name, out_val, out) ||
                       writePropertyVectorBinary<int>(partitioned_properties,
                                                      name, out_val, out) ||
                       writePropertyVectorBinary<long>(partitioned_properties,
                                                       name, out_val, out) ||
                       writePropertyVectorBinary<unsigned>(
                           partitioned_properties, name, out_val, out) ||
                       writePropertyVectorBinary<unsigned long>(
                           partitioned_properties, name, out_val, out) ||
                       writePropertyVectorBinary<std::size_t>(
                           partitioned_properties, name, out_val, out);
        if (!success)
        {
            OGS_FATAL(
                "writeCellPropertiesBinary: Could not write PropertyVector "
                "'%s'.",
                name.c_str());
        }
    }

    unsigned long offset = 0;
    for (const auto& partition : partitions)
    {
        MeshLib::IO::PropertyVectorPartitionMetaData pvpmd{};
        pvpmd.offset = offset;
        pvpmd.number_of_tuples =
            partition.regular_elements.size() + partition.ghost_elements.size();
        DBUG(
            "Write meta data for cell-based PropertyVector: global offset %d, "
            "number of tuples %d",
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

    std::ostream& writeConfigBinary(std::ostream& os) const;
};

std::ostream& ConfigOffsets::writeConfigBinary(std::ostream& os) const
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

PartitionOffsets
computePartitionElementOffsets(Partition const& partition)
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
        // Offset the ending entry of the element integer variales of
        // the non-ghost elements of this partition in the vector of elem_info.
        static_cast<long>(oldConfig.element_rank_offset +
                          offsets.regular_elements * sizeof(long)),

        // Offset the ending entry of the element integer variales of
        // the ghost elements of this partition in the vector of elem_info.
        static_cast<long>(oldConfig.ghost_element_rank_offset +
                          offsets.ghost_elements * sizeof(long))};
}

/// Write the configuration data of the partition data in binary files.
/// \return a pair of vectors for:
///  1. The number of all non-ghost element integer variables for each
///     partition.
///  2. The number of all ghost element integer variables for each partition.
std::tuple<std::vector<long>, std::vector<long>> writeConfigDataBinary(
    const std::string& file_name_base,
    std::vector<Partition> const& partitions)
{
    std::ofstream of_bin_cfg =
        BaseLib::createBinaryFile(file_name_base + "_partitioned_msh_cfg" +
                                  std::to_string(partitions.size()) + ".bin");

    std::vector<long> num_elem_integers;
    num_elem_integers.reserve(partitions.size());
    std::vector<long> num_g_elem_integers;
    num_g_elem_integers.reserve(partitions.size());

    ConfigOffsets config_offsets = {0, 0, 0};  // 0 for first partition.
    for (const auto& partition : partitions)
    {
        partition.writeConfigBinary(of_bin_cfg);

        config_offsets.writeConfigBinary(of_bin_cfg);
        auto const& new_offsets = computePartitionElementOffsets(partition);
        config_offsets = incrementConfigOffsets(config_offsets, new_offsets);

        num_elem_integers.push_back(new_offsets.regular_elements);
        num_g_elem_integers.push_back(new_offsets.ghost_elements);
    }

    return std::make_tuple(num_elem_integers, num_g_elem_integers);
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
    unsigned mat_id = 0;  // TODO: Material ID to be set from the mesh data
    const long nn = elem.getNumberOfNodes();
    elem_info[counter++] = mat_id;
    elem_info[counter++] = static_cast<long>(elem.getCellType());
    elem_info[counter++] = nn;

    for (long i = 0; i < nn; i++)
    {
        elem_info[counter++] = local_node_ids.at(elem.getNodeIndex(i));
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
/// \param num_elem_integers    The numbers of all non-ghost element
///                             integer variables of each partitions.
/// \param num_g_elem_integers  The numbers of all ghost element
void writeElementsBinary(std::string const& file_name_base,
                         std::vector<Partition> const& partitions,
                         std::vector<long> const& num_elem_integers,
                         std::vector<long> const& num_g_elem_integers)
{
    const std::string npartitions_str = std::to_string(partitions.size());

    std::ofstream element_info_os = BaseLib::createBinaryFile(
        file_name_base + "_partitioned_msh_ele" + npartitions_str + ".bin");
    std::ofstream ghost_element_info_os = BaseLib::createBinaryFile(
        file_name_base + "_partitioned_msh_ele_g" + npartitions_str + ".bin");

    for (std::size_t i = 0; i < partitions.size(); i++)
    {
        const auto& partition = partitions[i];
        auto const local_node_ids = enumerateLocalNodeIds(partition.nodes);

        // A vector contians all element integer variables of
        // the non-ghost elements of this partition
        std::vector<long> ele_info(num_elem_integers[i]);

        // Non-ghost elements.
        long counter = partition.regular_elements.size();

        for (std::size_t j = 0; j < partition.regular_elements.size(); j++)
        {
            const auto* elem = partition.regular_elements[j];
            ele_info[j] = counter;
            getElementIntegerVariables(*elem, local_node_ids, ele_info,
                                       counter);
        }
        // Write vector data of non-ghost elements
        element_info_os.write(reinterpret_cast<const char*>(ele_info.data()),
                              ele_info.size() * sizeof(long));

        // Ghost elements
        ele_info.resize(num_g_elem_integers[i]);

        counter = partition.ghost_elements.size();

        for (std::size_t j = 0; j < partition.ghost_elements.size(); j++)
        {
            const auto* elem = partition.ghost_elements[j];
            ele_info[j] = counter;
            getElementIntegerVariables(*elem, local_node_ids, ele_info,
                                       counter);
        }
        // Write vector data of ghost elements
        ghost_element_info_os.write(
            reinterpret_cast<const char*>(ele_info.data()),
            ele_info.size() * sizeof(long));
    }
}

/// Write the nodes of all partitions into a binary file.
/// \param file_name_base The prefix of the file name.
/// \param partitions the list of partitions
/// \param global_node_ids global numbering of nodes
void writeNodesBinary(const std::string& file_name_base,
                      std::vector<Partition> const& partitions,
                      std::vector<std::size_t> const& global_node_ids)
{
    std::ofstream os =
        BaseLib::createBinaryFile(file_name_base + "_partitioned_msh_nod" +
                                  std::to_string(partitions.size()) + ".bin");

    for (const auto& partition : partitions)
    {
        partition.writeNodesBinary(os, global_node_ids);
    }
}

void NodeWiseMeshPartitioner::writeBinary(const std::string& file_name_base)
{
    writeNodePropertiesBinary(file_name_base, _partitioned_properties,
                              _partitions);
    writeCellPropertiesBinary(file_name_base, _partitioned_properties,
                              _partitions);

    const auto elem_integers =
        writeConfigDataBinary(file_name_base, _partitions);

    const std::vector<IntegerType>& num_elem_integers =
        std::get<0>(elem_integers);
    const std::vector<IntegerType>& num_g_elem_integers =
        std::get<1>(elem_integers);
    writeElementsBinary(file_name_base, _partitions, num_elem_integers,
                        num_g_elem_integers);

    writeNodesBinary(file_name_base, _partitions, _nodes_global_ids);
}

void NodeWiseMeshPartitioner::writeConfigDataASCII(
    const std::string& file_name_base)
{
    const std::string fname = file_name_base + "_partitioned_cfg" +
                              std::to_string(_npartitions) + ".msh";
    std::fstream os_subd_head(fname, std::ios::out | std::ios::trunc);
    const std::string mesh_info =
        "Subdomain mesh ("
        "Number of nodes; Number of base nodes;"
        " Number of regular elements; Number of ghost elements;"
        " Number of non-ghost base nodes; Number of non-ghost nodes"
        " Number of base nodes of the global mesh;"
        " Number of nodes of the global mesh;"
        " Number of integer variables to define non-ghost elements;"
        " Number of integer variables to define ghost elements.)";
    os_subd_head << mesh_info << "\n";
    os_subd_head << _npartitions << "\n";

    for (const auto& partition : _partitions)
    {
        os_subd_head << partition.nodes.size();
        os_subd_head << " " << partition.number_of_base_nodes;
        os_subd_head << " " << partition.regular_elements.size();
        os_subd_head << " " << partition.ghost_elements.size();
        os_subd_head << " " << partition.number_of_non_ghost_base_nodes;
        os_subd_head << " " << partition.number_of_non_ghost_nodes;
        os_subd_head << " " << _mesh->getNumberOfBaseNodes();
        os_subd_head << " " << _mesh->getNumberOfNodes();
        os_subd_head << " "
                     << getNumberOfIntegerVariablesOfElements(
                            partition.regular_elements);
        os_subd_head << " "
                     << getNumberOfIntegerVariablesOfElements(
                            partition.ghost_elements)
                     << " 0\n";
    }
}

void NodeWiseMeshPartitioner::writeElementsASCII(
    const std::string& file_name_base)
{
    const std::string fname = file_name_base + "_partitioned_elems_" +
                              std::to_string(_npartitions) + ".msh";
    std::fstream os_subd(fname, std::ios::out | std::ios::trunc);
    for (const auto& partition : _partitions)
    {
        // Set the local node indices of the current partition.
        IntegerType node_local_id_offset = 0;
        std::vector<IntegerType> nodes_local_ids(_mesh->getNumberOfNodes(), -1);
        for (const auto* node : partition.nodes)
        {
            nodes_local_ids[node->getID()] = node_local_id_offset;
            node_local_id_offset++;
        }

        for (const auto* elem : partition.regular_elements)
        {
            writeLocalElementNodeIndices(os_subd, *elem, nodes_local_ids);
        }
        for (const auto* elem : partition.ghost_elements)
        {
            writeLocalElementNodeIndices(os_subd, *elem, nodes_local_ids);
        }
        os_subd << std::endl;
    }
}

void NodeWiseMeshPartitioner::writeNodesASCII(const std::string& file_name_base)
{
    const std::string fname = file_name_base + "_partitioned_nodes_" +
                              std::to_string(_npartitions) + ".msh";
    std::fstream os_subd_node(fname, std::ios::out | std::ios::trunc);
    os_subd_node.precision(std::numeric_limits<double>::digits10);
    os_subd_node.setf(std::ios::scientific);

    for (const auto& partition : _partitions)
    {
        for (const auto* node : partition.nodes)
        {
            double const* coords = node->getCoords();
            os_subd_node << _nodes_global_ids[node->getID()] << " " << coords[0]
                         << " " << coords[1] << " " << coords[2] << "\n";
        }
        os_subd_node << std::endl;
    }
}

void NodeWiseMeshPartitioner::writeASCII(const std::string& file_name_base)
{
    writeConfigDataASCII(file_name_base);
    writeElementsASCII(file_name_base);
    writeNodesASCII(file_name_base);
}

void NodeWiseMeshPartitioner::writeLocalElementNodeIndices(
    std::ostream& os,
    const MeshLib::Element& elem,
    const std::vector<IntegerType>& local_node_ids)
{
    unsigned mat_id = 0;  // TODO: Material ID to be set from the mesh data
    os << mat_id << " " << static_cast<unsigned>(elem.getCellType()) << " "
       << elem.getNumberOfNodes() << " ";
    for (unsigned i = 0; i < elem.getNumberOfNodes(); i++)
    {
        os << " " << local_node_ids[elem.getNodeIndex(i)];
    }
    os << "\n";
}

}  // namespace ApplicationUtils
