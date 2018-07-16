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

void splitOffHigherOrderNode(std::vector<MeshLib::Node*> const& nodes,
                             bool const is_mixed_high_order_linear_elems,
                             unsigned const node_id,
                             unsigned const n_base_nodes,
                             std::vector<MeshLib::Node*>& base_nodes,
                             std::vector<MeshLib::Node*>& extra_nodes)
{
    if (!is_mixed_high_order_linear_elems || node_id > n_base_nodes)
    {
        base_nodes.push_back(nodes[node_id]);
    }
    else
    {
        extra_nodes.push_back(nodes[node_id]);
    }
}

/// 1 copy pointers to nodes belonging to the partition part_id into base nodes
/// vector, and
/// 2 collect non-linear element nodes belonging to the partition part_id in
/// extra nodes vector.
/// \return a pair of base node and extra nodes.
std::pair<std::vector<MeshLib::Node*>, std::vector<MeshLib::Node*>>
findNonGhostNodesInPartition(std::size_t const part_id,
                             const bool is_mixed_high_order_linear_elems,
                             std::size_t const n_base_nodes,
                             std::vector<MeshLib::Node*> const& nodes,
                             std::vector<std::size_t> const& partition_ids)
{
    auto const n_nodes = nodes.size();

    std::vector<MeshLib::Node*> base_nodes;
    std::vector<MeshLib::Node*> extra_nodes;
    for (std::size_t i = 0; i < n_nodes; i++)
    {
        if (partition_ids[i] == part_id)
        {
            splitOffHigherOrderNode(nodes, is_mixed_high_order_linear_elems, i,
                                    n_base_nodes, base_nodes, extra_nodes);
        }
    }

    return {base_nodes, extra_nodes};
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
                                        node_id, _mesh->getNumberOfBaseNodes(),
                                        partition.nodes, extra_nodes);
                nodes_reserved[node_id] = true;
            }
        }
    }
}

void NodeWiseMeshPartitioner::processPartition(
    std::size_t const part_id, const bool is_mixed_high_order_linear_elems)
{
    std::vector<MeshLib::Node*> extra_nodes;
    std::tie(_partitions[part_id].nodes, extra_nodes) =
        findNonGhostNodesInPartition(part_id, is_mixed_high_order_linear_elems,
                                     _mesh->getNumberOfBaseNodes(),
                                     _mesh->getNodes(), _nodes_partition_ids);

    _partitions[part_id].number_of_non_ghost_base_nodes =
        _partitions[part_id].nodes.size();
    _partitions[part_id].number_of_non_ghost_nodes =
        _partitions[part_id].number_of_non_ghost_base_nodes +
        extra_nodes.size();

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
    for (std::size_t i = 0; i < nnodes; ++i)
    {
        const auto global_id = nodes[i]->getID();
        partitioned_pv[offset + i] = pv[global_id];
    }
    return nnodes;
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
    for (std::size_t i = 0; i < n_regular; ++i)
    {
        const auto id = p.regular_elements[i]->getID();
        partitioned_pv[offset + i] = pv[id];
    }

    std::size_t const n_ghost(p.ghost_elements.size());
    for (std::size_t i = 0; i < n_ghost; ++i)
    {
        const auto id = p.ghost_elements[i]->getID();
        partitioned_pv[offset + n_regular + i] = pv[id];
    }
    return n_regular + n_ghost;
}

template <typename T>
bool copyPropertyVector(MeshLib::Properties const& original_properties,
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

    auto copy_property_vector_values = [&](Partition const& p,
                                           std::size_t offset) {
        if (pv->getMeshItemType() == MeshLib::MeshItemType::Node)
        {
            return copyNodePropertyVectorValues(p, offset, *pv,
                                                *partitioned_pv);
        }
        if (pv->getMeshItemType() == MeshLib::MeshItemType::Cell)
        {
            return copyCellPropertyVectorValues(p, offset, *pv,
                                                *partitioned_pv);
        }
        OGS_FATAL(
            "Copying of property vector values for mesh item type %s is "
            "not implemented.",
            pv->getMeshItemType());
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
void applyToPropertyVectors(std::vector<std::string> const& property_names,
                            Function f)
{
    for (auto const& name : property_names)
    {
        // Open question, why is the 'unsigned long' case not compiling giving
        // an error "expected '(' for function-style cast or type construction"
        // with clang-7, and "error C4576: a parenthesized type followed by an
        // initializer list is a non-standard explicit type conversion syntax"
        // with MSVC-15.
        bool success =
            f(double{}, name) || f(float{}, name) || f(int{}, name) ||
            f(long{}, name) || f(unsigned{}, name) ||
            f(static_cast<unsigned long>(0), name) || f(std::size_t{}, name);
        if (!success)
        {
            OGS_FATAL("Could not apply function to PropertyVector '%s'.",
                      name.c_str());
        }
    }
}

void NodeWiseMeshPartitioner::processProperties(
    MeshLib::MeshItemType const mesh_item_type)
{
    std::size_t const total_number_of_tuples =
        std::accumulate(std::begin(_partitions), std::end(_partitions), 0,
                        [&](std::size_t const sum, Partition const& p) {
                            return sum + p.numberOfMeshItems(mesh_item_type);
                        });

    DBUG(
        "total number of tuples define on mesh item type '%d' after "
        "partitioning: %d ",
        mesh_item_type, total_number_of_tuples);

    // 1 create new PV
    // 2 resize the PV with total_number_of_tuples
    // 3 copy the values according to the partition info
    auto const& original_properties(_mesh->getProperties());

    applyToPropertyVectors(
        original_properties.getPropertyVectorNames(mesh_item_type),
        [&](auto type, std::string const& name) {
            return copyPropertyVector<decltype(type)>(
                original_properties, _partitioned_properties, _partitions, name,
                total_number_of_tuples);
        });
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

    processProperties(MeshLib::MeshItemType::Node);
    processProperties(MeshLib::MeshItemType::Cell);
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

void writePropertiesBinary(const std::string& file_name_base,
                           MeshLib::Properties const& partitioned_properties,
                           std::vector<Partition> const& partitions,
                           MeshLib::MeshItemType const mesh_item_type)
{
    auto const& property_names =
        partitioned_properties.getPropertyVectorNames(mesh_item_type);
    if (property_names.empty())
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
        OGS_FATAL("Could not open file '%s' for output.",
                  file_name_cfg.c_str());
    }

    auto const file_name_val = file_name_base + "_partitioned_" +
                               file_name_infix + "_properties_val" +
                               std::to_string(partitions.size()) + ".bin";
    std::ofstream out_val(file_name_val, std::ios::binary);
    if (!out_val)
    {
        OGS_FATAL("Could not open file '%s' for output.",
                  file_name_val.c_str());
    }

    std::size_t const number_of_properties(property_names.size());
    BaseLib::writeValueBinary(out, number_of_properties);

    applyToPropertyVectors(property_names,
                         [&](auto type, std::string const& name) {
                             return writePropertyVectorBinary<decltype(type)>(
                                 partitioned_properties, name, out_val, out);
                         });

    unsigned long offset = 0;
    for (const auto& partition : partitions)
    {
        MeshLib::IO::PropertyVectorPartitionMetaData pvpmd{
            offset, static_cast<unsigned long>(
                        partition.numberOfMeshItems(mesh_item_type))};
        DBUG(
            "Write meta data for node-based PropertyVector: global offset %d, "
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
    auto const file_name_cfg = file_name_base + "_partitioned_msh_cfg" +
                               std::to_string(partitions.size()) + ".bin";
    std::ofstream of_bin_cfg(file_name_cfg, std::ios::binary);
    if (!of_bin_cfg)
    {
        OGS_FATAL("Could not open file '%s' for output.",
                  file_name_cfg.c_str());
    }

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

    auto const file_name_ele =
        file_name_base + "_partitioned_msh_ele" + npartitions_str + ".bin";
    std::ofstream element_info_os(file_name_ele, std::ios::binary);
    if (!element_info_os)
    {
        OGS_FATAL("Could not open file '%s' for output.",
                  file_name_ele.c_str());
    }

    auto const file_name_ele_g =
        file_name_base + "_partitioned_msh_ele_g" + npartitions_str + ".bin";
    std::ofstream ghost_element_info_os(file_name_ele_g, std::ios::binary);
    if (!ghost_element_info_os)
    {
        OGS_FATAL("Could not open file '%s' for output.",
                  file_name_ele_g.c_str());
    }

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
    auto const file_name = file_name_base + "_partitioned_msh_nod" +
                           std::to_string(partitions.size()) + ".bin";
    std::ofstream os(file_name, std::ios::binary);
    if (!os)
    {
        OGS_FATAL("Could not open file '%s' for output.", file_name.c_str());
    }

    for (const auto& partition : partitions)
    {
        partition.writeNodesBinary(os, global_node_ids);
    }
}

void NodeWiseMeshPartitioner::writeBinary(const std::string& file_name_base)
{
    writePropertiesBinary(file_name_base, _partitioned_properties, _partitions,
                          MeshLib::MeshItemType::Node);
    writePropertiesBinary(file_name_base, _partitioned_properties, _partitions,
                          MeshLib::MeshItemType::Cell);

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
