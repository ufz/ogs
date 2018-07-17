/*!
  \file NodeWiseMeshPartitioner.h
  \date   2016.05

  \brief  Declare a class to perform node wise mesh partitioning

  \copyright
  Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/

#pragma once

#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/MPI_IO/PropertyVectorMetaData.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

namespace ApplicationUtils
{
///  A subdomain mesh.
struct Partition
{
    std::vector<MeshLib::Node*> nodes;  ///< nodes.
    std::size_t number_of_non_ghost_base_nodes;
    std::size_t number_of_non_ghost_nodes;
    std::size_t number_of_base_nodes;
    std::size_t number_of_mesh_base_nodes;
    std::size_t number_of_mesh_all_nodes;
    /// Non ghost elements
    std::vector<const MeshLib::Element*> regular_elements;
    std::vector<const MeshLib::Element*> ghost_elements;

    std::size_t numberOfMeshItems(MeshLib::MeshItemType const item_type) const;

    std::ostream& writeNodesBinary(
        std::ostream& os,
        std::vector<std::size_t> const& global_node_ids) const;

    std::ostream& writeConfigBinary(std::ostream& os) const;
};

/// Mesh partitioner.
class NodeWiseMeshPartitioner
{
public:
    using IntegerType = long;

public:
    /*!
     * \param num_partitions Number of partitions,
     * \param mesh           Pointer to a mesh object.
     */
    NodeWiseMeshPartitioner(const IntegerType num_partitions,
                            std::unique_ptr<MeshLib::Mesh>&& mesh)
        : _npartitions(num_partitions),
          _partitions(num_partitions),
          _partitioned_properties(),
          _mesh(std::move(mesh)),
          _nodes_global_ids(_mesh->getNumberOfNodes()),
          _nodes_partition_ids(_mesh->getNumberOfNodes())
    {
    }

    /// Partition by node.
    /// \param is_mixed_high_order_linear_elems Flag to indicate whether the
    /// elements of a mesh can be used for both linear and high order
    /// interpolation
    void partitionByMETIS(const bool is_mixed_high_order_linear_elems);

    /// Write the partitions into ASCII files
    /// \param file_name_base The prefix of the file name.
    void writeASCII(const std::string& file_name_base);

    /// Write the partitions into binary files
    /// \param file_name_base The prefix of the file name.
    void writeBinary(const std::string& file_name_base);

    void resetPartitionIdsForNodes(
        std::vector<std::size_t>&& node_partition_ids)
    {
        _nodes_partition_ids = std::move(node_partition_ids);
    }

    MeshLib::Mesh const& mesh() const { return *_mesh; }

private:
    /// Number of partitions.
    IntegerType _npartitions;

    /// Data for all  partitions.
    std::vector<Partition> _partitions;

    /// Properties where values at ghost nodes and extra nodes are inserted.
    MeshLib::Properties _partitioned_properties;

    /// Pointer to a mesh object.
    std::unique_ptr<MeshLib::Mesh> _mesh;

    /// Global IDs of all nodes after partitioning.
    std::vector<std::size_t> _nodes_global_ids;

    /// Partition IDs of each nodes.
    std::vector<std::size_t> _nodes_partition_ids;

    // Renumber the global indices of nodes,
    /// \param is_mixed_high_order_linear_elems Flag to indicate whether the
    /// elements of a mesh can be used for both linear and high order
    /// interpolation
    void renumberNodeIndices(const bool is_mixed_high_order_linear_elems);

    /// Prerequisite: the ghost elements has to be found (using
    /// findElementsInPartition).
    /// Finds ghost nodes and non-linear element ghost nodes by walking over
    /// ghost elements.
    void findGhostNodesInPartition(std::size_t const part_id,
                                   const bool is_mixed_high_order_linear_elems,
                                   std::vector<MeshLib::Node*>& extra_nodes);

    void processPartition(std::size_t const part_id,
                          const bool is_mixed_high_order_linear_elems);

    void processProperties(MeshLib::MeshItemType const mesh_item_type);

    /// Write the configuration data of the partition data in ASCII files.
    /// \param file_name_base The prefix of the file name.
    void writeConfigDataASCII(const std::string& file_name_base);

    ///  Write the element integer variables of all partitions into
    ///  ASCII files.
    /// \param file_name_base The prefix of the file name.
    void writeElementsASCII(const std::string& file_name_base);

    ///  Write the nodes of all partitions into a ASCII file.
    /// \param file_name_base The prefix of the file name.
    void writeNodesASCII(const std::string& file_name_base);

    /*!
        \brief Write local indices of element nodes to a ASCII file
        \param os              Output stream
        \param elem            Element
        \param local_node_ids  Local node indices of a partition
    */
    void writeLocalElementNodeIndices(
        std::ostream& os,
        const MeshLib::Element& elem,
        const std::vector<IntegerType>& local_node_ids);
};

}  // namespace ApplicationUtils
