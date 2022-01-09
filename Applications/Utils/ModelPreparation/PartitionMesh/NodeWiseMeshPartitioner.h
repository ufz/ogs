/*!
  \file
  \date   2016.05

  \brief  Declare a class to perform node wise mesh partitioning

  \copyright
  Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
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
    std::size_t number_of_regular_base_nodes;
    std::size_t number_of_regular_nodes;
    std::size_t number_of_base_nodes;
    std::size_t number_of_mesh_base_nodes;
    std::size_t number_of_mesh_all_nodes;
    /// Non ghost elements
    std::vector<const MeshLib::Element*> regular_elements;
    std::vector<const MeshLib::Element*> ghost_elements;
    std::vector<bool> duplicate_ghost_cell;

    std::size_t numberOfMeshItems(MeshLib::MeshItemType const item_type) const;

    std::ostream& writeNodes(
        std::ostream& os,
        std::vector<std::size_t> const& global_node_ids) const;

    std::ostream& writeConfig(std::ostream& os) const;
};

/// Creates partitioned mesh properties for nodes and cells.
MeshLib::Properties partitionProperties(
    MeshLib::Properties const& properties,
    std::vector<Partition> const& partitions);

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
        : _partitions(num_partitions),
          _partitioned_properties(),
          _mesh(std::move(mesh)),
          _nodes_global_ids(_mesh->getNumberOfNodes()),
          _nodes_partition_ids(_mesh->getNumberOfNodes())
    {
    }

    /// Partition by node.
    void partitionByMETIS();

    std::vector<Partition> partitionOtherMesh(MeshLib::Mesh const& mesh) const;

    /// Renumber the bulk_node_ids property for each partition to match the
    /// partitioned bulk mesh nodes.
    void renumberBulkNodeIdsProperty(
        MeshLib::PropertyVector<std::size_t>* const bulk_node_ids,
        std::vector<Partition> const& local_partitions) const;

    /// Renumber the bulk_element_ids property for each partition to match the
    /// partitioned bulk mesh elements.
    void renumberBulkElementIdsProperty(
        MeshLib::PropertyVector<std::size_t>* const bulk_element_ids_pv,
        std::vector<Partition> const& local_partitions) const;

    /// Write the partitions into binary files
    /// \param file_name_base The prefix of the file name.
    void write(const std::string& file_name_base);

    void writeOtherMesh(
        std::string const& output_filename_base,
        std::vector<Partition> const& partitions,
        MeshLib::Properties const& partitioned_properties) const;

    void resetPartitionIdsForNodes(
        std::vector<std::size_t>&& node_partition_ids)
    {
        _nodes_partition_ids = std::move(node_partition_ids);
    }

    MeshLib::Mesh const& mesh() const { return *_mesh; }

private:
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
    void renumberNodeIndices();

    void processPartition(std::size_t const part_id);
};

}  // namespace ApplicationUtils
