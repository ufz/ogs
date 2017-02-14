/*!
  \file NodeWiseMeshPartitioner.h
  \date   2016.05

  \brief  Declare a class to perform node wise mesh partitioning

  \copyright
  Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/

#pragma once

#include <memory>
#include <tuple>
#include <vector>
#include <string>
#include <fstream>

#include "MeshLib/Mesh.h"

namespace ApplicationUtils
{
///  A subdomain mesh.
struct Partition
{
    std::vector<MeshLib::Node*> nodes;  ///< nodes.
    std::size_t number_of_non_ghost_base_nodes;
    std::size_t number_of_non_ghost_nodes;
    std::size_t number_of_base_nodes;
    /// Non ghost elements
    std::vector<const MeshLib::Element*> regular_elements;
    std::vector<const MeshLib::Element*> ghost_elements;
};

/// Mesh partitioner.
class NodeWiseMeshPartitioner
{
public:
    typedef long IntegerType;

public:
    /*!
     * \param num_partitions Number of partitions,
     * \param mesh           Pointer to a mesh object.
     */
    NodeWiseMeshPartitioner(const IntegerType num_partitions,
                            std::unique_ptr<MeshLib::Mesh>&& mesh)
        : _npartitions(num_partitions),
          _partitions(num_partitions),
          _mesh(std::move(mesh)),
          _nodes_global_ids(_mesh->getNumberOfNodes()),
          _nodes_partition_ids(_mesh->getNumberOfNodes()),
          _elements_status(_mesh->getNumberOfElements(), false)
    {
    }

    /// Partition by node.
    /// \param is_mixed_high_order_linear_elems Flag to indicate whether the
    /// elements of a mesh can be used for both linear and high order
    /// interpolation
    void partitionByMETIS(const bool is_mixed_high_order_linear_elems);

    /// Read metis data
    /// \param file_name_base The prefix of the file name.
    void readMetisData(const std::string& file_name_base);

    /// Write mesh to METIS input file
    /// \param file_name File name with an extension of mesh.
    void writeMETIS(const std::string& file_name);

    /// Write the partitions into ASCII files
    /// \param file_name_base The prefix of the file name.
    void writeASCII(const std::string& file_name_base);

    /// Write the partitions into binary files
    /// \param file_name_base The prefix of the file name.
    void writeBinary(const std::string& file_name_base);

private:
    /// Number of partitions.
    IntegerType _npartitions;

    /// Data for all  partitions.
    std::vector<Partition> _partitions;

    /// Pointer to a mesh object.
    std::unique_ptr<MeshLib::Mesh> _mesh;

    /// Global IDs of all nodes after partitioning.
    std::vector<std::size_t> _nodes_global_ids;

    /// Partition IDs of each nodes.
    std::vector<std::size_t> _nodes_partition_ids;

    /// Flags to indicate that whether elements are processed or not.
    std::vector<bool> _elements_status;

    // Renumber the global indices of nodes,
    /// \param is_mixed_high_order_linear_elems Flag to indicate whether the
    /// elements of a mesh can be used for both linear and high order
    /// interpolation
    void renumberNodeIndices(const bool is_mixed_high_order_linear_elems);

    /*!
       Calculate the total number of integer variables of an element
       vector. Each element has three integer variables for element ID,
       element type, number of nodes of the element. Therefore
       the total number of the integers in an element vector is
        3 * vector size + sum (number of nodes of each element)
    */
    IntegerType getNumberOfIntegerVariablesOfElements(
        const std::vector<const MeshLib::Element*>& elements) const;

    /*!
         \brief Get integer variables, which are used to define an element
         \param elem            Element
         \param local_node_ids  Local node indices of a partition
         \param elem_info       A vector holds all integer variables of
                                element definitions
         \param counter         Recorder of the number of integer variables.
    */
    void getElementIntegerVariables(const MeshLib::Element& elem,
                                    const std::vector<IntegerType>& local_node_ids,
                                    std::vector<IntegerType>& elem_info,
                                    IntegerType& counter);

    void writePropertiesBinary(std::string const& file_name_base) const;

    /*!
         \brief Write the configuration data of the partition data in
                binary files.
         \param file_name_base The prefix of the file name.
         \return element 1: The numbers of all non-ghost element integer
                            variables of each partitions.
                 element 2: The numbers of all ghost element integer
                            variables of each partitions.
    */
    std::tuple<std::vector<IntegerType>, std::vector<IntegerType>>
    writeConfigDataBinary(const std::string& file_name_base);

    /*!
         \brief Write the element integer variables of all partitions
                into binary files.
         \param file_name_base      The prefix of the file name.
         \param num_elem_integers   The numbers of all non-ghost element
                                    integer variables of each partitions.
         \param num_g_elem_integers The numbers of all ghost element
                                    integer variables of each partitions.
    */
    void writeElementsBinary(const std::string& file_name_base,
                      const std::vector<IntegerType>& num_elem_integers,
                      const std::vector<IntegerType>& num_g_elem_integers);

    ///  Write the nodes of all partitions into a binary file.
    ///  \param file_name_base The prefix of the file name.
    void writeNodesBinary(const std::string& file_name_base);


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

}  // namespace MeshLib
