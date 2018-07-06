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

#include <fstream>
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
    /// Non ghost elements
    std::vector<const MeshLib::Element*> regular_elements;
    std::vector<const MeshLib::Element*> ghost_elements;
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
    void getElementIntegerVariables(
        const MeshLib::Element& elem,
        const std::vector<IntegerType>& local_node_ids,
        std::vector<IntegerType>& elem_info,
        IntegerType& counter);

    void writeNodePropertiesBinary(std::string const& file_name_base) const;
    void writeCellPropertiesBinary(std::string const& file_name_base) const;

    /// 1 copy pointers to nodes belonging to the partition part_id
    /// 2 collect non-linear element nodes belonging to the partition part_id in
    /// extra_nodes
    void findNonGhostNodesInPartition(
        std::size_t const part_id,
        const bool is_mixed_high_order_linear_elems,
        std::vector<MeshLib::Node*>& extra_nodes);

    /// 1 find elements belonging to the partition part_id:
    /// fills vector partition.regular_elements
    /// 2 find ghost elements belonging to the partition part_id
    /// fills vector partition.ghost_elements
    void findElementsInPartition(std::size_t const part_id);

    /// Prerequisite: the ghost elements has to be found (using
    /// findElementsInPartition).
    /// Finds ghost nodes and non-linear element ghost nodes by walking over
    /// ghost elements.
    void findGhostNodesInPartition(std::size_t const part_id,
                                   const bool is_mixed_high_order_linear_elems,
                                   std::vector<MeshLib::Node*>& extra_nodes);

    void splitOfHigherOrderNode(std::vector<MeshLib::Node*> const& nodes,
                                bool const is_mixed_high_order_linear_elems,
                                unsigned const node_id,
                                std::vector<MeshLib::Node*>& base_nodes,
                                std::vector<MeshLib::Node*>& extra_nodes);

    void processPartition(std::size_t const part_id,
                          const bool is_mixed_high_order_linear_elems);

    void processNodeProperties();
    void processCellProperties();

    template <typename T>
    bool copyNodePropertyVector(std::string const& name,
                                std::size_t const total_number_of_tuples)
    {
        auto const& original_properties(_mesh->getProperties());
        if (!original_properties.existsPropertyVector<T>(name))
            return false;

        auto const& pv(original_properties.getPropertyVector<T>(name));
        auto partitioned_pv =
            _partitioned_properties.createNewPropertyVector<T>(
                name, pv->getMeshItemType(), pv->getNumberOfComponents());
        partitioned_pv->resize(total_number_of_tuples *
                               pv->getNumberOfComponents());
        std::size_t position_offset(0);
        for (auto p : _partitions)
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
    bool copyCellPropertyVector(std::string const& name,
                                std::size_t const total_number_of_tuples)
    {
        auto const& original_properties(_mesh->getProperties());
        if (!original_properties.existsPropertyVector<T>(name))
            return false;

        auto const& pv(original_properties.getPropertyVector<T>(name));
        auto partitioned_pv =
            _partitioned_properties.createNewPropertyVector<T>(
                name, pv->getMeshItemType(), pv->getNumberOfComponents());
        partitioned_pv->resize(total_number_of_tuples *
                               pv->getNumberOfComponents());
        std::size_t position_offset(0);
        for (auto const& p : _partitions)
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

    template <typename T>
    void writePropertyVectorValuesBinary(
        std::ostream& os, MeshLib::PropertyVector<T> const& pv) const
    {
        std::size_t number_of_components(pv.getNumberOfComponents());
        std::size_t number_of_tuples(pv.getNumberOfTuples());
        std::vector<T> property_vector_buffer;
        property_vector_buffer.resize(number_of_tuples * number_of_components);
        for (std::size_t i = 0; i < pv.getNumberOfTuples(); ++i)
        {
            for (std::size_t c(0); c < number_of_components; ++c)
                property_vector_buffer[i * number_of_components + c] =
                    pv.getComponent(i, c);
        }
        os.write(reinterpret_cast<char*>(property_vector_buffer.data()),
                 number_of_components * number_of_tuples * sizeof(T));
    }

    template <typename T>
    bool writePropertyVectorBinary(std::string const& name,
                                   std::ostream& out_val,
                                   std::ostream& out_meta) const
    {
        if (!_partitioned_properties.existsPropertyVector<T>(name))
            return false;

        MeshLib::IO::PropertyVectorMetaData pvmd;
        pvmd.property_name = name;
        auto* pv = _partitioned_properties.getPropertyVector<T>(name);
        pvmd.fillPropertyVectorMetaDataTypeInfo<T>();
        pvmd.number_of_components = pv->getNumberOfComponents();
        pvmd.number_of_tuples = pv->getNumberOfTuples();
        writePropertyVectorValuesBinary(out_val, *pv);
        MeshLib::IO::writePropertyVectorMetaDataBinary(out_meta, pvmd);
        return true;
    }

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
    void writeElementsBinary(
        const std::string& file_name_base,
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

}  // namespace ApplicationUtils
