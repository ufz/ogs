/*!
  \file NodeWiseMeshPartitioner.h
  \date   2016.05

  \brief  Declare a class to perform node wise mesh partitioning

  \copyright
  Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/

#ifndef NODE_WISE_MESH_PARTITIONER_H_
#define NODE_WISE_MESH_PARTITIONER_H_

#include <memory>
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
    typedef long PetscInt;

public:
    /*!
     * \param file_name_base The prefix of the file name.
     * \param num_partitions Number of partitions,
     * \param mesh           Pointer to a mesh object.
     */
    NodeWiseMeshPartitioner(const PetscInt num_partitions,
                            std::unique_ptr<MeshLib::Mesh>& mesh)
        : _npartitions(num_partitions),
          _partitions(num_partitions),
          _mesh(std::move(mesh)),
          _nodes_global_ids(_mesh->getNumberOfNodes()),
          _nodes_partition_ids(_mesh->getNumberOfNodes()),
          _elements_status(_mesh->getNumberOfElements(), false)
    {
    }
    ~NodeWiseMeshPartitioner() = default;

    /// Partition by node.
    /// \param is_mixed_hl_elem Flag to indicate whether the elements of
    /// a mesh can be used for both linear and high order interpolation
    void partitionByMETIS(const bool is_mixed_hl_elem);

    void resetGlobalNodeIndecis();

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

    /// Write the global mesh into a VTU file
    /// \param file_name_base The prefix of the file name.
    void writeGlobalMeshVTU(const std::string& file_name_base);

private:
    /// Number of partitions.
    PetscInt _npartitions;

    /// Data for all  partitions.
    std::vector<Partition> _partitions;

    /// Pointer to a mesh object.
    std::unique_ptr<MeshLib::Mesh> _mesh;

    /// Global IDs of all nodes after partitioning.
    std::vector<std::size_t> _nodes_global_ids;

    /// Partition IDs of each nodes.
    std::vector<std::size_t> _nodes_partition_ids;

    /// Flags to indicate the status of all elements.
    std::vector<bool> _elements_status;

    void renumberNodeIndecies();

    /*!
       Calculate the totoal number of integer variables of an element
       vector.
           Each element has three integer variables for material ID,
       element type, number of nodes of the element. Therefore
       the total number of the integers in an element vector is
        3 * vector size + sum (number of nodes of each element)
     \param is_ghost Flag to indicate ghost elements or not
   */
    PetscInt getNumberOfIntegerVariablesOfElements(
        const std::vector<const MeshLib::Element*>& elements) const;

    /*!
         \brief get integer variables, which are used to define an element
         \param elem            Element
         \param local_node_ids  Local node indicies of a partition
         \param elem_info       An vector holds all integer variables of element
       definitions
         \param counter         Recorder of the number of integer variables.
    */
    void getElementIntegerVariables(const MeshLib::Element& elem,
                                    const std::vector<PetscInt>& local_node_ids,
                                    std::vector<PetscInt>& elem_info,
                                    PetscInt& counter);

    /*!
        \brief Write local indicies of element nodes to a ASCII file
        \param os              Output stream
        \param elem            Element
        \param local_node_ids  Local node indicies of a partition
    */
    void writeLocalElementNodeIndicies(
        std::ostream& os,
        const MeshLib::Element& elem,
        const std::vector<PetscInt>& local_node_ids);
};

}  // namespace MeshLib

#endif  // NODE_WISE_MESH_PARTITIONER_H_
