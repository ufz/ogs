/*!
  \file
  \author Wenqing Wang
  \date   2014.08
  \brief  Declare a class to read node-wise partitioned mesh with MPI functions.

  \copyright
  Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
 */

#pragma once

#include <iosfwd>
#include <string>
#include <vector>

#include <mpi.h>

#include "MeshLib/NodePartitionedMesh.h"
#include "MeshLib/Properties.h"
#include "MeshLib/IO/MPI_IO/PropertyVectorMetaData.h"

namespace MeshLib
{
class Node;
class Element;
class Properties;

namespace IO
{
/// Class for parallel reading of binary partitioned mesh files into a
/// NodePartitionedMesh via MPI.
class NodePartitionedMeshReader final
{
public:
    ///  \param comm   MPI communicator.
    explicit NodePartitionedMeshReader(MPI_Comm comm);

    ~NodePartitionedMeshReader();

    /*!
         \brief Create a NodePartitionedMesh object, read data to it, and return
                a pointer to it. Data files are in binary format.
         \param file_name_base  Name of file to be read, and it must be base name without name extension.
         \return                Pointer to Mesh object. If the creation of mesh object
                                fails, return a null pointer.
     */
    MeshLib::NodePartitionedMesh* read(const std::string &file_name_base);

private:
    /// Pointer to MPI communicator.
    MPI_Comm _mpi_comm;

    /// Number of processes in the communicator: _mpi_comm.
    int _mpi_comm_size;

    /// Rank of compute core.
    int _mpi_rank;

    /// MPI data type for struct NodeData.
    MPI_Datatype _mpi_node_type;

    /// Node data only for parallel reading.
    struct NodeData
    {
        std::size_t index;     ///< Global node index.
        double x;
        double y;
        double z;
    };

    /// Define MPI data type for NodeData struct.
    void registerNodeDataMpiType();

    /// A collection of integers that configure the partitioned mesh data.
    struct PartitionedMeshInfo
    {
        unsigned long nodes;  ///< 0: Number of all nodes of a partition,
        unsigned long base_nodes;        ///< 1: Number of nodes for linear
                                         ///elements of a partition,
        unsigned long regular_elements;  ///< 2: Number of non-ghost elements
                                         ///of a partition,
        unsigned long
            ghost_elements;  ///< 3: Number of ghost element of a partition,
        unsigned long active_base_nodes;  ///< 4: Number of active nodes for
                                          /// linear element of a partition,
        unsigned long
            active_nodes;  ///< 5: Number of all active nodes a partition,
        unsigned long global_base_nodes;  ///< 6: unused, previously number of
                                          /// nodes for linear element of global
                                          /// mesh,
        unsigned long global_nodes;  ///< 7: Number of all nodes of global mesh,
        unsigned long offset[5];   ///< 8~12: Offsets of positions of partitions
                                   /// in the data arrays.
        unsigned long extra_flag;  ///< 13: Reserved for extra flag.

        std::size_t size() const { return 14; }
        unsigned long* data() { return &nodes; }
    } _mesh_info;

    /*!
        \brief Create a new mesh of NodePartitionedMesh after reading and
       processing the data.
        \param mesh_name    Name assigned to the new mesh.
        \param mesh_nodes   Node data.
        \param glb_node_ids Global IDs of nodes.
        \param mesh_elems   Element data.
        \param properties Collection of PropertyVector's assigned to the mesh.
        \return Returns a pointer to a NodePartitionedMesh
     */
    MeshLib::NodePartitionedMesh* newMesh(
        std::string const& mesh_name,
        std::vector<MeshLib::Node*> const& mesh_nodes,
        std::vector<unsigned long> const& glb_node_ids,
        std::vector<MeshLib::Element*> const& mesh_elems,
        MeshLib::Properties const& properties) const;

    /*!
        Parallel reading of a binary file via MPI_File_read, reading mesh data
        head, nodes, non-ghost elements and ghost elements.

        \note In case of failure during opening of the file, an error message is
              printed.
        \note If the number of elements in container is larger than
              MPI_file_read() supports (maximum of current \c int type), an
              error is printed.

        \param filename File name containing data.
        \param offset   Displacement of the data accessible from the view.
                        see MPI_File_set_view() documentation.
        \param type     Type of data.
        \param data     A container to be filled with data. Its size is used
                        to determine how many values should be read.
        \tparam DATA    A homogeneous contaner type supporting data() and size().

        \return         True on success and false otherwise.
     */
    template <typename DATA>
    bool readDataFromFile(std::string const& filename, MPI_Offset offset,
        MPI_Datatype type, DATA& data) const;

    /*!
         \brief Create a NodePartitionedMesh object, read binary mesh data
                in the manner of parallel, and return a pointer to it.
                Four binary files have to been read in this function named as:
                file_name_base+_partitioned_msh_cfg[number of partitions].bin
                file_name_base+_partitioned_msh_nod[number of partitions].bin
                file_name_base+_partitioned_msh_ele[number of partitions].bin
                file_name_base+_partitioned_msh_ele_g[number of partitions].bin
                in which, the first file contains an array of integers for the
                PartitionMeshInfo for all partitions

                the second file contains a struct type (long, double double double) array of
                nodes information of global IDs and coordinates of all partitions.

                the third file contains a long type integer array of element information of
                material ID, element type and node IDs of each non-ghost element of all partitoions.

                the forth file contains a long type integer array of element information of
                material ID, element type and node IDs of each ghost element of all partitoions.
         \param file_name_base  Name of file to be read, which must be a name with the
                           path to the file and without file extension.
         \return           Pointer to Mesh object.
     */
    MeshLib::NodePartitionedMesh* readMesh(const std::string &file_name_base);

    MeshLib::Properties readProperties(
        const std::string& file_name_base) const;

    void readProperties(const std::string& file_name_base,
                              MeshLib::MeshItemType t,
                              MeshLib::Properties& p) const;

    void readDomainSpecificPartOfPropertyVectors(
        std::vector<std::optional<MeshLib::IO::PropertyVectorMetaData>> const&
            vec_pvmd,
        MeshLib::IO::PropertyVectorPartitionMetaData const& pvpmd,
        MeshLib::MeshItemType t,
        std::istream& is,
        MeshLib::Properties& p) const;

    template <typename T>
    void createPropertyVectorPart(
        std::istream& is, MeshLib::IO::PropertyVectorMetaData const& pvmd,
        MeshLib::IO::PropertyVectorPartitionMetaData const& pvpmd,
        MeshLib::MeshItemType t, unsigned long global_offset,
        MeshLib::Properties& p) const
    {
        MeshLib::PropertyVector<T>* pv = p.createNewPropertyVector<T>(
            pvmd.property_name, t, pvmd.number_of_components);
        pv->resize(pvpmd.number_of_tuples * pvmd.number_of_components);
        // jump to the place for reading the specific part of the
        // PropertyVector
        is.seekg(global_offset +
                 pvpmd.offset * pvmd.number_of_components * sizeof(T));
        // read the values
        unsigned long const number_of_bytes = pvmd.data_type_size_in_bytes *
                                              pvpmd.number_of_tuples *
                                              pvmd.number_of_components;
        if (!is.read(reinterpret_cast<char*>(pv->data()), number_of_bytes))
            OGS_FATAL(
                "Error in NodePartitionedMeshReader::readProperties: "
                "Could not read part {:d} of the PropertyVector.",
                _mpi_rank);
    }

    /*!
         \brief Set mesh nodes from a temporary array containing node data read from file.
         \param node_data     Vector containing node data read from file.
         \param mesh_node     Vector of mesh nodes to be set.
         \param glb_node_ids  Global IDs of nodes of a partition.
     */
    void setNodes(const std::vector<NodeData> &node_data,
        std::vector<MeshLib::Node*> &mesh_node,
        std::vector<unsigned long> &glb_node_ids) const;

    /*!
         \brief Set mesh elements from a temporary array containing node data read from file.
         \param mesh_nodes        Vector of mesh nodes used to set element nodes.
         \param elem_data         Vector containing element data read from file.
         \param mesh_elems        Vector of mesh elements to be set.
         \param ghost             Flag of processing ghost elements.
     */
    void setElements(const std::vector<MeshLib::Node*> &mesh_nodes,
        const std::vector<unsigned long> &elem_data,
        std::vector<MeshLib::Element*> &mesh_elems,
        const bool ghost = false) const;
};
}   // end namespace IO
}   // end namespace MeshLib
