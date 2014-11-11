/*!
  \file NodePartitionedMeshReader.h
  \author Wenqing Wang
  \date   2014.08
  \brief  Declare a class to read node-wise partitioned mesh with MPI functions.

  \copyright
  Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/

#ifndef NODE_PARTITIONED_MESH_READER_H
#define NODE_PARTITIONED_MESH_READER_H

#include <string>
#include <fstream>

#include <mpi.h>

#include "MeshLib/NodePartitionedMesh.h"

namespace MeshLib
{
class Node;
class Element;
}

namespace FileIO
{
/// Class to handle reading data of partitioned mesh.
class NodePartitionedMeshReader
{
    public:
        NodePartitionedMeshReader() = default;

        /*!
             \brief Create a NodePartitionedMesh object, read data to it,
                    and return a pointer to it.
             \param comm            MPI Communicator.
             \param file_name_base  Name of file to be read, and it must be base name without name extension.
             \return           Pointer to Mesh object. If the creation of mesh object
                               fails, return a null pointer.
        */
        MeshLib::NodePartitionedMesh* read(MPI_Comm comm, const std::string &file_name_base);

    private:
        /// Number of all nodes of a partition.
        long _num_nodes_part;

        /// Number of regular (non-ghost) elements of a partition.
        long _num_regular_elems_part;

        /// Number of ghost elements of a partition.
        long _num_ghost_elems_part;

        /// Number of MPI processes
        int _size;

        /// _size converted to string
        std::string _size_str;

        /// Rank of compute core
        int _rank;

        /*!
             \brief Create a NodePartitionedMesh object, read binary mesh data to it,
                    and return a pointer to it.
             \param comm  MPI  Communicator.
             \param file_name  Name of file to be read, which must be a name with the
                               path to the file and without file extension.
             \return           Pointer to Mesh object.
        */
        MeshLib::NodePartitionedMesh* readBinary(MPI_Comm comm, const std::string &file_name);

        /*!
             \brief Create a NodePartitionedMesh object, read ASCII mesh data to it,
                    and return a pointer to it.
             \param comm  MPI  Communicator.
             \param file_name  Name of file to be read, which must be a name with the
                               path to the file and without file extension.
             \return           Pointer to Mesh object.
        */
        MeshLib::NodePartitionedMesh* readASCII(MPI_Comm comm, const std::string &file_name);

        /*!
             \brief Read elements data from ASCII file.
             \param ins       Input stream.
             \param elem_info Pointer to array that contains element data, which to be filled.
             \param ghost     Flag to read ghost elements.
        */
        void readElementASCII(std::ifstream &ins, long *elem_info,
                              const bool ghost = false);

        /// Node data only for parallel reading.
        struct NodeData
        {
            long index; ///< Global node index.
            double x;
            double y;
            double z;
        };

        /*! Define MPI data type, MPI_Node_ptr, for struct MeshNode for palllel reading of nodes
              \param anode        a NodeData variable.
              \param MPI_Node_ptr Defined MPI data type of struct NodeData.
        */
        void buildNodeStrucTypeMPI(NodeData *anode, MPI_Datatype *MPI_Node_ptr);
        /*!
             \brief Set mesh nodes from a tempory array containing node data read from file.
             \param node_data  Array containing node data read from file.
             \param mesh_node  Vector of mesh nodes to be set.
             \param glb_node_ids  Global IDs of nodes of a partition.
        */
        void setNodes(const NodeData *node_data, std::vector<MeshLib::Node*> &mesh_node,
                      std::vector<unsigned> &glb_node_ids);

        /*!
             \brief Set mesh elements from a tempory array containing node data read from file.
             \param mesh_nodes        Vector of mesh nodes used to set element nodes.
             \param elem_data         Array containing element data read from file.
             \param mesh_elems        Vector of mesh elements to be set.
             \param ghost             Flag of processing ghost elements.
        */
        void setElements(const std::vector<MeshLib::Node*> &mesh_nodes, const long *elem_data,
                         std::vector<MeshLib::Element*> &mesh_elems, const bool ghost = false);

        /// Terminate programm due to failed in file opening or due to mismatch between two requested numbers
        void printMessage(const std::string & err_message, const bool for_fileopen = true);
};

} // End of namespace

#endif // end of #ifndef READ_NODE_PARTITIONED_MESH_H
