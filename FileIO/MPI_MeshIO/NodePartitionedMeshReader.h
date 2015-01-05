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

#include "NodePartitionedMesh.h"

namespace MeshLib
{
class Node;
class Element;
}

namespace FileIO
{
/// Node data only for parallel reading.
struct NodeData
{
    long index; ///< Global node index.
    double x;
    double y;
    double z;
};

/// Class to handle reading data of partitioned mesh.
class NodePartitionedMeshReader
{
    public:

        NodePartitionedMeshReader()
        {
        }

        /*!
             \brief Create a NodePartitionedMesh object, read data to it,
                    and return a pointer to it.
             \param comm  MPI  communicator.
             \param file_name  Name of file to be read.
             \return           Pointer to Mesh object. If the creation of mesh object
                               fails, return a null pointer.
        */
        MeshLib::NodePartitionedMesh* read(MPI_Comm comm, const std::string &file_name);

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

        /*!
             \brief Set mesh nodes from a tempory array containing node data read from file.
             \param node_data  Array containing node data read from file.
             \param mesh_node  Vector of mesh nodes to be set.
        */
        void setNodes(const NodeData *node_data, std::vector<MeshLib::Node*> &mesh_node);

        /*!
             \brief Set mesh elements from a tempory array containing node data read from file.
             \param mesh_nodes        Vector of mesh nodes used to set element nodes.
             \param elem_data         Array containing element data read from file.
             \param mesh_elems        Vector of mesh elements to be set.
             \param mesh_ghost_elems  Local IDs of active element nodes.
             \param ghost             Flag of processing ghost elements.
        */
        void setElements(const std::vector<MeshLib::Node*> &mesh_nodes, const long *elem_data,
                         std::vector<MeshLib::Element*> &mesh_elems,
                         std::vector<short*> &mesh_ghost_elems, const bool ghost = false);

        /// Teminate programm due to failed to open the file or mismatch between two requested numbers
        void printMessage(const std::string & err_message, const bool for_fileopen = true);
};

} // End of namespace

#endif // end of #ifndef READ_NODE_PARTITIONED_MESH_H
