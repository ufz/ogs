/*!
  \file MeshPartitioning.h
  \date   2016.05

  \brief  Declare a class to perform mesh partitioning

  \copyright
  Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/

#ifndef MESH_PARTITIONING_H_
#define MESH_PARTITIONING_H_

#include <vector>
#include <string>

#include "MeshLib/Mesh.h"

namespace MeshLib
{

/// A subdomain mesh.
class MeshPartitioning : public Mesh
{
    public:
        /// Copy constructor
        MeshPartitioning(const MeshLib::Mesh &mesh) : Mesh(mesh)
        {}

        /// Write mesh to METIS input file
        void write2METIS(const std::string& file_name);

        /// Partition by element.
        void partitionByElement()
        {
            /* to be decided */
        }

        /// Partition by node.
        void partitionByNode();

        /// Write the partitioned mesh into ASCII files.
        void writeNodePartitionedMeshASCII();

        /// Write the partitioned mesh into binary files.
        void writeNodePartitionedMeshBinary();

    private:

};

}   // namespace MeshLib

#endif // MESH_PARTITIONING_H_
