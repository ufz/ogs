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
#include <fstream>

#include "MeshLib/Mesh.h"

namespace MeshLib
{

enum class ElementType : unsigned {LINE2, QUAD4, HEX8, TRI3, TET4, PRISM6, PYRAMID5, INVALID};

/// A subdomain mesh.
class MeshPartitioning : public Mesh
{
        typedef long MyInt; // for PetscInt
    public:
        /// Copy constructor
        MeshPartitioning(const MeshLib::Mesh &mesh) : Mesh(mesh)
        {}

        /// Write mesh to METIS input file
        void write2METIS(const std::string& file_name);

        /// Partition by element.
        void partitionByElementMETIS(const std::string& /*file_name_base*/, const unsigned /*npartitions*/)
        {
            /* to be decided */
        }

        /// Partition by node.
        void partitionByNodeMETIS(const std::string& file_name_base, const unsigned npartitions,
                                  const bool output_binary);

    private:
        ElementType getElementType(const Element& elem);

        /*!
             \brief get integer variables, which are used to define an element
             \param elem            Element
             \param local_node_ids  Local node indicies of a partition
             \param elem_info       An vector holds all integer variables of element definitions
             \param counter         Recorder of the number of integer variables.
        */
        void getElementIntegerVariables(const Element& elem,
                                        const std::vector<unsigned>& local_node_ids,
                                        std::vector<MyInt>& elem_info,
                                        MyInt& counter);

        /*!
            \brief Write local indicies of element nodes to a ASCII file
            \param os              Output stream
            \param elem            Element
            \param local_node_ids  Local node indicies of a partition
        */
        void writeLocalElementNodeIndicies(std::ostream& os, const Element& elem,
                                           const std::vector<unsigned>& local_node_ids);

        struct NodeStruct
        {
            MyInt id;
            double x;
            double y;
            double z;
        };
};

}   // namespace MeshLib

#endif // MESH_PARTITIONING_H_
