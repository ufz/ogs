/*!
  \file MeshPartitioning.cpp
  \date   2016.05

  \brief  Define the members of class MeshPartitioning

  \copyright
  Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license

*/

#include "MeshPartitioning.h"

#include <fstream>

#include "MeshLib/Elements/Element.h"

namespace MeshLib
{
void MeshPartitioning :: write2METIS(const std::string& file_name)
{
    std::ofstream os(file_name.data(), std::ios::trunc);
    os << _elements.size() <<" \n";
    for (const auto elem : _elements)
    {
        for (unsigned j=0; j<elem->getNNodes(); j++)
        {
            os << elem->getNodeIndex(j) + 1 <<" ";
        }
        os << std::endl;
    }
}

void MeshPartitioning :: partitionByNode()
{

}

void MeshPartitioning :: writeNodePartitionedMeshASCII()
{
}

void MeshPartitioning :: writeNodePartitionedMeshBinary()
{
}

}   // namespace MeshLib

