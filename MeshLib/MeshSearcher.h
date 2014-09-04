/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHSEARCHER_H_
#define MESHSEARCHER_H_

#include <vector>


namespace MeshLib
{
// forward declarations
class Mesh;

/**
 * get a list of elements connected to given nodes
 * @param msh       a mesh object
 * @param node_ids  a vector of mesh node ids
 * @return a vector of element ids which connect to the given nodes
 */
std::vector<std::size_t> getConnectedElementIDs(MeshLib::Mesh const& msh, const std::vector<std::size_t> &node_ids);

} // end namespace MeshLib

#endif //MESHSEARCHER_H_
