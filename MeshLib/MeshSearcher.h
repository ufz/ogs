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
class Element;

/**
 * get a vector of elements connected to given nodes
 * @param msh       a mesh object
 * @param node_ids  a vector of mesh node ids
 * @return a vector of element ids which connect to the given nodes
 */
std::vector<std::size_t> getConnectedElementIDs(MeshLib::Mesh const& msh, const std::vector<std::size_t> &node_ids);

/**
 * get a vector of node ID connected to given elements
 * @param elements  a vector of a pointer to a mesh element object
 * @return a vector of node ID
 */
std::vector<std::size_t> getConnectedNodeIDs(const std::vector<MeshLib::Element*> &elements);

} // end namespace MeshLib

#endif //MESHSEARCHER_H_
