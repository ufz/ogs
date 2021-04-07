/**
 * \file
 * \author Karsten Rink
 * \date   2014-03-25
 * \brief  Definition of Duplicate functions
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

namespace MeshLib
{
class Mesh;
class Node;
class Element;

/// Creates a deep copy of a Node vector
std::vector<Node*> copyNodeVector(const std::vector<Node*>& nodes);

/** Creates a deep copy of an element vector using the given Node vector.
 * @param elements The element vector that should be duplicated.
 * @param new_nodes The new node vector used for the duplicated element vector.
 * @param node_id_map An optional mapping from the nodes the 'old' elements
 * based on to the new nodes. This should be consistent with the original node
 * vector.
 * @return A deep copy of the elements vector using the new nodes vector.
 */
std::vector<Element*> copyElementVector(
    std::vector<Element*> const& elements,
    std::vector<Node*> const& new_nodes,
    std::vector<std::size_t> const* const node_id_map = nullptr);

/// Copies an element without change, using the nodes vector from the result
/// mesh and an optional mapping from the nodes the 'old' elements.
Element* copyElement(Element const* const element,
                     const std::vector<Node*>& nodes,
                     std::vector<std::size_t> const* const id_map = nullptr);

/// Clones a vector of elements using the Element::clone() function.
std::vector<Element*> cloneElements(std::vector<Element*> const& elements);
}  // end namespace MeshLib
