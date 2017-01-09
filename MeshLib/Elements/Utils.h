/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHLIB_ELEMENTS_UTILS_H_
#define MESHLIB_ELEMENTS_UTILS_H_

#include <algorithm>
#include <vector>

#include "BaseLib/makeVectorUnique.h"
#include "MeshLib/Node.h"

#include "Element.h"

namespace MeshLib
{
/// Returns a vector of node pointers containing the base nodes of the elements
/// input vector.
inline std::vector<Node*> getBaseNodes(std::vector<Element*> const& elements)
{
    std::vector<Node*> base_nodes;
    base_nodes.reserve(elements.size() * 2);  // Save some of the realloctions.

    for (auto* const e : elements)
    {
        std::copy(e->getNodes(), e->getNodes() + e->getNumberOfBaseNodes(),
                  std::back_inserter(base_nodes));
    }

    BaseLib::makeVectorUnique(base_nodes, [](Node const* a, Node* b) {
        return a->getID() < b->getID();
    });

    return base_nodes;
}

}  // namespace MeshLib

#endif  // MESHLIB_ELEMENTS_UTILS_H_
