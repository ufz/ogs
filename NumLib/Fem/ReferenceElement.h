/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "BaseLib/DynamicSpan.h"
#include "NumLib/Fem/CoordinatesMapping/NaturalNodeCoordinates.h"

namespace NumLib
{
// Returns the coordinates as a span of dynamic size.
template <typename MeshElementType>
BaseLib::DynamicSpan<const std::array<double, 3>>
getNodeCoordsOfReferenceElement()
{
    auto const& coords =
        NumLib::NaturalCoordinates<MeshElementType>::coordinates;

    return BaseLib::DynamicSpan<const std::array<double, 3>>(coords.begin(),
                                                             coords.size());
}

template <typename MeshElementType>
class ReferenceElement
{
    static std::vector<MeshLib::Node> createReferenceElementNodes()
    {
        auto const natural_node_coordss =
            getNodeCoordsOfReferenceElement<MeshElementType>();

        return {natural_node_coordss.begin(), natural_node_coordss.end()};
    }

    static MeshElementType createElement(
        std::vector<MeshLib::Node> const& nodes)
    {
        constexpr unsigned num_nodes = MeshElementType::n_all_nodes;

        std::array<MeshLib::Node*, num_nodes> node_ptrs{};

        for (std::size_t i = 0; i < num_nodes; ++i)
        {
            node_ptrs[i] = const_cast<MeshLib::Node*>(&(nodes[i]));
        }

        return MeshElementType(node_ptrs);
    }

    std::vector<MeshLib::Node> const nodes = createReferenceElementNodes();

public:
    MeshElementType const element = createElement(nodes);
};
}  // namespace NumLib
