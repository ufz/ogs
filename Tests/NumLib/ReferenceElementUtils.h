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

#include <gtest/gtest.h>

#include "BaseLib/DynamicSpan.h"
#include "MeshLib/Elements/Elements.h"
#include "NumLib/Fem/CoordinatesMapping/NaturalNodeCoordinates.h"
#include "Tests/Utils.h"

namespace ReferenceElementUtils
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

// Returns the coordinates as a span of dynamic size.
BaseLib::DynamicSpan<const std::array<double, 3>>
getNodeCoordsOfReferenceElement(MeshLib::CellType const cell_type);

template <typename MeshElementType>
class ReferenceElement
{
    static std::vector<MeshLib::Node> createReferenceElementNodes()
    {
        auto const natural_node_coordss =
            getNodeCoordsOfReferenceElement<MeshElementType>();

        std::vector<MeshLib::Node> nodes;
        nodes.reserve(natural_node_coordss.size());

        for (auto const& coords : natural_node_coordss)
        {
            nodes.emplace_back(coords);
        }

        return nodes;
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

std::shared_ptr<MeshLib::Element const> getReferenceElement(
    MeshLib::CellType const cell_type);

std::vector<std::array<double, 3>> getCoordsInReferenceElementForTest(
    MeshLib::Element const& element);

}  // namespace ReferenceElementUtils
