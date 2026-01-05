// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>

#include "NumLib/Fem/CoordinatesMapping/NaturalNodeCoordinates.h"

namespace NumLib
{
template <typename MeshElementType>
class ReferenceElement
{
    static std::vector<MeshLib::Node> createReferenceElementNodes()
    {
        constexpr auto natural_node_coordss =
            NumLib::NaturalCoordinates<MeshElementType>::coordinates;

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
