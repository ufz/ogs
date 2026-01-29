// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>

#include "MeshLib/Elements/Pyramid.h"
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

        // natural and real node coordinates coincide on most reference elements
        std::vector<MeshLib::Node> real_node_coords{
            natural_node_coordss.begin(), natural_node_coordss.end()};

        if constexpr (std::is_same_v<MeshElementType, MeshLib::Pyramid13>)
        {
            // on pyramid reference elements some natural and real coordinates
            // differ
            real_node_coords[9] = {-0.5, -0.5, 0};
            real_node_coords[10] = {0.5, -0.5, 0};
            real_node_coords[11] = {0.5, 0.5, 0};
            real_node_coords[12] = {-0.5, 0.5, 0};
        }

        return real_node_coords;
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
