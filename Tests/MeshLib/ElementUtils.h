/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <gtest/gtest.h>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Point.h"

namespace ElementUtils
{
// Returns the number of faces of a mesh element.
//
// Here, faces of a d dimensional element are always considered d-1 dimensional.
// This differs from other places in OGS.
inline std::size_t getNumberOfFaces(MeshLib::Element const& e)
{
    switch (e.getDimension())
    {
        case 3:
            return e.getNumberOfFaces();
        case 2:
            return e.getNumberOfEdges();
        case 1:
            return 2;
    }

    return 0;
}

// Returns a d-1 dimensional face of a d dimensional mesh element.
//
// Here, faces of a d dimensional element are always considered d-1 dimensional.
// This differs from other places in OGS.
inline std::unique_ptr<MeshLib::Element const> getFace(
    MeshLib::Element const& bulk_element, unsigned face_id)
{
    auto const dim = bulk_element.getDimension();

    switch (dim)
    {
        case 3:
            return std::unique_ptr<MeshLib::Element const>{
                bulk_element.getFace(face_id)};
        case 2:
            return std::unique_ptr<MeshLib::Element const>{
                bulk_element.getEdge(face_id)};
        case 1:
        {
            auto* node =
                const_cast<MeshLib::Node*>(bulk_element.getNode(face_id));
            return std::make_unique<MeshLib::Point>(std::array{node});
        }
    }

    OGS_FATAL("Unsupported element dimension: {}", dim);
}

template <typename BulkElementRule>
unsigned getLocalBulkNodeIndexOfFaceNode(
    MeshLib::TemplateElement<BulkElementRule> const& bulk_element,
    unsigned const face_id,
    unsigned const i_face_node)
{
    EXPECT_GT(getNumberOfFaces(bulk_element), face_id);
    EXPECT_GT(getFace(bulk_element, face_id)->getNumberOfNodes(), i_face_node);

    auto constexpr dim = MeshLib::TemplateElement<BulkElementRule>::dimension;

    if constexpr (dim == 0)
    {
        return 0;
    }
    if constexpr (dim == 1)
    {
        EXPECT_EQ(0, i_face_node)
            << "Faces of line elements have only one node.";
        return BulkElementRule::edge_nodes[0][face_id];
    }
    if constexpr (dim == 2)
    {
        return BulkElementRule::edge_nodes[face_id][i_face_node];
    }
    if constexpr (dim == 3)
    {
        return BulkElementRule::face_nodes[face_id][i_face_node];
    }
}
}  // namespace ElementUtils
