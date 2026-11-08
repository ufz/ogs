// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <memory>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/transform.hpp>
#include <stdexcept>
#include <vector>

#include "MeshLib/Elements/Hex.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Prism.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Node.h"
#include "MeshToolsLib/MeshEditing/ReverseNodeOrder.h"

namespace
{
std::vector<std::unique_ptr<MeshLib::Node>> makeNodes(
    std::vector<std::array<double, 3>> const& coords)
{
    return coords |
           ranges::views::transform(
               [](auto const& c)
               { return std::make_unique<MeshLib::Node>(c); }) |
           ranges::to<std::vector>;
}

// Nodes of the unit square in the node order OGS6 expects, i.e., the Jacobian
// determinant is positive.
std::vector<std::unique_ptr<MeshLib::Node>> makeSquareNodes()
{
    return makeNodes(
        {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {1.0, 1.0, 0.0}, {0.0, 1.0, 0.0}});
}

std::vector<MeshLib::Node*> asPointers(
    std::vector<std::unique_ptr<MeshLib::Node>> const& nodes)
{
    return nodes |
           ranges::views::transform([](auto const& n) { return n.get(); }) |
           ranges::to<std::vector>;
}
}  // namespace

// An element whose node order is already correct is neither touched nor
// counted, so the caller's count is the number of reordered elements.
TEST(ReverseNodeOrder, LeavesCorrectlyOrderedQuadUntouched)
{
    auto const nodes = makeSquareNodes();
    auto const pointers = asPointers(nodes);

    std::array<MeshLib::Node*, 4> const quad_nodes{pointers[0], pointers[1],
                                                   pointers[2], pointers[3]};
    MeshLib::Quad quad(quad_nodes);
    std::vector<MeshLib::Element*> const elements{&quad};

    EXPECT_EQ(MeshToolsLib::reverseNodeOrder(elements, 2, false, true), 0u);

    for (std::size_t i = 0; i < quad_nodes.size(); ++i)
    {
        EXPECT_EQ(quad.getNode(i), quad_nodes[i]);
    }
}

// The reversed quad is put into the OGS6 order by swapping its nodes 0 and 2
// and counted once.
TEST(ReverseNodeOrder, RestoresReversedQuadOrderAndCountsIt)
{
    auto const nodes = makeSquareNodes();
    auto const pointers = asPointers(nodes);

    std::array<MeshLib::Node*, 4> const quad_nodes{pointers[0], pointers[3],
                                                   pointers[2], pointers[1]};
    MeshLib::Quad quad(quad_nodes);
    std::vector<MeshLib::Element*> const elements{&quad};

    EXPECT_EQ(MeshToolsLib::reverseNodeOrder(elements, 2, false, true), 1u);

    std::array<MeshLib::Node*, 4> const expected{pointers[2], pointers[3],
                                                 pointers[0], pointers[1]};
    for (std::size_t i = 0; i < expected.size(); ++i)
    {
        EXPECT_EQ(quad.getNode(i), expected[i]);
    }
}

// forced additionally verifies the Jacobian determinant of the reordered
// element, which a corrected quad passes.
TEST(ReverseNodeOrder, ForcedCheckAcceptsReorderedQuad)
{
    auto const nodes = makeSquareNodes();
    auto const pointers = asPointers(nodes);

    std::array<MeshLib::Node*, 4> const quad_nodes{pointers[0], pointers[3],
                                                   pointers[2], pointers[1]};
    MeshLib::Quad quad(quad_nodes);
    std::vector<MeshLib::Element*> const elements{&quad};

    EXPECT_EQ(MeshToolsLib::reverseNodeOrder(elements, 2, true, true), 1u);
}

// Only elements with a wrong node order are counted.
TEST(ReverseNodeOrder, CountsOnlyTheReorderedElements)
{
    auto const nodes = makeSquareNodes();
    auto const pointers = asPointers(nodes);

    std::array<MeshLib::Node*, 4> const correct_nodes{pointers[0], pointers[1],
                                                      pointers[2], pointers[3]};
    MeshLib::Quad correct(correct_nodes);

    std::array<MeshLib::Node*, 4> const reversed_nodes{
        pointers[0], pointers[3], pointers[2], pointers[1]};
    MeshLib::Quad reversed(reversed_nodes);

    std::vector<MeshLib::Element*> const elements{&correct, &reversed};

    EXPECT_EQ(MeshToolsLib::reverseNodeOrder(elements, 2, false, false), 1u);
}

// The one-dimensional case uses the same contract.
TEST(ReverseNodeOrder, RestoresReversedLineOrderAndCountsIt)
{
    MeshLib::Node a(0.0, 0.0, 0.0);
    MeshLib::Node b(1.0, 0.0, 0.0);

    std::array<MeshLib::Node*, 2> const line_nodes{&b, &a};
    MeshLib::Line line(line_nodes);
    std::vector<MeshLib::Element*> const elements{&line};

    EXPECT_EQ(MeshToolsLib::reverseNodeOrder(elements, 1, true, true), 1u);
    EXPECT_EQ(line.getNode(0), &a);
    EXPECT_EQ(line.getNode(1), &b);
}

// A warped hex whose Jacobian determinant is still negative after the
// reordering terminates via OGS_FATAL, as documented in the header.
TEST(ReverseNodeOrder, ForcedCheckRejectsStillInvertedHex)
{
    auto const nodes = makeNodes({{0.0, 2.0, 2.0},
                                  {1.0, 0.0, 0.0},
                                  {2.0, 2.0, 1.0},
                                  {1.0, 1.0, 0.0},
                                  {0.0, 1.0, 0.0},
                                  {2.0, 2.0, 0.0},
                                  {2.0, 2.0, 2.0},
                                  {0.0, 0.0, 0.0}});
    auto const pointers = asPointers(nodes);

    std::array<MeshLib::Node*, 8> const hex_nodes{
        pointers[0], pointers[1], pointers[2], pointers[3],
        pointers[4], pointers[5], pointers[6], pointers[7]};
    MeshLib::Hex hex(hex_nodes);
    std::vector<MeshLib::Element*> const elements{&hex};

    EXPECT_THROW(MeshToolsLib::reverseNodeOrder(elements, 3, true, false),
                 std::runtime_error);
}

// A warped prism whose reordering changes its volume terminates via
// OGS_FATAL, as documented in the header.
TEST(ReverseNodeOrder, VolumeCheckRejectsChangedVolume)
{
    auto const nodes = makeNodes({{1.0, 0.0, 0.0},
                                  {2.0, 2.0, 1.0},
                                  {2.0, 2.0, 0.0},
                                  {0.0, 0.0, 1.0},
                                  {0.0, 2.0, 1.0},
                                  {1.0, 2.0, 1.0}});
    auto const pointers = asPointers(nodes);

    std::array<MeshLib::Node*, 6> const prism_nodes{pointers[0], pointers[1],
                                                    pointers[2], pointers[3],
                                                    pointers[4], pointers[5]};
    MeshLib::Prism prism(prism_nodes);
    std::vector<MeshLib::Element*> const elements{&prism};

    EXPECT_THROW(MeshToolsLib::reverseNodeOrder(elements, 3, false, true),
                 std::runtime_error);
}
