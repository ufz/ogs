/**
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "NumLib/Fem/CoordinatesMapping/NaturalNodeCoordinates.h"

#include <gtest/gtest.h>

#include <vector>

#include "MeshLib/ElementCoordinatesMappingLocal.h"
#include "NumLib/Fem/CoordinatesMapping/ShapeMatrices.h"
#include "NumLib/Fem/FiniteElement/C0IsoparametricElements.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/Integration/GaussLegendreIntegrationPolicy.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

using Node = MeshLib::Node;

template <typename ShapeFunction, int GlobalDim>
bool test(MeshLib::Element const& element)
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;

    // Test if the cast is possible.
    bool const dynamic_cast_failed =
        dynamic_cast<typename ShapeFunction::MeshElement const*>(&element) ==
        nullptr;
    EXPECT_FALSE(dynamic_cast_failed);
    if (dynamic_cast_failed)
    {
        return false;
    }

    // For each node evaluate the shape functions at natural coordinates and
    // test if only the corresponding shape function has value 1 and all others
    // must return 0.
    int const number_nodes = element.getNumberOfNodes();
    std::vector<MathLib::Point3d> points;
    points.reserve(number_nodes);
    for (int n = 0; n < number_nodes; ++n)
    {
        // Evaluate shape matrices at natural coordinates.
        // Compute only N, because for pyramid the detJ becomes zero at the
        // apex, and we only use N in the following test anyway.
        points.emplace_back(
            NumLib::NaturalCoordinates<
                typename ShapeFunction::MeshElement>::coordinates[n]);
    }
    auto const shape_matrices =
        NumLib::computeShapeMatrices<ShapeFunction, ShapeMatricesType,
                                     GlobalDim, NumLib::ShapeMatrixType::N>(
            element, false /*is_axially_symmetric*/, points);

    bool result = true;
    for (int n = 0; n < number_nodes; ++n)
    {
        auto const& N = shape_matrices[n].N;
        for (int p = 0; p < static_cast<int>(ShapeFunction::NPOINTS); ++p)
        {
            if (p == n)
            {
                if (N[p] != 1)
                {
                    EXPECT_EQ(1, N[p]) << "for n = " << n << ", p = " << p
                                       << " and dimension " << GlobalDim;
                    result = false;
                }
            }
            else
            {
                if (N[p] != 0)
                {
                    EXPECT_EQ(0, N[p]) << "for n = " << n << ", p = " << p
                                       << " and dimension " << GlobalDim;
                    result = false;
                }
            }
        }
    }
    return result;
}

TEST(NumLibNaturalCoordinates, Line2)
{
    using ShapeFunction = NumLib::ShapeLine2;
    using MeshElement = ShapeFunction::MeshElement;

    auto const& xs = NumLib::NaturalCoordinates<MeshElement>::coordinates;
    std::array<Node, MeshElement::n_all_nodes> ns{{Node{xs[0]}, Node{xs[1]}}};
    std::array<Node*, MeshElement::n_all_nodes> nodes{{&ns[0], &ns[1]}};
    MeshElement element(nodes);

    bool const test_result1 = test<ShapeFunction, 1>(element);
    ASSERT_TRUE(test_result1) << "Test failed for dimension 1.";

    bool const test_result2 = test<ShapeFunction, 2>(element);
    ASSERT_TRUE(test_result2) << "Test failed for dimension 2.";

    bool const test_result3 = test<ShapeFunction, 3>(element);
    ASSERT_TRUE(test_result3) << "Test failed for dimension 3.";
}

TEST(NumLibNaturalCoordinates, Line3)
{
    using ShapeFunction = NumLib::ShapeLine3;
    using MeshElement = ShapeFunction::MeshElement;

    auto const& xs = NumLib::NaturalCoordinates<MeshElement>::coordinates;
    std::array<Node, MeshElement::n_all_nodes> ns{
        {Node{xs[0]}, Node{xs[1]}, Node{xs[2]}}};
    std::array<Node*, MeshElement::n_all_nodes> nodes{{&ns[0], &ns[1], &ns[2]}};
    MeshElement element(nodes);

    bool const test_result1 = test<ShapeFunction, 1>(element);
    ASSERT_TRUE(test_result1) << "Test failed for dimension 1.";

    bool const test_result2 = test<ShapeFunction, 2>(element);
    ASSERT_TRUE(test_result2) << "Test failed for dimension 2.";

    bool const test_result3 = test<ShapeFunction, 3>(element);
    ASSERT_TRUE(test_result3) << "Test failed for dimension 3.";
}

TEST(NumLibNaturalCoordinates, Tri3)
{
    using ShapeFunction = NumLib::ShapeTri3;
    using MeshElement = ShapeFunction::MeshElement;

    auto const& xs = NumLib::NaturalCoordinates<MeshElement>::coordinates;
    std::array<Node, MeshElement::n_all_nodes> ns{
        {Node{xs[0]}, Node{xs[1]}, Node{xs[2]}}};
    std::array<Node*, MeshElement::n_all_nodes> nodes{{&ns[0], &ns[1], &ns[2]}};
    MeshElement element(nodes);

    bool const test_result2 = test<ShapeFunction, 2>(element);
    ASSERT_TRUE(test_result2) << "Test failed for dimension 2.";

    bool const test_result3 = test<ShapeFunction, 3>(element);
    ASSERT_TRUE(test_result3) << "Test failed for dimension 3.";
}

TEST(NumLibNaturalCoordinates, Tri6)
{
    using ShapeFunction = NumLib::ShapeTri6;
    using MeshElement = ShapeFunction::MeshElement;

    auto const& xs = NumLib::NaturalCoordinates<MeshElement>::coordinates;
    std::array<Node, MeshElement::n_all_nodes> ns{{Node{xs[0]}, Node{xs[1]},
                                                   Node{xs[2]}, Node{xs[3]},
                                                   Node{xs[4]}, Node{xs[5]}}};
    std::array<Node*, MeshElement::n_all_nodes> nodes{
        {&ns[0], &ns[1], &ns[2], &ns[3], &ns[4], &ns[5]}};
    MeshElement element(nodes);

    bool const test_result2 = test<ShapeFunction, 2>(element);
    ASSERT_TRUE(test_result2) << "Test failed for dimension 2.";

    bool const test_result3 = test<ShapeFunction, 3>(element);
    ASSERT_TRUE(test_result3) << "Test failed for dimension 3.";
}

TEST(NumLibNaturalCoordinates, Quad4)
{
    using ShapeFunction = NumLib::ShapeQuad4;
    using MeshElement = ShapeFunction::MeshElement;

    auto const& xs = NumLib::NaturalCoordinates<MeshElement>::coordinates;
    std::array<Node, MeshElement::n_all_nodes> ns{
        {Node{xs[0]}, Node{xs[1]}, Node{xs[2]}, Node{xs[3]}}};
    std::array<Node*, MeshElement::n_all_nodes> nodes{
        {&ns[0], &ns[1], &ns[2], &ns[3]}};
    MeshElement element(nodes);

    bool const test_result2 = test<ShapeFunction, 2>(element);
    ASSERT_TRUE(test_result2) << "Test failed for dimension 2.";

    bool const test_result3 = test<ShapeFunction, 3>(element);
    ASSERT_TRUE(test_result3) << "Test failed for dimension 3.";
}

TEST(NumLibNaturalCoordinates, Quad8)
{
    using ShapeFunction = NumLib::ShapeQuad8;
    using MeshElement = ShapeFunction::MeshElement;

    auto const& xs = NumLib::NaturalCoordinates<MeshElement>::coordinates;
    std::array<Node, MeshElement::n_all_nodes> ns{
        {Node{xs[0]}, Node{xs[1]}, Node{xs[2]}, Node{xs[3]}, Node{xs[4]},
         Node{xs[5]}, Node{xs[6]}, Node{xs[7]}}};
    std::array<Node*, MeshElement::n_all_nodes> nodes{
        {&ns[0], &ns[1], &ns[2], &ns[3], &ns[4], &ns[5], &ns[6], &ns[7]}};
    MeshElement element(nodes);

    bool const test_result2 = test<ShapeFunction, 2>(element);
    ASSERT_TRUE(test_result2) << "Test failed for dimension 2.";

    bool const test_result3 = test<ShapeFunction, 3>(element);
    ASSERT_TRUE(test_result3) << "Test failed for dimension 3.";
}

TEST(NumLibNaturalCoordinates, Quad9)
{
    using ShapeFunction = NumLib::ShapeQuad9;
    using MeshElement = ShapeFunction::MeshElement;

    auto const& xs = NumLib::NaturalCoordinates<MeshElement>::coordinates;
    std::array<Node, MeshElement::n_all_nodes> ns{
        {Node{xs[0]}, Node{xs[1]}, Node{xs[2]}, Node{xs[3]}, Node{xs[4]},
         Node{xs[5]}, Node{xs[6]}, Node{xs[7]}, Node{xs[8]}}};
    std::array<Node*, MeshElement::n_all_nodes> nodes{{&ns[0], &ns[1], &ns[2],
                                                       &ns[3], &ns[4], &ns[5],
                                                       &ns[6], &ns[7], &ns[8]}};
    MeshElement element(nodes);

    bool const test_result2 = test<ShapeFunction, 2>(element);
    ASSERT_TRUE(test_result2) << "Test failed for dimension 2.";

    bool const test_result3 = test<ShapeFunction, 3>(element);
    ASSERT_TRUE(test_result3) << "Test failed for dimension 3.";
}

TEST(NumLibNaturalCoordinates, Tet4)
{
    using ShapeFunction = NumLib::ShapeTet4;
    using MeshElement = ShapeFunction::MeshElement;

    auto const& xs = NumLib::NaturalCoordinates<MeshElement>::coordinates;
    std::array<Node, MeshElement::n_all_nodes> ns{
        {Node{xs[0]}, Node{xs[1]}, Node{xs[2]}, Node{xs[3]}}};
    std::array<Node*, MeshElement::n_all_nodes> nodes{
        {&ns[0], &ns[1], &ns[2], &ns[3]}};
    MeshElement element(nodes);

    bool const test_result3 = test<ShapeFunction, 3>(element);
    ASSERT_TRUE(test_result3) << "Test failed for dimension 3.";
}

TEST(NumLibNaturalCoordinates, Tet10)
{
    using ShapeFunction = NumLib::ShapeTet10;
    using MeshElement = ShapeFunction::MeshElement;

    auto const& xs = NumLib::NaturalCoordinates<MeshElement>::coordinates;
    std::array<Node, MeshElement::n_all_nodes> ns{
        {Node{xs[0]}, Node{xs[1]}, Node{xs[2]}, Node{xs[3]}, Node{xs[4]},
         Node{xs[5]}, Node{xs[6]}, Node{xs[7]}, Node{xs[8]}, Node{xs[9]}}};
    std::array<Node*, MeshElement::n_all_nodes> nodes{
        {&ns[0], &ns[1], &ns[2], &ns[3], &ns[4], &ns[5], &ns[6], &ns[7], &ns[8],
         &ns[9]}};
    MeshElement element(nodes);

    bool const test_result3 = test<ShapeFunction, 3>(element);
    ASSERT_TRUE(test_result3) << "Test failed for dimension 3.";
}

TEST(NumLibNaturalCoordinates, Prism6)
{
    using ShapeFunction = NumLib::ShapePrism6;
    using MeshElement = ShapeFunction::MeshElement;

    auto const& xs = NumLib::NaturalCoordinates<MeshElement>::coordinates;
    std::array<Node, MeshElement::n_all_nodes> ns{{Node{xs[0]}, Node{xs[1]},
                                                   Node{xs[2]}, Node{xs[3]},
                                                   Node{xs[4]}, Node{xs[5]}}};
    std::array<Node*, MeshElement::n_all_nodes> nodes{
        {&ns[0], &ns[1], &ns[2], &ns[3], &ns[4], &ns[5]}};
    MeshElement element(nodes);

    bool const test_result3 = test<ShapeFunction, 3>(element);
    ASSERT_TRUE(test_result3) << "Test failed for dimension 3.";
}

TEST(NumLibNaturalCoordinates, Prism15)
{
    using ShapeFunction = NumLib::ShapePrism15;
    using MeshElement = ShapeFunction::MeshElement;

    auto const& xs = NumLib::NaturalCoordinates<MeshElement>::coordinates;
    std::array<Node, MeshElement::n_all_nodes> ns{
        {Node{xs[0]}, Node{xs[1]}, Node{xs[2]}, Node{xs[3]}, Node{xs[4]},
         Node{xs[5]}, Node{xs[6]}, Node{xs[7]}, Node{xs[8]}, Node{xs[9]},
         Node{xs[10]}, Node{xs[11]}, Node{xs[12]}, Node{xs[13]}, Node{xs[14]}}};
    std::array<Node*, MeshElement::n_all_nodes> nodes{
        {&ns[0], &ns[1], &ns[2], &ns[3], &ns[4], &ns[5], &ns[6], &ns[7], &ns[8],
         &ns[9], &ns[10], &ns[11], &ns[12], &ns[13], &ns[14]}};
    MeshElement element(nodes);

    bool const test_result3 = test<ShapeFunction, 3>(element);
    ASSERT_TRUE(test_result3) << "Test failed for dimension 3.";
}

TEST(NumLibNaturalCoordinates, Pyramid5)
{
    using ShapeFunction = NumLib::ShapePyra5;
    using MeshElement = ShapeFunction::MeshElement;

    auto const& xs = NumLib::NaturalCoordinates<MeshElement>::coordinates;
    std::array<Node, MeshElement::n_all_nodes> ns{
        {Node{xs[0]}, Node{xs[1]}, Node{xs[2]}, Node{xs[3]}, Node{xs[4]}}};
    std::array<Node*, MeshElement::n_all_nodes> nodes{
        {&ns[0], &ns[1], &ns[2], &ns[3], &ns[4]}};
    MeshElement element(nodes);

    bool const test_result3 = test<ShapeFunction, 3>(element);
    ASSERT_TRUE(test_result3) << "Test failed for dimension 3.";
}

TEST(NumLibNaturalCoordinates, Pyramid13)
{
    using ShapeFunction = NumLib::ShapePyra13;
    using MeshElement = ShapeFunction::MeshElement;

    auto const& xs = NumLib::NaturalCoordinates<MeshElement>::coordinates;
    std::array<Node, MeshElement::n_all_nodes> ns{
        {Node{xs[0]}, Node{xs[1]}, Node{xs[2]}, Node{xs[3]}, Node{xs[4]},
         Node{xs[5]}, Node{xs[6]}, Node{xs[7]}, Node{xs[8]}, Node{xs[9]},
         Node{xs[10]}, Node{xs[11]}, Node{xs[12]}}};

    std::array<Node*, MeshElement::n_all_nodes> nodes{
        {&ns[0], &ns[1], &ns[2], &ns[3], &ns[4], &ns[5], &ns[6], &ns[7], &ns[8],
         &ns[9], &ns[10], &ns[11], &ns[12]}};
    MeshElement element(nodes);

    bool const test_result3 = test<ShapeFunction, 3>(element);
    ASSERT_TRUE(test_result3) << "Test failed for dimension 3.";
}

TEST(NumLibNaturalCoordinates, Hex8)
{
    using ShapeFunction = NumLib::ShapeHex8;
    using MeshElement = ShapeFunction::MeshElement;

    auto const& xs = NumLib::NaturalCoordinates<MeshElement>::coordinates;
    std::array<Node, MeshElement::n_all_nodes> ns{
        {Node{xs[0]}, Node{xs[1]}, Node{xs[2]}, Node{xs[3]}, Node{xs[4]},
         Node{xs[5]}, Node{xs[6]}, Node{xs[7]}}};
    std::array<Node*, MeshElement::n_all_nodes> nodes{
        {&ns[0], &ns[1], &ns[2], &ns[3], &ns[4], &ns[5], &ns[6], &ns[7]}};
    MeshElement element(nodes);

    bool const test_result3 = test<ShapeFunction, 3>(element);
    ASSERT_TRUE(test_result3) << "Test failed for dimension 3.";
}

TEST(NumLibNaturalCoordinates, Hex20)
{
    using ShapeFunction = NumLib::ShapeHex20;
    using MeshElement = ShapeFunction::MeshElement;

    auto const& xs = NumLib::NaturalCoordinates<MeshElement>::coordinates;
    std::array<Node, MeshElement::n_all_nodes> ns{
        {Node{xs[0]},  Node{xs[1]},  Node{xs[2]},  Node{xs[3]},  Node{xs[4]},
         Node{xs[5]},  Node{xs[6]},  Node{xs[7]},  Node{xs[8]},  Node{xs[9]},
         Node{xs[10]}, Node{xs[11]}, Node{xs[12]}, Node{xs[13]}, Node{xs[14]},
         Node{xs[15]}, Node{xs[16]}, Node{xs[17]}, Node{xs[18]}, Node{xs[19]}}};

    std::array<Node*, MeshElement::n_all_nodes> nodes{
        {&ns[0],  &ns[1],  &ns[2],  &ns[3],  &ns[4],  &ns[5],  &ns[6],
         &ns[7],  &ns[8],  &ns[9],  &ns[10], &ns[11], &ns[12], &ns[13],
         &ns[14], &ns[15], &ns[16], &ns[17], &ns[18], &ns[19]}};

    MeshElement element(nodes);

    bool const test_result3 = test<ShapeFunction, 3>(element);
    ASSERT_TRUE(test_result3) << "Test failed for dimension 3.";
}
