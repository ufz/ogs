/**
 * \author Norihiro Watanabe
 * \date   2013-09-03
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <gtest/gtest.h>

#include <algorithm>
#include <type_traits>
#include <limits>
#include <numeric>
#include <valarray>

#include "MeshLib/Elements/Elements.h"
#include "NumLib/Fem/ShapeFunction/ShapeHex20.h"
#include "NumLib/Fem/ShapeFunction/ShapeHex8.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine2.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine3.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad4.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad8.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad9.h"
#include "NumLib/Fem/ShapeFunction/ShapeTet10.h"
#include "NumLib/Fem/ShapeFunction/ShapeTet4.h"
#include "NumLib/Fem/ShapeFunction/ShapeTri3.h"
#include "NumLib/Fem/ShapeFunction/ShapeTri6.h"
#include "Tests/AutoCheckTools.h"
#include "Tests/TestTools.h"

namespace autocheck
{

template <typename ShapeFunction,
          typename Gen = randomTupleGenerator<double, ShapeFunction::DIM>>
struct NaturalPointGenerator
{
    Gen generator;  // The default interval is fixed to [-1, 1].

    using result_type = std::array<double, ShapeFunction::DIM>;

    result_type intervalMap(result_type const& tuple) const
    {
        if (std::is_same<ShapeFunction, NumLib::ShapeLine2>::value ||
            std::is_same<ShapeFunction, NumLib::ShapeLine3>::value ||
            std::is_same<ShapeFunction, NumLib::ShapeQuad4>::value ||
            std::is_same<ShapeFunction, NumLib::ShapeQuad8>::value ||
            std::is_same<ShapeFunction, NumLib::ShapeQuad9>::value ||
            std::is_same<ShapeFunction, NumLib::ShapeHex8>::value ||
            std::is_same<ShapeFunction, NumLib::ShapeHex20>::value)
            return tuple;
        if (std::is_same<ShapeFunction, NumLib::ShapeTri3>::value ||
            std::is_same<ShapeFunction, NumLib::ShapeTri6>::value ||
            std::is_same<ShapeFunction, NumLib::ShapeTet4>::value ||
            std::is_same<ShapeFunction, NumLib::ShapeTet10>::value)
        {
            // Map square (x, y) \in [-1, 1]^2 to a triangle such that x,y
            // \in [0, 1], and x+y \in [0, 2/2].
            // Same for three-d case with x + y + z \in [0, 3/2].
            result_type mapped_tuple;
            std::transform(std::begin(tuple), std::end(tuple),
                           std::begin(mapped_tuple),
                           [](double const& v) { return (v + 1.) / 2.; });

            if (std::accumulate(std::begin(mapped_tuple),
                                std::end(mapped_tuple),
                                0.) > ShapeFunction::DIM / 2.)
                std::transform(std::begin(mapped_tuple), std::end(mapped_tuple),
                               std::begin(mapped_tuple),
                               [](double const& v) { return 1. - v; });
            return mapped_tuple;
        }
    }

    result_type operator()(std::size_t /*size*/ = 0)
    {
        return intervalMap(fix(1, generator)());
    }
};

}

namespace ac = autocheck;
using namespace NumLib;

template <typename ShapeFunction>
struct ShapeFunctionTest : public ::testing::Test
{
    ac::NaturalPointGenerator<ShapeFunction> natural_point_generator;
    ac::gtest_reporter gtest_reporter;
};

using ShapeFunctionTestTypes =
    ::testing::Types<ShapeLine2, ShapeLine3, ShapeTri3, ShapeTri6, ShapeTet4,
                     ShapeTet10, ShapeQuad4, ShapeQuad8, ShapeQuad9, ShapeHex8,
                     ShapeHex20>;

TYPED_TEST_CASE(ShapeFunctionTest, ShapeFunctionTestTypes);

// TypeParam is the type of the ShapeFunction.
// Access private members via this pointer or TestFixture:: for types
TYPED_TEST(ShapeFunctionTest, PartitionOfUnity)
{
    auto isPartitionOfUnity =
        [](std::array<double, TypeParam::DIM>& natural_coordinates_point)
        -> bool {

        // compute shape functions
        std::array<double, TypeParam::NPOINTS> N;
        TypeParam::computeShapeFunction(natural_coordinates_point, N);
        double const sum = std::accumulate(std::begin(N), std::end(N), 0.);

        return std::abs(1. - sum) < 1e-15;
    };

    ac::check<std::array<double, TypeParam::DIM>>(
        isPartitionOfUnity,
        1000,
        ac::make_arbitrary(this->natural_point_generator),
        this->gtest_reporter);
}

TYPED_TEST(ShapeFunctionTest, SumOfGradientsIsZero)
{
    auto isSumOfGradientsZero =
        [](std::array<double, TypeParam::DIM>& natural_coordinates_point)
        -> bool {

        // compute shape functions
        std::array<double, TypeParam::DIM * TypeParam::NPOINTS> dNdr;
        TypeParam::computeGradShapeFunction(natural_coordinates_point, dNdr);
        double const sum =
            std::accumulate(std::begin(dNdr), std::end(dNdr), 0.);

        return std::abs(sum) <= 5e-15;
    };

    ac::check<std::array<double, TypeParam::DIM>>(
        isSumOfGradientsZero,
        1000,
        ac::make_arbitrary(this->natural_point_generator),
        this->gtest_reporter);
}

TEST(NumLib, FemShapeQuad4)
{
    static const double eps = std::numeric_limits<double>::epsilon();
    static const unsigned NNodes = 4;
    static const unsigned dim = 2;
    std::valarray<double> r(dim);
    std::valarray<double> N(NNodes);
    std::valarray<double> dN(NNodes*dim);

    // check N, dN at specific location
    {
        r = .5; // r = (0,5, 0.5)
        ShapeQuad4::computeShapeFunction(r, N);
        ShapeQuad4::computeGradShapeFunction(r, dN);
        double exp_N[]= {0.5625, 0.1875, 0.0625, 0.1875};
        double exp_dN[]= {0.375, -0.375, -0.125, 0.125, 0.375, 0.125, -0.125, -0.375};
        ASSERT_ARRAY_NEAR(exp_N, N, N.size(), eps);
        ASSERT_ARRAY_NEAR(exp_dN, dN, dN.size(), eps);
    }

    std::valarray<double> exp_N(NNodes);
    // check N_i(r_j)= {i==j: 1, i!=j: 0}
    for (unsigned i=0; i<NNodes; i++)
    {
        r[0] = (i==0 || i==3) ? 1 : -1;
        r[1] = i<2 ? 1 : -1;
        exp_N = .0;
        exp_N[i] = 1.0;
        ShapeQuad4::computeShapeFunction(r, N);
        ASSERT_ARRAY_NEAR(exp_N, N, NNodes, eps);
    }

    // check sum_i[N_i(r)] = 1
    const double dist = 0.5;
    for (unsigned i=0; i<NNodes; i++)
    {
        r[0] = (i==0 || i==3) ? dist : -dist;
        r[1] = i<2 ? dist : -dist;
        ShapeQuad4::computeShapeFunction(r, N);
        ASSERT_NEAR(1.0, N.sum(), eps);
//        for(auto v: N)
//            std::cout << v << " ";
//        std::cout << "\n";
    }

}


