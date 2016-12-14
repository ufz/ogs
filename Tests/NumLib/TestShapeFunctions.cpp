/**
 * \author Norihiro Watanabe
 * \date   2013-09-03
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <gtest/gtest.h>

#include <limits>
#include <valarray>
#include <algorithm>

#include "NumLib/Fem/ShapeFunction/ShapeQuad4.h"

#include "Tests/TestTools.h"

using namespace NumLib;

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


