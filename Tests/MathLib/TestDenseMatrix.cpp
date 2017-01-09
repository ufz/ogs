/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "gtest/gtest.h"

#include "MathLib/LinAlg/Dense/DenseMatrix.h"
#include "Tests/TestTools.h"

using namespace MathLib;

TEST(MathLib, DenseMatrixTransposeInPlace)
{
    const double eps(std::numeric_limits<double>::epsilon());

    // square matrix
    DenseMatrix<double> m1(3,3);
    unsigned cnt = 0;
    for (unsigned i=0; i<m1.getNumberOfRows(); i++)
        for (unsigned j=0; j<m1.getNumberOfColumns(); j++)
            m1(i,j) = cnt++;
    double expected_m1[] = {0, 1, 2, 3, 4, 5, 6, 7, 8};
    ASSERT_EQ(3u, m1.getNumberOfRows());
    ASSERT_EQ(3u, m1.getNumberOfColumns());
    ASSERT_ARRAY_NEAR(expected_m1, m1.data(), m1.size(), eps);
    m1.transposeInPlace();
    ASSERT_EQ(3u, m1.getNumberOfRows());
    ASSERT_EQ(3u, m1.getNumberOfColumns());
    double expected_m1t[] = {0, 3, 6, 1, 4, 7, 2, 5, 8};
    ASSERT_ARRAY_NEAR(expected_m1t, m1.data(), m1.size(), eps);

    // non-square matrix 1
    DenseMatrix<double> m2(2,3);
    cnt = 0;
    for (unsigned i=0; i<m2.getNumberOfRows(); i++)
        for (unsigned j=0; j<m2.getNumberOfColumns(); j++)
            m2(i,j) = cnt++;
    ASSERT_EQ(2u, m2.getNumberOfRows());
    ASSERT_EQ(3u, m2.getNumberOfColumns());
    double expected_m2[] = {0, 1, 2, 3, 4, 5};
    ASSERT_ARRAY_NEAR(expected_m2, m2.data(), m2.size(), eps);
    m2.transposeInPlace();
    ASSERT_EQ(3u, m2.getNumberOfRows());
    ASSERT_EQ(2u, m2.getNumberOfColumns());
    double expected_m2t[] = {0, 3, 1, 4, 2, 5};
    ASSERT_ARRAY_NEAR(expected_m2t, m2.data(), m2.size(), eps);

    // non-square matrix 2
    DenseMatrix<double> m3(3,2);
    cnt = 0;
    for (unsigned i=0; i<m3.getNumberOfRows(); i++)
        for (unsigned j=0; j<m3.getNumberOfColumns(); j++)
            m3(i,j) = cnt++;
    ASSERT_EQ(3u, m3.getNumberOfRows());
    ASSERT_EQ(2u, m3.getNumberOfColumns());
    double expected_m3[] = {0, 1, 2, 3, 4, 5};
    ASSERT_ARRAY_NEAR(expected_m3, m3.data(), m3.size(), eps);
    m3.transposeInPlace();
    ASSERT_EQ(2u, m3.getNumberOfRows());
    ASSERT_EQ(3u, m3.getNumberOfColumns());
    double expected_m3t[] = {0, 2, 4, 1, 3, 5};
    ASSERT_ARRAY_NEAR(expected_m3t, m3.data(), m3.size(), eps);
}
