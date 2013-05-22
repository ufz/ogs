/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-05-20
 * \brief  Implementation tests.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include "MathLib/LinAlg/Dense/Matrix.h"

namespace
{

inline void ASSERT_DOUBLE_ARRAY_EQ(const double* Expected, const double* Actual, size_t N, double epsilon=1.0e-8) {
    for (size_t i=0; i<N; i++) \
        ASSERT_NEAR(Expected[i], Actual[i], epsilon);
}

}


TEST(Math, CheckInterface_DenseMatrix)
{
    MathLib::Matrix<double> mat(4,4);

    ASSERT_EQ(4, mat.getNRows());
    ASSERT_EQ(4, mat.getNCols());
    mat.setZero();
    mat.setValue(0, 0, 1.0);
    mat.addValue(0, 0, 1.0);

    mat.finishAssembly();
    ASSERT_TRUE(mat.isAssembled());

    MathLib::Matrix<double> sub(2,2);
    std::vector<std::size_t> vec_pos(2);
    mat.addSubMatrix(vec_pos, vec_pos, sub);

    MathLib::Matrix<double> mat2(mat);

}

