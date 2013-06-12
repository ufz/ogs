/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-05-15
 * \brief  Implementation tests of Matrix classes.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#include "MathLib/LinAlg/Dense/DenseMatrix.h"
#include "MathLib/LinAlg/Dense/GlobalDenseMatrix.h"
#ifdef USE_LIS
#include "MathLib/LinAlg/Lis/LisMatrix.h"
#endif
#include "MathLib/LinAlg/FinishMatrixAssembly.h"
#include "MathLib/LinAlg/IsMatrixAssembled.h"

namespace
{
template <class T_MATRIX>
void checkInterface(T_MATRIX &m)
{
    ASSERT_EQ(10u, m.getNRows());
    ASSERT_EQ(10u, m.getNCols());
    ASSERT_EQ(0u,  m.getRangeBegin());
    ASSERT_EQ(10u, m.getRangeEnd());

    m.setValue(0, 0, 1.0);
    m.addValue(0, 0, 1.0);
    m.setZero();

    MathLib::DenseMatrix<double> local_m(2,2, 1.0);
    std::vector<std::size_t> vec_pos(2);
    vec_pos[0] = 1;
    vec_pos[1] = 3;
    m.addSubMatrix(vec_pos, vec_pos, local_m);

    finishMatrixAssembly(m);
    ASSERT_TRUE(isMatrixAssembled(m));
}
} // end namespace

TEST(Math, CheckInterface_GlobalDenseMatrix)
{
	MathLib::GlobalDenseMatrix<double> m(10, 10);
	checkInterface(m);
}

#ifdef USE_LIS
TEST(Math, CheckInterface_LisMatrix)
{
	MathLib::LisMatrix m(10);
	checkInterface(m);
}
#endif
