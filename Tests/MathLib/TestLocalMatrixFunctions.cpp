/**
 * @copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <gtest/gtest.h>

#ifdef OGS_USE_EIGEN
#include <Eigen/Eigen>
#endif

#include "MathLib/LinAlg/MatrixTools.h"

#include "Tests/TestTools.h"

#ifdef OGS_USE_EIGEN
TEST(MathLib, LocalMatrixDeterminantInverse_Eigen)
{
    Eigen::Matrix3d fMat, fInv;
    fMat << 1, 2, 3,
            0, 1, 4,
            5, 6, 0;
    double fMat_det = MathLib::determinant(fMat);
    MathLib::inverse(fMat, fMat_det, fInv);

    Eigen::MatrixXd dMat(3,3), dInv(3,3);
    dMat = fMat;
    double dMat_det = MathLib::determinant(dMat);
    MathLib::inverse(dMat, dMat_det, dInv);

    ASSERT_NEAR(fMat_det, dMat_det, std::numeric_limits<double>::epsilon());
    ASSERT_ARRAY_NEAR(fInv.data(), dInv.data(), fInv.size(), std::numeric_limits<double>::epsilon());
}
#endif
