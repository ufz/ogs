/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>

#ifdef OGS_USE_EIGEN
#include <Eigen/Eigen>

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> EigenMatrixd;

TEST(NumLib, TestEigenMultiplicationWithAuto)
{
    std::cout << "# Eigen version: " << EIGEN_WORLD_VERSION << "." << EIGEN_MAJOR_VERSION << "." << EIGEN_MINOR_VERSION << std::endl;

    EigenMatrixd c(2,1);
    c << 0, 1;
    EigenMatrixd m(1,1);
    m << 0.5;
    EigenMatrixd r(1,2);
    r << -0.5, 0.5;

    // Expression 1 : c*(m*r)
    auto exp1_auto = c * (m * r);
    EigenMatrixd exp1_auto_mat = exp1_auto;
    EigenMatrixd exp1_mat = c * (m * r);
//    std::cout << "== One step multiplications ==" << std::endl;
//    std::cout << "c*(m*r) -> auto =\n" << exp1_auto << std::endl;
//    std::cout << "c*(m*r) auto -> Mat =\n" << exp1_auto_mat << std::endl;
//    std::cout << "c*(m*r) -> Mat =\n" << exp1_mat << std::endl;
//    std::cout << "~~ c*(m*r) -> auto =\n" << exp1_auto << std::endl;

    ASSERT_TRUE(exp1_auto != exp1_mat);
    ASSERT_TRUE(exp1_auto_mat != exp1_mat);

    // Expression 2: v=m*r -> c*v
    auto m_r = m * r;
    auto exp2_auto = c * m_r;
    EigenMatrixd exp2_auto_mat = exp2_auto;
    EigenMatrixd exp2_mat = c * m_r;
//    std::cout << "\n== Two step multiplications ==" << std::endl;
//    std::cout << "c*m_r -> auto =\n" << exp2_auto << std::endl;
//    std::cout << "~~ c*(m*r) -> auto =\n" << exp1_auto << std::endl;
//    std::cout << "auto -> mat =\n" << exp2_auto_mat << std::endl;
//    std::cout << "~~ c*(m*r) -> auto =\n" << exp2_auto_mat << std::endl;
//    std::cout << "c*m_r -> Mat =\n" << exp2_mat << std::endl;
//    std::cout << "~~ c*(m*r) -> auto =\n" << exp1_auto << std::endl;

    ASSERT_TRUE(exp1_auto == exp1_mat);
    ASSERT_TRUE(exp2_auto == exp2_mat);
    ASSERT_TRUE(exp2_auto_mat == exp2_mat);

}

#endif // OGS_USE_EIGEN
