/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */


#include <limits>

#include <gtest/gtest.h>

#include <Eigen/Eigen>

#include "BaseLib/ConfigTree.h"

#include "MaterialLib/FractureModels/LinearElasticIsotropic.h"
#include "MaterialLib/FractureModels/MohrCoulomb.h"

#include "ProcessLib/Parameter/ConstantParameter.h"

using namespace MaterialLib::Fracture;

static const double eps_sigma = 1e6*1e-5;
static const double eps_C = 1e10*1e-5;

TEST(MaterialLib_Fracture, LinearElasticIsotropic)
{
    ProcessLib::ConstantParameter<double> const kn("", 1e11);
    ProcessLib::ConstantParameter<double> const ks("", 1e9);
    LinearElasticIsotropic<2>::MaterialProperties const mp{kn, ks};

    LinearElasticIsotropic<2> fractureModel{mp};

    Eigen::Vector2d w_prev, w, sigma_prev, sigma;
    Eigen::Matrix2d C;

    ProcessLib::SpatialPosition x;
    w << -1e-5, -1e-5;
    fractureModel.computeConstitutiveRelation(0, x, w_prev, w, sigma_prev, sigma, C);

    ASSERT_NEAR(-1e4, sigma[0], eps_sigma);
    ASSERT_NEAR(-1e6, sigma[1], eps_sigma);
    ASSERT_NEAR(1e9, C(0,0), eps_C);
    ASSERT_NEAR(0, C(0,1), eps_C);
    ASSERT_NEAR(0, C(1,0), eps_C);
    ASSERT_NEAR(1e11, C(1,1), eps_C);


    w << -1e-5, 1e-5;
    fractureModel.computeConstitutiveRelation(0, x, w_prev, w, sigma_prev, sigma, C);

    ASSERT_NEAR(0, sigma[0], eps_sigma);
    ASSERT_NEAR(0, sigma[1], eps_sigma);
    ASSERT_NEAR(0, C(0,0), eps_C);
    ASSERT_NEAR(0, C(0,1), eps_C);
    ASSERT_NEAR(0, C(1,0), eps_C);
    ASSERT_NEAR(0, C(1,1), eps_C);
}

TEST(MaterialLib_Fracture, MohrCoulomb)
{
    ProcessLib::ConstantParameter<double> const kn("", 50e9);
    ProcessLib::ConstantParameter<double> const ks("", 20e9);
    ProcessLib::ConstantParameter<double> const phi("", 15);
    ProcessLib::ConstantParameter<double> const psi("", 5);
    ProcessLib::ConstantParameter<double> const c("", 3e6);
    MohrCoulomb<2>::MaterialProperties const mp{kn, ks, phi, psi, c};

    MohrCoulomb<2> fractureModel{mp};

    Eigen::Vector2d w_prev, w, sigma_prev, sigma;
    Eigen::Matrix2d C;

    ProcessLib::SpatialPosition x;
    sigma_prev << -3.46e6, -2e6;
    w << -1.08e-5, -0.25e-5;
    fractureModel.computeConstitutiveRelation(0, x, w_prev, w, sigma_prev, sigma, C);

    ASSERT_NEAR(-3.50360e6, sigma[0], eps_sigma);
    ASSERT_NEAR(-2.16271e6, sigma[1], eps_sigma);
    ASSERT_NEAR(1.10723e+09, C(0,0), eps_C);
    ASSERT_NEAR(1.26558e+10, C(0,1), eps_C);
    ASSERT_NEAR(4.13226e+09, C(1,0), eps_C);
    ASSERT_NEAR(4.72319e+10, C(1,1), eps_C);
}

