/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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

    double const aperture0 = 1e-5;
    double const penalty_aperture_cutoff = aperture0;
    bool const tension_cutoff = true;
    LinearElasticIsotropic<2> fractureModel{penalty_aperture_cutoff,
                                            tension_cutoff, mp};
    std::unique_ptr<FractureModelBase<2>::MaterialStateVariables> state(
        fractureModel.createMaterialStateVariables());

    Eigen::Vector2d const w_prev = Eigen::Vector2d::Zero();
    Eigen::Vector2d const sigma_prev = Eigen::Vector2d::Zero();
    Eigen::Vector2d w(-1e-5, -1e-5);

    // Result vectors, not initialized.
    Eigen::Vector2d sigma;
    Eigen::Matrix2d C;

    ProcessLib::SpatialPosition x;
    fractureModel.computeConstitutiveRelation(0, x, aperture0, w_prev, w,
                                              sigma_prev, sigma, C, *state);

    ASSERT_NEAR(-1e4, sigma[0], eps_sigma);
    ASSERT_NEAR(-1e6, sigma[1], eps_sigma);
    ASSERT_NEAR(1e9, C(0,0), eps_C);
    ASSERT_NEAR(0, C(0,1), eps_C);
    ASSERT_NEAR(0, C(1,0), eps_C);
    ASSERT_NEAR(1e11, C(1,1), eps_C);


    w << -1e-5, 1e-5;
    fractureModel.computeConstitutiveRelation(0, x, aperture0, w_prev, w,
                                              sigma_prev, sigma, C, *state);

    ASSERT_NEAR(0, sigma[0], eps_sigma);
    ASSERT_NEAR(0, sigma[1], eps_sigma);
    ASSERT_NEAR(0, C(0,0), eps_C);
    ASSERT_NEAR(0, C(0,1), eps_C);
    ASSERT_NEAR(0, C(1,0), eps_C);
    ASSERT_NEAR(0, C(1,1), eps_C);
}

TEST(MaterialLib_Fracture, MohrCoulomb2D_elastic)
{
    ProcessLib::ConstantParameter<double> const kn("", 1e11);
    ProcessLib::ConstantParameter<double> const ks("", 1e9);
    ProcessLib::ConstantParameter<double> const phi("", 15);
    ProcessLib::ConstantParameter<double> const psi("", 5);
    ProcessLib::ConstantParameter<double> const c("", 3e16);
    MohrCoulomb<2>::MaterialProperties const mp{kn, ks, phi, psi, c};

    double const aperture0 = 1;
    double const penalty_aperture_cutoff = aperture0;
    bool const tension_cutoff = true;
    MohrCoulomb<2> fractureModel{penalty_aperture_cutoff, tension_cutoff, mp};
    std::unique_ptr<FractureModelBase<2>::MaterialStateVariables> state(
        fractureModel.createMaterialStateVariables());

    Eigen::Vector2d const w_prev = Eigen::Vector2d::Zero();
    Eigen::Vector2d const sigma0 = Eigen::Vector2d::Zero();
    Eigen::Vector2d const sigma_prev = sigma0;
    Eigen::Vector2d w(-1e-5, -1e-5);

    // Result vectors, not initialized.
    Eigen::Vector2d sigma;
    Eigen::Matrix2d C;

    ProcessLib::SpatialPosition x;
    fractureModel.computeConstitutiveRelation(0, x, aperture0, sigma0, w_prev,
                                              w, sigma_prev, sigma, C, *state);

    ASSERT_NEAR(-1e4, sigma[0], eps_sigma);
    ASSERT_NEAR(-1e6, sigma[1], eps_sigma);
    ASSERT_NEAR(1e9, C(0, 0), eps_C);
    ASSERT_NEAR(0, C(0, 1), eps_C);
    ASSERT_NEAR(0, C(1, 0), eps_C);
    ASSERT_NEAR(1e11, C(1, 1), eps_C);

    w << -1e-5, 1e-5;
    fractureModel.computeConstitutiveRelation(0, x, aperture0, sigma0, w_prev,
                                              w, sigma_prev, sigma, C, *state);

    ASSERT_NEAR(0, sigma[0], eps_sigma);
    ASSERT_NEAR(0, sigma[1], eps_sigma);
    ASSERT_NEAR(0, C(0, 0), eps_C);
    ASSERT_NEAR(0, C(0, 1), eps_C);
    ASSERT_NEAR(0, C(1, 0), eps_C);
    ASSERT_NEAR(0, C(1, 1), eps_C);
}

TEST(MaterialLib_Fracture, MohrCoulomb2D_negative_t)
{
    ProcessLib::ConstantParameter<double> const kn("", 50e9);
    ProcessLib::ConstantParameter<double> const ks("", 20e9);
    ProcessLib::ConstantParameter<double> const phi("", 15);
    ProcessLib::ConstantParameter<double> const psi("", 5);
    ProcessLib::ConstantParameter<double> const c("", 3e6);
    MohrCoulomb<2>::MaterialProperties const mp{kn, ks, phi, psi, c};

    double const aperture0 = 1e-5;
    double const penalty_aperture_cutoff = aperture0;
    bool const tension_cutoff = true;
    MohrCoulomb<2> fractureModel{penalty_aperture_cutoff, tension_cutoff, mp};
    std::unique_ptr<FractureModelBase<2>::MaterialStateVariables> state(
        fractureModel.createMaterialStateVariables());

    Eigen::Vector2d const w_prev = Eigen::Vector2d::Zero();
    Eigen::Vector2d const sigma_prev(-3.46e6, -2e6);
    Eigen::Vector2d const w(-1.08e-5, -0.25e-5);

    // Result vectors, not initialized.
    Eigen::Vector2d sigma;
    Eigen::Matrix2d C;

    ProcessLib::SpatialPosition x;
    fractureModel.computeConstitutiveRelation(0, x, aperture0, w_prev, w,
                                              sigma_prev, sigma, C, *state);

    ASSERT_NEAR(-3.50360e6, sigma[0], eps_sigma);
    ASSERT_NEAR(-2.16271e6, sigma[1], eps_sigma);
    ASSERT_NEAR(1.10723e+09, C(0,0), eps_C);
    ASSERT_NEAR(1.26558e+10, C(0,1), eps_C);
    ASSERT_NEAR(4.13226e+09, C(1,0), eps_C);
    ASSERT_NEAR(4.72319e+10, C(1,1), eps_C);
    ASSERT_NEAR(1.06608E+5, state->getShearYieldFunctionValue(), eps_sigma);
}


TEST(MaterialLib_Fracture, MohrCoulomb2D_positive_t)
{
    ProcessLib::ConstantParameter<double> const kn("", 50e9);
    ProcessLib::ConstantParameter<double> const ks("", 20e9);
    ProcessLib::ConstantParameter<double> const phi("", 15);
    ProcessLib::ConstantParameter<double> const psi("", 5);
    ProcessLib::ConstantParameter<double> const c("", 3e6);
    MohrCoulomb<2>::MaterialProperties const mp{kn, ks, phi, psi, c};

    double const aperture0 = 1e-5;
    double const penalty_aperture_cutoff = aperture0;
    bool const tension_cutoff = true;
    MohrCoulomb<2> fractureModel{penalty_aperture_cutoff, tension_cutoff, mp};
    std::unique_ptr<FractureModelBase<2>::MaterialStateVariables> state(
        fractureModel.createMaterialStateVariables());

    Eigen::Vector2d const w_prev = Eigen::Vector2d::Zero();
    Eigen::Vector2d const sigma_prev(0, -2e6);
    Eigen::Vector2d const w(20.08e-5, -0.25e-5);

    // Result vectors, not initialized.
    Eigen::Vector2d sigma;
    Eigen::Matrix2d C;

    ProcessLib::SpatialPosition x;
    fractureModel.computeConstitutiveRelation(0, x, aperture0, w_prev, w,
                                              sigma_prev, sigma, C, *state);

    ASSERT_NEAR(253972.0, sigma[0], eps_sigma);
    ASSERT_NEAR(-2.94784e+06, sigma[1], eps_sigma);
    ASSERT_NEAR(1.10723e+09, C(0,0), eps_C);
    ASSERT_NEAR(-1.26558e+10, C(0,1), eps_C);
    ASSERT_NEAR(-4.13226e+09, C(1,0), eps_C);
    ASSERT_NEAR(4.72319e+10, C(1,1), eps_C);
    ASSERT_NEAR(446608.0, state->getShearYieldFunctionValue(), eps_sigma);
}

TEST(MaterialLib_Fracture, MohrCoulomb3D_negative_t1)
{
    ProcessLib::ConstantParameter<double> const kn("", 50e9);
    ProcessLib::ConstantParameter<double> const ks("", 20e9);
    ProcessLib::ConstantParameter<double> const phi("", 15);
    ProcessLib::ConstantParameter<double> const psi("", 5);
    ProcessLib::ConstantParameter<double> const c("", 3e6);
    MohrCoulomb<3>::MaterialProperties const mp{kn, ks, phi, psi, c};

    double const aperture0 = 1e-5;
    double const penalty_aperture_cutoff = aperture0;
    bool const tension_cutoff = true;
    MohrCoulomb<3> fractureModel{penalty_aperture_cutoff, tension_cutoff, mp};
    std::unique_ptr<FractureModelBase<3>::MaterialStateVariables> state(
        fractureModel.createMaterialStateVariables());

    Eigen::Vector3d const w_prev = Eigen::Vector3d::Zero();
    Eigen::Vector3d const sigma_prev(-3.46e6, 0, -2e6);
    Eigen::Vector3d const w(-1.08e-5, 0, -0.25e-5);

    // Result vectors, not initialized.
    Eigen::Vector3d sigma;
    Eigen::Matrix3d C;

    ProcessLib::SpatialPosition x;
    fractureModel.computeConstitutiveRelation(0, x, aperture0, w_prev, w,
                                              sigma_prev, sigma, C, *state);

    ASSERT_NEAR(-3.50360e6, sigma[0], eps_sigma);
    ASSERT_NEAR(0.0, sigma[1], eps_sigma);
    ASSERT_NEAR(-2.16271e6, sigma[2], eps_sigma);
    ASSERT_NEAR(1.10723e+09, C(0,0), eps_C);
    ASSERT_NEAR(0.0, C(0,1), eps_C);
    ASSERT_NEAR(1.26558e+10, C(0,2), eps_C);
    ASSERT_NEAR(0.0, C(1,0), eps_C);
    ASSERT_NEAR(20.e9, C(1,1), eps_C);
    ASSERT_NEAR(0.0, C(1,2), eps_C);
    ASSERT_NEAR(4.13226e+09, C(2,0), eps_C);
    ASSERT_NEAR(0.0, C(2,1), eps_C);
    ASSERT_NEAR(4.72319e+10, C(2,2), eps_C);
    ASSERT_NEAR(1.06608E+5, state->getShearYieldFunctionValue(), eps_sigma);
}

TEST(MaterialLib_Fracture, MohrCoulomb3D_positive_t1)
{
    ProcessLib::ConstantParameter<double> const kn("", 50e9);
    ProcessLib::ConstantParameter<double> const ks("", 20e9);
    ProcessLib::ConstantParameter<double> const phi("", 15);
    ProcessLib::ConstantParameter<double> const psi("", 5);
    ProcessLib::ConstantParameter<double> const c("", 3e6);
    MohrCoulomb<3>::MaterialProperties const mp{kn, ks, phi, psi, c};

    double const aperture0 = 1e-5;
    double const penalty_aperture_cutoff = aperture0;
    bool const tension_cutoff = true;
    MohrCoulomb<3> fractureModel{penalty_aperture_cutoff, tension_cutoff, mp};
    std::unique_ptr<FractureModelBase<3>::MaterialStateVariables> state(
        fractureModel.createMaterialStateVariables());

    Eigen::Vector3d const w_prev = Eigen::Vector3d::Zero();
    Eigen::Vector3d const sigma_prev(0, 0, -2e6);
    Eigen::Vector3d const w(20.08e-5, 0, -0.25e-5);

    // Result vectors, not initialized.
    Eigen::Vector3d sigma;
    Eigen::Matrix3d C;

    ProcessLib::SpatialPosition x;
    fractureModel.computeConstitutiveRelation(0, x, aperture0, w_prev, w,
                                              sigma_prev, sigma, C, *state);

    ASSERT_NEAR(253972.0, sigma[0], eps_sigma);
    ASSERT_NEAR(0.0, sigma[1], eps_sigma);
    ASSERT_NEAR(-2.94784e+06, sigma[2], eps_sigma);
    ASSERT_NEAR(1.10723e+09, C(0,0), eps_C);
    ASSERT_NEAR(0.0, C(0,1), eps_C);
    ASSERT_NEAR(-1.26558e+10, C(0,2), eps_C);
    ASSERT_NEAR(0.0, C(1,0), eps_C);
    ASSERT_NEAR(20.e9, C(1,1), eps_C);
    ASSERT_NEAR(0.0, C(1,2), eps_C);
    ASSERT_NEAR(-4.13226e+09, C(2,0), eps_C);
    ASSERT_NEAR(0.0, C(2,1), eps_C);
    ASSERT_NEAR(4.72319e+10, C(2,2), eps_C);
    ASSERT_NEAR(446608.0, state->getShearYieldFunctionValue(), eps_sigma);
}


TEST(MaterialLib_Fracture, MohrCoulomb3D_positive_t2)
{
    ProcessLib::ConstantParameter<double> const kn("", 50e9);
    ProcessLib::ConstantParameter<double> const ks("", 20e9);
    ProcessLib::ConstantParameter<double> const phi("", 15);
    ProcessLib::ConstantParameter<double> const psi("", 5);
    ProcessLib::ConstantParameter<double> const c("", 3e6);
    MohrCoulomb<3>::MaterialProperties const mp{kn, ks, phi, psi, c};

    double const aperture0 = 1e-5;
    double const penalty_aperture_cutoff = aperture0;
    bool const tension_cutoff = true;
    MohrCoulomb<3> fractureModel{penalty_aperture_cutoff, tension_cutoff, mp};
    std::unique_ptr<FractureModelBase<3>::MaterialStateVariables> state(
        fractureModel.createMaterialStateVariables());

    Eigen::Vector3d const w_prev = Eigen::Vector3d::Zero();
    Eigen::Vector3d const sigma_prev(0, 0, -2e6);
    Eigen::Vector3d const w(0, 20.08e-5, -0.25e-5);

    // Result vectors, not initialized.
    Eigen::Vector3d sigma;
    Eigen::Matrix3d C;

    ProcessLib::SpatialPosition x;
    fractureModel.computeConstitutiveRelation(0, x, aperture0, w_prev, w,
                                              sigma_prev, sigma, C, *state);

    ASSERT_NEAR(0.0, sigma[0], eps_sigma);
    ASSERT_NEAR(253972.0, sigma[1], eps_sigma);
    ASSERT_NEAR(-2.94784e+06, sigma[2], eps_sigma);
    ASSERT_NEAR(20.e9, C(0,0), eps_C);
    ASSERT_NEAR(0.0, C(0,1), eps_C);
    ASSERT_NEAR(0.0, C(0,2), eps_C);
    ASSERT_NEAR(0.0, C(1,0), eps_C);
    ASSERT_NEAR(1.10723e+09, C(1,1), eps_C);
    ASSERT_NEAR(-1.26558e+10, C(1,2), eps_C);
    ASSERT_NEAR(0.0, C(2,0), eps_C);
    ASSERT_NEAR(-4.13226e+09, C(2,1), eps_C);
    ASSERT_NEAR(4.72319e+10, C(2,2), eps_C);
    ASSERT_NEAR(446608.0, state->getShearYieldFunctionValue(), eps_sigma);
}


TEST(MaterialLib_Fracture, MohrCoulomb3D_negative_t1t2)
{
    ProcessLib::ConstantParameter<double> const kn("", 50e9);
    ProcessLib::ConstantParameter<double> const ks("", 20e9);
    ProcessLib::ConstantParameter<double> const phi("", 15);
    ProcessLib::ConstantParameter<double> const psi("", 5);
    ProcessLib::ConstantParameter<double> const c("", 3e6);
    MohrCoulomb<3>::MaterialProperties const mp{kn, ks, phi, psi, c};

    double const aperture0 = 1e-5;
    double const penalty_aperture_cutoff = aperture0;
    bool const tension_cutoff = true;
    MohrCoulomb<3> fractureModel{penalty_aperture_cutoff, tension_cutoff, mp};
    std::unique_ptr<FractureModelBase<3>::MaterialStateVariables> state(
        fractureModel.createMaterialStateVariables());

    Eigen::Vector3d const w_prev = Eigen::Vector3d::Zero();
    Eigen::Vector3d const sigma_prev(-3.46e6/std::sqrt(2), -3.46e6/std::sqrt(2), -2e6);
    Eigen::Vector3d const w(-1.08e-5/std::sqrt(2), -1.08e-5/std::sqrt(2), -0.25e-5);

    // Result vectors, not initialized.
    Eigen::Vector3d sigma;
    Eigen::Matrix3d C;

    ProcessLib::SpatialPosition x;
    fractureModel.computeConstitutiveRelation(0, x, aperture0, w_prev, w,
                                              sigma_prev, sigma, C, *state);

    ASSERT_NEAR(-3.50360e6/std::sqrt(2), sigma[0], eps_sigma);
    ASSERT_NEAR(-3.50360e6/std::sqrt(2), sigma[1], eps_sigma);
    ASSERT_NEAR(-2.16271e6, sigma[2], eps_sigma);
    ASSERT_NEAR(1.05536e+10, C(0,0), eps_C);
    ASSERT_NEAR(-9.44638e+09, C(0,1), eps_C);
    ASSERT_NEAR(8.94897e+09, C(0,2), eps_C);
    ASSERT_NEAR(-9.44638e+09, C(1,0), eps_C);
    ASSERT_NEAR(1.05536e+10, C(1,1), eps_C);
    ASSERT_NEAR(8.94897e+09, C(1,2), eps_C);
    ASSERT_NEAR(2.92195e+09, C(2,0), eps_C);
    ASSERT_NEAR(2.92195e+09, C(2,1), eps_C);
    ASSERT_NEAR(4.72319e+10, C(2,2), eps_C);
    ASSERT_NEAR(1.06608E+5, state->getShearYieldFunctionValue(), eps_sigma);
}


TEST(MaterialLib_Fracture, MohrCoulomb3D_negative_t1_positive_t2)
{
    ProcessLib::ConstantParameter<double> const kn("", 50e9);
    ProcessLib::ConstantParameter<double> const ks("", 20e9);
    ProcessLib::ConstantParameter<double> const phi("", 15);
    ProcessLib::ConstantParameter<double> const psi("", 5);
    ProcessLib::ConstantParameter<double> const c("", 3e6);
    MohrCoulomb<3>::MaterialProperties const mp{kn, ks, phi, psi, c};

    double const aperture0 = 1e-5;
    double const penalty_aperture_cutoff = aperture0;
    bool const tension_cutoff = true;
    MohrCoulomb<3> fractureModel{penalty_aperture_cutoff, tension_cutoff, mp};
    std::unique_ptr<FractureModelBase<3>::MaterialStateVariables> state(
        fractureModel.createMaterialStateVariables());

    Eigen::Vector3d const w_prev = Eigen::Vector3d::Zero();
    Eigen::Vector3d const sigma_prev(-3.46e6/std::sqrt(2), 3.46e6/std::sqrt(2), -2e6);
    Eigen::Vector3d const w(-1.08e-5/std::sqrt(2), 1.08e-5/std::sqrt(2), -0.25e-5);

    // Result vectors, not initialized.
    Eigen::Vector3d sigma;
    Eigen::Matrix3d C;

    ProcessLib::SpatialPosition x;
    fractureModel.computeConstitutiveRelation(0, x, aperture0, w_prev, w,
                                              sigma_prev, sigma, C, *state);

    ASSERT_NEAR(-3.50360e6/std::sqrt(2), sigma[0], eps_sigma);
    ASSERT_NEAR(3.50360e6/std::sqrt(2), sigma[1], eps_sigma);
    ASSERT_NEAR(-2.16271e6, sigma[2], eps_sigma);
    ASSERT_NEAR(1.05536e+10, C(0,0), eps_C);
    ASSERT_NEAR(9.44638e+09, C(0,1), eps_C);
    ASSERT_NEAR(8.94897e+09, C(0,2), eps_C);
    ASSERT_NEAR(9.44638e+09, C(1,0), eps_C);
    ASSERT_NEAR(1.05536e+10, C(1,1), eps_C);
    ASSERT_NEAR(-8.94897e+09, C(1,2), eps_C);
    ASSERT_NEAR(2.92195e+09, C(2,0), eps_C);
    ASSERT_NEAR(-2.92195e+09, C(2,1), eps_C);
    ASSERT_NEAR(4.72319e+10, C(2,2), eps_C);
    ASSERT_NEAR(1.06608E+5, state->getShearYieldFunctionValue(), eps_sigma);
}

