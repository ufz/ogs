// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <numbers>
#include <range/v3/all.hpp>
#include <string>
#include <tuple>
#include <vector>

#include "CreateTestConstitutiveRelation.h"
#include "MaterialLib/SolidModels/CreateLinearElasticIsotropic.h"
#include "MaterialLib/SolidModels/CreateLinearElasticOrthotropic.h"
#include "MaterialLib/SolidModels/CreateLinearElasticTransverseIsotropic.h"
#include "MathLib/KelvinVector.h"
#include "ParameterLib/ConstantParameter.h"
#include "ParameterLib/CoordinateSystem.h"

constexpr char xml_ti[] =
    "<constitutive_relation> "
    "    <type>LinearElasticTransverseIsotropic</type>"
    "    <youngs_modulus_i>E_i</youngs_modulus_i>"
    "    <youngs_modulus_a>E_a</youngs_modulus_a>"
    "    <poissons_ratio_ii>nu_i</poissons_ratio_ii>"
    "    <poissons_ratio_ia>nu_ia</poissons_ratio_ia>"
    "    <shear_modulus_ia>G_a</shear_modulus_ia>"
    "</constitutive_relation> ";

constexpr char xml_orth[] =
    "<constitutive_relation> "
    "    <type>LinearElasticOrthotropic</type>"
    "    <youngs_moduli>E0</youngs_moduli>"
    "    <shear_moduli>G0</shear_moduli>"
    "    <poissons_ratios>nu0</poissons_ratios>"
    "</constitutive_relation> ";

constexpr char xml_iso[] =
    "<constitutive_relation> "
    "    <type>LinearElasticIsotropic</type>"
    "    <youngs_modulus>E</youngs_modulus>"
    "    <poissons_ratio>nu</poissons_ratio>"
    "</constitutive_relation> ";

using Tensor4R =
    std::array<std::array<std::array<std::array<double, 3>, 3>, 3>, 3>;

constexpr double phi = 35 * std::numbers::pi / 180.0;

std::vector<std::unique_ptr<ParameterLib::ParameterBase>>
setParametersForLinearElasticTransverseIsotropic(double const E_i,
                                                 double const E_a,
                                                 double const nu_i,
                                                 double const nu_ia,
                                                 double const G_a)
{
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;
    parameters.push_back(
        std::make_unique<ParameterLib::ConstantParameter<double>>("E_i", E_i));
    parameters.push_back(
        std::make_unique<ParameterLib::ConstantParameter<double>>("E_a", E_a));
    parameters.push_back(
        std::make_unique<ParameterLib::ConstantParameter<double>>("nu_i",
                                                                  nu_i));
    parameters.push_back(
        std::make_unique<ParameterLib::ConstantParameter<double>>("nu_ia",
                                                                  nu_ia));
    parameters.push_back(
        std::make_unique<ParameterLib::ConstantParameter<double>>("G_a", G_a));
    return parameters;
}

std::vector<std::unique_ptr<ParameterLib::ParameterBase>>
setParametersForLinearElasticOrthotropic(std::vector<double> const& E0,
                                         std::vector<double> const& nu0,
                                         std::vector<double> const& G0)
{
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;
    parameters.push_back(
        std::make_unique<ParameterLib::ConstantParameter<double>>("E0", E0));
    parameters.push_back(
        std::make_unique<ParameterLib::ConstantParameter<double>>("nu0", nu0));
    parameters.push_back(
        std::make_unique<ParameterLib::ConstantParameter<double>>("G0", G0));
    return parameters;
}

std::vector<std::unique_ptr<ParameterLib::ParameterBase>>
setParametersForLinearElasticIsotropic(double const E, double const nu)
{
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;
    parameters.push_back(
        std::make_unique<ParameterLib::ConstantParameter<double>>("E", E));
    parameters.push_back(
        std::make_unique<ParameterLib::ConstantParameter<double>>("nu", nu));
    return parameters;
}

class LinearElasticTransverseIsotropic : public ::testing::Test
{
public:
    template <int Dim>
    void compareWithElasticOrthotropic(
        std::optional<ParameterLib::CoordinateSystem> const& coordinate_system)
    {
        // Create a LinearElasticTransverseIsotropic instance:
        auto const parameters_ti =
            setParametersForLinearElasticTransverseIsotropic(8.0e9, 4.0e9, 0.35,
                                                             0.25, 1.2e9);
        auto const elastic_model_transverse_isotropy =
            Tests::createTestConstitutiveRelation<
                MaterialLib::Solids::LinearElasticTransverseIsotropic<Dim>>(
                xml_ti, parameters_ti, coordinate_system, false,
                MaterialLib::Solids::createLinearElasticTransverseIsotropic<
                    Dim>);

        ParameterLib::SpatialPosition const pos;
        auto const C = elastic_model_transverse_isotropy->getElasticTensor(
            t_, pos, T_ref_);

        // Create a LinearElasticOrthotropic instance:
        auto const parameters_orth =
            getParametersForLinearElasticOrthotropic(Dim);

        auto const elastic_model_orthotropic =
            Tests::createTestConstitutiveRelation<
                MaterialLib::Solids::LinearElasticOrthotropic<Dim>>(
                xml_orth, parameters_orth, coordinate_system, false,
                MaterialLib::Solids::createLinearElasticOrthotropic<Dim>);
        auto const C_oth =
            elastic_model_orthotropic->getElasticTensor(t_, pos, T_ref_);

        // Compare the elastic tensor obtained by two models, respectively:
        ASSERT_LE((C - C_oth).norm() / C.norm(), 1e-14);

        MathLib::KelvinVector::KelvinMatrixType<Dim> Cel;

        // Compare the bulk modulus obtained by the transverse isotropic elastic
        // model with the expected value:
        double const k_ti =
            elastic_model_transverse_isotropy->getBulkModulus(t_, pos, &Cel);
        ASSERT_LE(std::abs(k_ti - 4301075268.8172045), 1e-14)
            << "Calculated bulk modulus by the transverse isotropy model: "
            << k_ti << "\n"
            << "Expected Bulk modulus: 4301075268.8172045";

        // The definitions of the bulk modulus by LinearElasticOrthotropic and
        // by LinearElasticTransverseIsotropic are different. Therefore, the
        // comparison of the bulk modulus values by the models are not
        // presented.
    }

    template <int Dimension>
    void compareWithLinearElasticIsotropic(
        std::optional<ParameterLib::CoordinateSystem> const& coordinate_system)
    {
        // Create an isotropic elastic model by using
        // LinearElasticTransverseIsotropic:
        double const E = 8.0e9;
        double const nu = 0.25;
        double const Ga = 0.5 * E / (1 + nu);
        auto const parameters_ti =
            setParametersForLinearElasticTransverseIsotropic(E, E, nu, nu, Ga);
        auto const elastic_model_transverse_isotropy =
            Tests::createTestConstitutiveRelation<
                MaterialLib::Solids::LinearElasticTransverseIsotropic<
                    Dimension>>(
                xml_ti, parameters_ti, coordinate_system, false,
                MaterialLib::Solids::createLinearElasticTransverseIsotropic<
                    Dimension>);
        ParameterLib::SpatialPosition const pos;
        auto const C_a = elastic_model_transverse_isotropy->getElasticTensor(
            t_, pos, T_ref_);

        // Create an isotropic elastic model by using LinearElasticIsotropic:
        auto parameters_iso = setParametersForLinearElasticIsotropic(E, nu);
        auto const elastic_model_linear_elastic_isotropic =
            Tests::createTestConstitutiveRelation<
                MaterialLib::Solids::LinearElasticIsotropic<Dimension>>(
                xml_iso, parameters_iso, false,
                MaterialLib::Solids::createLinearElasticIsotropic<Dimension>);
        auto const C_i =
            elastic_model_linear_elastic_isotropic->getElasticTensor(t_, pos,
                                                                     T_ref_);
        // Compare the elastic tensor obtained by two models, respectively:
        ASSERT_LE((C_a - C_i).norm() / C_a.norm(), 1e-14);

        MathLib::KelvinVector::KelvinMatrixType<Dimension> Cel;

        // Compare the bulk modulus obtained by obtained by two models,
        // respectively:
        double const k_ti =
            elastic_model_transverse_isotropy->getBulkModulus(t_, pos, &Cel);
        double const k_i =
            elastic_model_linear_elastic_isotropic->getBulkModulus(t_, pos,
                                                                   &Cel);
        ASSERT_LE(std::abs(k_ti - k_i), 1e-14)
            << "Calculated bulk modulus by the transverse isotropic elastic "
               "model: "
            << k_ti << "\n"
            << "Calculated bulk modulus by the isotropic elastic model:" << k_i;
    }

    template <int Dim>
    void checkRotationElasticOrthotropic(
        std::optional<ParameterLib::CoordinateSystem> const& coordinate_system)
    {
        // read in moduli
        ParameterLib::SpatialPosition const pos;
        auto const parameters_orth =
            getParametersForLinearElasticOrthotropic(Dim);
        double const E1 = static_cast<ParameterLib::Parameter<double> const&>(
            *parameters_orth[0])(t_, pos)[0];
        double const E2 = static_cast<ParameterLib::Parameter<double> const&>(
            *parameters_orth[0])(t_, pos)[1];
        double const E3 = static_cast<ParameterLib::Parameter<double> const&>(
            *parameters_orth[0])(t_, pos)[2];
        double const nu12 = static_cast<ParameterLib::Parameter<double> const&>(
            *parameters_orth[1])(t_, pos)[0];
        double const nu23 = static_cast<ParameterLib::Parameter<double> const&>(
            *parameters_orth[1])(t_, pos)[1];
        double const nu13 = static_cast<ParameterLib::Parameter<double> const&>(
            *parameters_orth[1])(t_, pos)[2];
        double const G12 = static_cast<ParameterLib::Parameter<double> const&>(
            *parameters_orth[2])(t_, pos)[0];
        double const G23 = static_cast<ParameterLib::Parameter<double> const&>(
            *parameters_orth[2])(t_, pos)[1];
        double const G13 = static_cast<ParameterLib::Parameter<double> const&>(
            *parameters_orth[2])(t_, pos)[2];

        double const nu21 = nu12 * E2 / E1;
        double const nu32 = nu23 * E3 / E2;
        double const nu31 = nu13 * E3 / E1;
        double const D = (1 - nu12 * nu21 - nu23 * nu32 - nu31 * nu13 -
                          2 * nu12 * nu23 * nu31) /
                         (E1 * E2 * E3);

        // set rotation matrix for 2D and 3D
        Eigen::Matrix3d R_c = Eigen::Matrix3d::Zero();
        if (Dim == 3)
        {
            R_c = coordinate_system->transformation<3>(pos);
        }
        else
        {
            R_c.topLeftCorner<2, 2>() =
                coordinate_system->transformation<2>(pos);
            R_c(2, 2) = 1.0;
        }

        std::optional<ParameterLib::CoordinateSystem> const
            local_coordinate_system{};

        // set orthotropic model in material coordinates
        auto const elastic_model_orthotropic_local_material =
            Tests::createTestConstitutiveRelation<
                MaterialLib::Solids::LinearElasticOrthotropic<Dim>>(
                xml_orth, parameters_orth, local_coordinate_system, false,
                MaterialLib::Solids::createLinearElasticOrthotropic<Dim>);

        auto const C_orthotropic_local_material =
            elastic_model_orthotropic_local_material->getElasticTensor(t_, pos,
                                                                       T_ref_);

        // set orthotropic model in global coordinates
        auto const elastic_model_orthotropic_global =
            Tests::createTestConstitutiveRelation<
                MaterialLib::Solids::LinearElasticOrthotropic<Dim>>(
                xml_orth, parameters_orth, coordinate_system, false,
                MaterialLib::Solids::createLinearElasticOrthotropic<Dim>);
        auto const C_orthotropic_global =
            elastic_model_orthotropic_global->getElasticTensor(t_, pos, T_ref_);

        Tensor4R tensor_local_material = {};

        // set tensor elements
        tensor_local_material[0][0][0][0] = (1 - nu23 * nu32) / (E2 * E3 * D);
        tensor_local_material[0][0][1][1] =
            (nu21 + nu31 * nu23) / (E2 * E3 * D);
        tensor_local_material[0][0][2][2] =
            (nu31 + nu21 * nu32) / (E2 * E3 * D);
        tensor_local_material[1][1][0][0] =
            (nu12 + nu13 * nu32) / (E1 * E3 * D);
        tensor_local_material[1][1][1][1] = (1 - nu31 * nu13) / (E1 * E3 * D);
        tensor_local_material[1][1][2][2] =
            (nu32 + nu31 * nu12) / (E1 * E3 * D);
        tensor_local_material[2][2][0][0] =
            (nu13 + nu12 * nu23) / (E1 * E2 * D);
        tensor_local_material[2][2][1][1] =
            (nu23 + nu13 * nu21) / (E1 * E2 * D);
        tensor_local_material[2][2][2][2] = (1 - nu12 * nu21) / (E1 * E2 * D);
        tensor_local_material[1][2][1][2] = G23;
        tensor_local_material[2][1][1][2] = tensor_local_material[1][2][1][2];
        tensor_local_material[1][2][2][1] = tensor_local_material[1][2][1][2];
        tensor_local_material[2][1][2][1] = tensor_local_material[1][2][1][2];
        tensor_local_material[0][2][0][2] = G13;
        tensor_local_material[2][0][0][2] = tensor_local_material[0][2][0][2];
        tensor_local_material[0][2][2][0] = tensor_local_material[0][2][0][2];
        tensor_local_material[2][0][2][0] = tensor_local_material[0][2][0][2];
        tensor_local_material[0][1][0][1] = G12;
        tensor_local_material[1][0][0][1] = tensor_local_material[0][1][0][1];
        tensor_local_material[0][1][1][0] = tensor_local_material[0][1][0][1];
        tensor_local_material[1][0][1][0] = tensor_local_material[0][1][0][1];

        checkTensors<Dim>(C_orthotropic_local_material, tensor_local_material);
        auto const tensor_global = rotateTensor(tensor_local_material, R_c);
        checkTensors<Dim>(C_orthotropic_global, tensor_global);
    }

private:
    double const t_ = std::numeric_limits<double>::quiet_NaN();
    double const T_ref_ = std::numeric_limits<double>::quiet_NaN();

    std::vector<std::unique_ptr<ParameterLib::ParameterBase>>
    getParametersForLinearElasticOrthotropic(int const dim)
    {
        if (dim == 2)
        {
            // Since the model LinearElasticTransverseIsotropic defines the unit
            // direction of anisotropy as the second base (base2) of the 2D
            // local coordinate system, the following data are computed from
            // that for the  model LinearElasticTransverseIsotropic:
            std::vector<double> E0{8.0e9, 4.0e9, 8.0e9};
            std::vector<double> nu0{0.25, 0.125, 0.35};
            std::vector<double> G0{1.2e9, 0.0, 0.0};
            return setParametersForLinearElasticOrthotropic(E0, nu0, G0);
        }

        std::vector<double> E0{8.0e9, 8.0e9, 4.0e9};
        std::vector<double> nu0{0.35, 0.25, 0.25};
        std::vector<double> G0{2.962962962963e+09, 1.2e9, 1.2e9};
        return setParametersForLinearElasticOrthotropic(E0, nu0, G0);
    }

    Tensor4R rotateTensor(Tensor4R c, Eigen::Matrix<double, 3, 3> R)
    {
        Tensor4R c_rot = {};

        auto outer_indices =
            ranges::views::cartesian_product(ranges::views::iota(0, 3),
                                             ranges::views::iota(0, 3),
                                             ranges::views::iota(0, 3),
                                             ranges::views::iota(0, 3));
        auto inner_indices = outer_indices;

        for (const auto& [i, j, k, l] : outer_indices)
        {
            for (const auto& [r, s, t, u] : inner_indices)
            {
                c_rot[i][j][k][l] +=
                    R(i, r) * R(j, s) * R(k, t) * R(l, u) * c[r][s][t][u];
            }
        }
        return c_rot;
    }

    std::tuple<int, int> kelvinVectorIndexTo2ndRankTensorIndices(int i)
    {
        switch (i)
        {
            case 3:
                return {0, 1};
            case 4:
                return {1, 2};
            case 5:
                return {0, 2};
            default:
                return {i, i};
        }
    }

    template <int Dimension>
    void checkTensors(MaterialLib::Solids::LinearElasticOrthotropic<
                          Dimension>::KelvinMatrix C,
                      Tensor4R c)
    {
        for (int i = 0; i < C.cols(); i++)
        {
            for (int j = 0; j < C.rows(); j++)
            {
                double fac = 1.0;
                auto const [k, l] = kelvinVectorIndexTo2ndRankTensorIndices(i);
                auto const [m, n] = kelvinVectorIndexTo2ndRankTensorIndices(j);
                if ((i < 3 && j > 2) || (j < 3 && i > 2))
                {
                    fac = 1 / std::sqrt(2);
                }
                else if ((i > 2) && (j > 2))
                {
                    fac = 0.5;
                }
                ASSERT_NEAR(C(i, j) * fac, c[k][l][m][n], 5e-6);
            }
        }
    }
};

TEST_F(LinearElasticTransverseIsotropic, test_agaist_ElasticOrthotropic)
{
    {  // 2D
        ParameterLib::ConstantParameter<double> const e1{
            "e1", {std::cos(phi), std::sin(phi)}};
        ParameterLib::ConstantParameter<double> const e2{
            "e2", {-std::sin(phi), std::cos(phi)}};
        ParameterLib::CoordinateSystem const coordinate_system{e1, e2};

        compareWithElasticOrthotropic<2>(coordinate_system);
    }
    {  // 3D
        ParameterLib::ConstantParameter<double> const e1{
            "e1", {std::cos(phi), 0.0, -std::sin(phi)}};
        ParameterLib::ConstantParameter<double> const e2{"e2", {0.0, 1.0, 0.0}};
        ParameterLib::ConstantParameter<double> const e3{
            "e3", {std::sin(phi), 0.0, std::cos(phi)}};
        ParameterLib::CoordinateSystem const coordinate_system{e1, e2, e3};

        compareWithElasticOrthotropic<3>(coordinate_system);
    }
}

TEST_F(LinearElasticTransverseIsotropic, test_Rotation_ElasticOrthotropic)
{
    {  // 2D
        ParameterLib::ConstantParameter<double> const e1{
            "e1", {std::cos(phi), std::sin(phi)}};
        ParameterLib::ConstantParameter<double> const e2{
            "e2", {-std::sin(phi), std::cos(phi)}};
        ParameterLib::CoordinateSystem const coordinate_system{e1, e2};

        checkRotationElasticOrthotropic<2>(coordinate_system);
    }
    {  // 3D
       // x: 15.3° y: 35° z: 44.9°
        ParameterLib::ConstantParameter<double> const e1{
            "e1", {0.5802380261233806, 0.5782161401269231, -0.573576436351046}};
        ParameterLib::ConstantParameter<double> const e2{
            "e2",
            {-0.5736454596106133, 0.7900690700491165, 0.21615214831190652}};
        ParameterLib::ConstantParameter<double> const e3{
            "e3",
            {0.5781476625470101, 0.20360982257358468, 0.7901191811638179}};
        ParameterLib::CoordinateSystem const coordinate_system{e1, e2, e3};

        checkRotationElasticOrthotropic<3>(coordinate_system);
    }
}

TEST_F(LinearElasticTransverseIsotropic,
       test_agaist_ElasticOrthotropicWithImplicitCoordinateBase)
{
    {  // 2D
        ParameterLib::ConstantParameter<double> const line_direction{
            "e3", {-std::sin(phi), std::cos(phi)}};
        ParameterLib::CoordinateSystem const coordinate_system{line_direction};

        compareWithElasticOrthotropic<2>(coordinate_system);
    }
    {  // 3D
        ParameterLib::ConstantParameter<double> const line_direction{
            "e3", {-std::sin(phi), 0.0, std::cos(phi)}};
        ParameterLib::CoordinateSystem const coordinate_system{line_direction};

        compareWithElasticOrthotropic<3>(coordinate_system);
    }
}

TEST_F(LinearElasticTransverseIsotropic, test_agaist_LinearElasticIsotropic)
{
    std::optional<ParameterLib::CoordinateSystem> coordinate_system = {};
    compareWithLinearElasticIsotropic<2>(coordinate_system);
    compareWithLinearElasticIsotropic<3>(coordinate_system);
}
