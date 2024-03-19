/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on November 7, 2023, 10:39 AM
 */

#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <string>

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
};

TEST_F(LinearElasticTransverseIsotropic, test_agaist_ElasticOrthotropic)
{
    {  // 2D
        ParameterLib::ConstantParameter<double> const e1{
            "e1", {0.8191520442889918, 0.573576436351046}};
        ParameterLib::ConstantParameter<double> const e2{
            "e2", {-0.573576436351046, 0.8191520442889918}};
        ParameterLib::CoordinateSystem const coordinate_system{e1, e2};

        compareWithElasticOrthotropic<2>(coordinate_system);
    }
    {  // 3D
        ParameterLib::ConstantParameter<double> const e1{
            "e1", {-0.8191520442889918, 0.0, -0.573576436351046}};
        ParameterLib::ConstantParameter<double> const e2{"e2",
                                                         {0.0, -1.0, 0.0}};
        ParameterLib::ConstantParameter<double> const e3{
            "e3", {-0.573576436351046, 0.0, 0.8191520442889918}};
        ParameterLib::CoordinateSystem const coordinate_system{e1, e2, e3};

        compareWithElasticOrthotropic<3>(coordinate_system);
    }
}

TEST_F(LinearElasticTransverseIsotropic,
       test_agaist_ElasticOrthotropicWithImplicitCoordinateBase)
{
    {  // 2D
        ParameterLib::ConstantParameter<double> const line_direction{
            "e3", {-0.573576436351046, 0.8191520442889918}};
        ParameterLib::CoordinateSystem const coordinate_system{line_direction};

        compareWithElasticOrthotropic<2>(coordinate_system);
    }
    {  // 3D
        ParameterLib::ConstantParameter<double> const line_direction{
            "e3", {-0.573576436351046, 0.0, 0.8191520442889918}};
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
