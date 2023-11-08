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

#include <memory>
#include <string>

#include "CreateTestConstitutiveRelation.h"
#include "MaterialLib/SolidModels/CreateLinearElasticOrthotropic.h"
#include "MaterialLib/SolidModels/CreateLinearElasticTransverseIsotropic.h"
#include "MaterialLib/SolidModels/LinearElasticOrthotropic.h"
#include "MaterialLib/SolidModels/LinearElasticTransverseIsotropic.h"
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

TEST(MaterialLib, LinearElasticTransverseIsotropic)
{
    ParameterLib::ConstantParameter<double> const e1{
        "e1", {-0.8191520442889918, 0.0, -0.573576436351046}};
    ParameterLib::ConstantParameter<double> const e2{"e2", {0.0, -1.0, 0.0}};
    ParameterLib::ConstantParameter<double> const e3{
        "e3", {-0.573576436351046, 0.0, 0.8191520442889918}};
    /*ParameterLib::ConstantParameter<double> const e1{"e1", {1, 0.0, 0.0}};
    ParameterLib::ConstantParameter<double> const e2{"e2", {0.0, 1.0, 0.0}};
    ParameterLib::ConstantParameter<double> const e3{"e3", {0, 0.0, 1.0}};*/
    ParameterLib::CoordinateSystem const coordinate_system{e1, e2, e3};

    double const t = std::numeric_limits<double>::quiet_NaN();
    ParameterLib::SpatialPosition const pos;
    double const T_ref = std::numeric_limits<double>::quiet_NaN();

    auto const parameters_ti = setParametersForLinearElasticTransverseIsotropic(
        8.0e9, 4.0e9, 0.35, 0.25, 1.2e9);

    auto const elastic_model_transverse_isotropy =
        Tests::createTestConstitutiveRelation<
            MaterialLib::Solids::LinearElasticTransverseIsotropic<3>>(
            xml_ti, parameters_ti, coordinate_system, false,
            MaterialLib::Solids::createLinearElasticTransverseIsotropic<3>);
    auto const E =
        elastic_model_transverse_isotropy->getElasticTensor(t, pos, T_ref);

    std::vector<double> E0{8.0e9, 8.0e9, 4.0e9};
    std::vector<double> nu0{0.35, 0.25, 0.25};
    std::vector<double> G0{2.962962962963e+09, 1.2e9, 1.2e9};
    auto const parameters_orth =
        setParametersForLinearElasticOrthotropic(E0, nu0, G0);

    auto const elastic_model_orthotropic =
        Tests::createTestConstitutiveRelation<
            MaterialLib::Solids::LinearElasticOrthotropic<3>>(
            xml_orth, parameters_orth, coordinate_system, false,
            MaterialLib::Solids::createLinearElasticOrthotropic<3>);

    auto const E_oth =
        elastic_model_orthotropic->getElasticTensor(t, pos, T_ref);

    ASSERT_LE((E - E_oth).norm() / E.norm(), 1e-14);

    MathLib::KelvinVector::KelvinMatrixType<3> Cel;

    double const k_ti =
        elastic_model_transverse_isotropy->getBulkModulus(t, pos, &Cel);
    ASSERT_LE(k_ti - 4301075268.8172045, 1e-14)
        << "Calculated bulk modulus by the transverse isotropy model: " << k_ti
        << "\n"
        << "Expected Bulk modulus: 4301075268.8172045";

    // The definitions of the bulk modulus by LinearElasticOrthotropic and by
    // LinearElasticTransverseIsotropic are different. Therefore, the comparison
    // of the bulk modulus values by the models are not presented.
}
