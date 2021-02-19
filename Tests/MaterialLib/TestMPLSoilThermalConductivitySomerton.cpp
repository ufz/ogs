/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on February 17, 2021, 4:22 PM
 */

#include <gtest/gtest.h>

#include <Eigen/Eigen>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <vector>

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/ThermalConductivity/CreateSoilThermalConductivitySomerton.h"
#include "MaterialLib/MPL/Properties/ThermalConductivity/SoilThermalConductivitySomerton.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "ParameterLib/ConstantParameter.h"
#include "ParameterLib/CoordinateSystem.h"
#include "Tests/TestTools.h"

namespace MPL = MaterialPropertyLib;

std::unique_ptr<MaterialPropertyLib::Property>
createTestSoilThermalConductivitySomertonProperty(
    const char xml[], int const geometry_dimension,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    ParameterLib::CoordinateSystem const* const local_coordinate_system,
    std::function<std::unique_ptr<MaterialPropertyLib::Property>(
        int const geometry_dimension,
        BaseLib::ConfigTree const& config,
        std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
            parameters,
        ParameterLib::CoordinateSystem const* const local_coordinate_system)>
        createProperty)
{
    auto const ptree = Tests::readXml(xml);
    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& sub_config = conf.getConfigSubtree("property");
    // Parsing the property name:
    auto const property_name =
        sub_config.getConfigParameter<std::string>("name");

    return createProperty(geometry_dimension, sub_config, parameters,
                          local_coordinate_system);
}

TEST(MaterialPropertyLib, SoilThermalConductivitySomerton1D)
{
    const char xml[] =
        "<property>"
        "   <name>thermal_conductivity</name>"
        "   <type>SoilThermalConductivitySomerton</type>"
        "   <dry_thermal_conductivity>k_T_dry</dry_thermal_conductivity>"
        "   <wet_thermal_conductivity>k_T_wet</wet_thermal_conductivity>"
        "</property>";

    double const k_T_dry = 0.1;
    double const k_T_wet = 0.3;
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;
    parameters.push_back(
        std::make_unique<ParameterLib::ConstantParameter<double>>("k_T_dry",
                                                                  k_T_dry));
    parameters.push_back(
        std::make_unique<ParameterLib::ConstantParameter<double>>("k_T_wet",
                                                                  k_T_wet));

    const int dimemsion = 1;
    std::unique_ptr<MPL::Property> const k_T_ptr =
        createTestSoilThermalConductivitySomertonProperty(
            xml, dimemsion, parameters, nullptr,
            MPL::createSoilThermalConductivitySomerton);

    MPL::Property const& k_T_property = *k_T_ptr;

    MPL::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    // For S <= 0.0
    {
        double const S_less_than_zero[] = {-1, 0.0};
        for (int i = 0; i < 2; i++)
        {
            variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
                S_less_than_zero[i];

            double const computed_k_T =
                k_T_property.template value<double>(variable_array, pos, t, dt);
            ASSERT_LE(std::fabs(k_T_dry - computed_k_T), 1e-10)
                << "for expected thermal conductivity " << k_T_dry
                << " and for computed thermal conductivity." << computed_k_T;

            double const computed_dk_T_dS =
                k_T_property.template dValue<double>(
                    variable_array,
                    MaterialPropertyLib::Variable::liquid_saturation, pos, t,
                    dt);

            ASSERT_LE(std::fabs(0.0 - computed_dk_T_dS), 1e-10)
                << "for expected derivative of thermal conductivity with "
                   "respect to saturation "
                << 0.0
                << " and for computed derivative of thermal conductivity "
                   "with respect to saturation."
                << computed_dk_T_dS;
        }
    }
    // For S > 1.0
    {
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
            2.0;

        double const computed_k_T =
            k_T_property.template value<double>(variable_array, pos, t, dt);
        ASSERT_LE(std::fabs(k_T_wet - computed_k_T), 1e-10)
            << "for expected thermal conductivity " << k_T_wet
            << " and for computed thermal conductivity." << computed_k_T;

        double const computed_dk_T_dS = k_T_property.template dValue<double>(
            variable_array, MaterialPropertyLib::Variable::liquid_saturation,
            pos, t, dt);

        ASSERT_LE(std::fabs(0.0 - computed_dk_T_dS), 1e-10)
            << "for expected derivative of thermal conductivity with "
               "respect to saturation "
            << 0.0 << " and for computed derivative of thermal conductivity "
                      "with respect to saturation."
            << computed_dk_T_dS;
    }

    // For S in (0, 1)
    std::vector<double> const S_L = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
    std::vector<double> const expected_k_T = {
        1.632456e-01, 1.894427e-01, 2.095445e-01, 2.264911e-01,
        2.414214e-01, 2.549193e-01, 2.673320e-01, 2.788854e-01};

    double const perturbation = 1.e-6;
    for (std::size_t i = 0; i < S_L.size(); i++)
    {
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S_L[i];

        double const k_T_i =
            k_T_property.template value<double>(variable_array, pos, t, dt);
        ASSERT_LE(std::fabs(expected_k_T[i] - k_T_i), 1e-7)
            << "for expected thermal conductivity " << expected_k_T[i]
            << " and for computed thermal conductivity." << k_T_i;

        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S_L[i] - perturbation;
        double const k_T_i_0 =
            k_T_property.template value<double>(variable_array, pos, t, dt);

        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S_L[i] + perturbation;

        double const k_T_i_1 =
            k_T_property.template value<double>(variable_array, pos, t, dt);

        double const approximated_dk_T_dS =
            0.5 * (k_T_i_1 - k_T_i_0) / perturbation;
        double const analytic_dk_T_dS = k_T_property.template dValue<double>(
            variable_array, MaterialPropertyLib::Variable::liquid_saturation,
            pos, t, dt);

        ASSERT_LE(std::fabs(analytic_dk_T_dS - approximated_dk_T_dS), 1e-5)
            << "for expected derivative of thermal conductivity with "
               "respect to saturation "
            << analytic_dk_T_dS
            << " and for computed derivative of thermal conductivity "
               "with respect to saturation."
            << approximated_dk_T_dS;
    }
}

TEST(MaterialPropertyLib, SoilThermalConductivitySomerton3D)
{
    const char xml[] =
        "<property>"
        "   <name>thermal_conductivity</name>"
        "   <type>SoilThermalConductivitySomerton</type>"
        "   <dry_thermal_conductivity>k_T_dry</dry_thermal_conductivity>"
        "   <wet_thermal_conductivity>k_T_wet</wet_thermal_conductivity>"
        "</property>";

    std::vector<double> k_T_dry{0.1, 0.1, 0.1};
    std::vector<double> k_T_wet{0.3, 0.3, 0.3};
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> parameters;
    parameters.push_back(
        std::make_unique<ParameterLib::ConstantParameter<double>>("k_T_dry",
                                                                  k_T_dry));
    parameters.push_back(
        std::make_unique<ParameterLib::ConstantParameter<double>>("k_T_wet",
                                                                  k_T_wet));

    const int dimemsion = 3;
    std::unique_ptr<MPL::Property> const k_T_ptr =
        createTestSoilThermalConductivitySomertonProperty(
            xml, dimemsion, parameters, nullptr,
            MPL::createSoilThermalConductivitySomerton);

    MPL::Property const& k_T_property = *k_T_ptr;

    MPL::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    // For S <= 0.0
    {
        double const S_less_than_zero[] = {-1, 0.0};
        for (int i = 0; i < 2; i++)
        {
            variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
                S_less_than_zero[i];

            auto const computed_k_T = MPL::formEigenTensor<3>(
                k_T_property.value(variable_array, pos, t, dt));

            for (std::size_t k = 0; k < k_T_dry.size(); k++)
            {
                ASSERT_LE(std::fabs(k_T_dry[k] - computed_k_T(k, k)), 1e-10)
                    << "for expected thermal conductivity " << k_T_dry[k]
                    << " and for computed thermal conductivity."
                    << computed_k_T(k, k);
            }

            auto const computed_dk_T_dS =
                MPL::formEigenTensor<3>(k_T_property.dValue(
                    variable_array,
                    MaterialPropertyLib::Variable::liquid_saturation, pos, t,
                    dt));

            for (std::size_t k = 0; k < k_T_dry.size(); k++)
            {
                ASSERT_LE(std::fabs(0.0 - computed_dk_T_dS(k, k)), 1e-10)
                    << "for expected derivative of thermal conductivity with "
                       "respect to saturation "
                    << 0.0
                    << " and for computed derivative of thermal conductivity "
                       "with respect to saturation."
                    << computed_dk_T_dS(k, k);
            }
        }
    }
    // For S > 1.0
    {
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
            2.0;

        auto const computed_k_T = MPL::formEigenTensor<3>(
            k_T_property.value(variable_array, pos, t, dt));

        for (std::size_t k = 0; k < k_T_dry.size(); k++)
        {
            ASSERT_LE(std::fabs(k_T_wet[k] - computed_k_T(k, k)), 1e-10)
                << "for expected thermal conductivity " << k_T_wet[k]
                << " and for computed thermal conductivity."
                << computed_k_T(k, k);
        }

        auto const computed_dk_T_dS =
            MPL::formEigenTensor<3>(k_T_property.dValue(
                variable_array,
                MaterialPropertyLib::Variable::liquid_saturation, pos, t, dt));

        for (std::size_t k = 0; k < k_T_dry.size(); k++)
        {
            ASSERT_LE(std::fabs(0.0 - computed_dk_T_dS(k, k)), 1e-10)
                << "for expected derivative of thermal conductivity with "
                   "respect to saturation "
                << 0.0
                << " and for computed derivative of thermal conductivity "
                   "with respect to saturation."
                << computed_dk_T_dS(k, k);
        }
    }

    // For S in (0, 1)
    std::vector<double> const S_L = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};
    std::vector<double> const expected_k_T = {
        1.632456e-01, 1.894427e-01, 2.095445e-01, 2.264911e-01,
        2.414214e-01, 2.549193e-01, 2.673320e-01, 2.788854e-01};

    double const perturbation = 1.e-6;
    for (std::size_t i = 0; i < S_L.size(); i++)
    {
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S_L[i];

        auto const k_T_i = MPL::formEigenTensor<3>(
            k_T_property.value(variable_array, pos, t, dt));

        for (std::size_t k = 0; k < k_T_dry.size(); k++)
        {
            ASSERT_LE(std::fabs(expected_k_T[i] - k_T_i(k, k)), 1e-7)
                << "for expected thermal conductivity " << expected_k_T[i]
                << " and for computed thermal conductivity." << k_T_i(k, k);
        }

        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S_L[i] - perturbation;
        auto const k_T_i_0 = MPL::formEigenTensor<3>(
            k_T_property.value(variable_array, pos, t, dt));

        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S_L[i] + perturbation;

        auto const k_T_i_1 = MPL::formEigenTensor<3>(
            k_T_property.value(variable_array, pos, t, dt));

        auto const analytic_dk_T_dS =
            MPL::formEigenTensor<3>(k_T_property.dValue(
                variable_array,
                MaterialPropertyLib::Variable::liquid_saturation, pos, t, dt));

        for (std::size_t k = 0; k < k_T_dry.size(); k++)
        {
            auto const approximated_dk_T_dS =
                0.5 * (k_T_i_1(k, k) - k_T_i_0(k, k)) / perturbation;
            ASSERT_LE(std::fabs(analytic_dk_T_dS(k, k) - approximated_dk_T_dS),
                      1e-5)
                << "for expected derivative of thermal conductivity with "
                   "respect to saturation "
                << analytic_dk_T_dS(k, k)
                << " and for computed derivative of thermal conductivity "
                   "with respect to saturation."
                << approximated_dk_T_dS;
        }
    }
}
