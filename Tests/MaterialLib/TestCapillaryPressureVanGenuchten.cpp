/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  Created on March 23, 2020, 9:29 AM
 */

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <random>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/CapillaryPressureVanGenuchten.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/SaturationVanGenuchten.h"
#include "TestMPL.h"
#include "Tests/TestTools.h"

namespace MPL = MaterialPropertyLib;

TEST(MaterialPropertyLib, CapillaryPressureVanGenuchten)
{
    double const residual_liquid_saturation = 0.1;
    double const residual_gas_saturation = 0.05;
    double const exponent = 0.79;
    double const entry_pressure = 5000;

    MPL::Property const& saturation_property = MPL::SaturationVanGenuchten{
        residual_liquid_saturation, residual_gas_saturation, exponent,
        entry_pressure};

    double const max_capillary_pressure = std::numeric_limits<double>::max();
    MPL::Property const& capillary_pressure_property =
        MPL::CapillaryPressureVanGenuchten{
            residual_liquid_saturation, residual_gas_saturation, exponent,
            entry_pressure, max_capillary_pressure};

    MPL::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    std::random_device rd;
    std::mt19937 mt(rd());
    const double offset = std::sqrt(std::numeric_limits<double>::epsilon());
    std::uniform_real_distribution<double> distributor(
        residual_liquid_saturation + offset,
        1.0 - residual_gas_saturation - offset);

    const int n = 20;
    for (int i = 0; i <= n; ++i)
    {
        double const S = distributor(mt);
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] = S;

        const double pc = capillary_pressure_property.template value<double>(
            variable_array, pos, t, dt);
        variable_array[static_cast<int>(MPL::Variable::capillary_pressure)] =
            pc;

        const double computed_saturation =
            saturation_property.template value<double>(variable_array, pos, t,
                                                       dt);

        ASSERT_LE(std::fabs(S - computed_saturation), 1e-9)
            << "for saturation " << S
            << " and re-computed saturation via capillary pressure"
            << computed_saturation;

        const double dPcdS =
            capillary_pressure_property.template dValue<double>(
                variable_array,
                MaterialPropertyLib::Variable::liquid_saturation, pos, t, dt);

        const double d2PcdS2 =
            capillary_pressure_property.template d2Value<double>(
                variable_array,
                MaterialPropertyLib::Variable::liquid_saturation,
                MaterialPropertyLib::Variable::liquid_saturation, pos, t, dt);

        const double S1 =
            std::clamp(S + offset, residual_liquid_saturation + 2 * offset,
                       1.0 - residual_gas_saturation - 2 * offset);
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] = S1;

        if (std::fabs(S1 - S) < std::numeric_limits<double>::epsilon())
        {
            continue;
        }

        const double pc1 = capillary_pressure_property.template value<double>(
            variable_array, pos, t, dt);
        const double numerical_dPcdS = (pc1 - pc) / (S1 - S);
        ASSERT_LE(std::fabs((dPcdS - numerical_dPcdS) / dPcdS), 1e-4)
            << "for the analytic derivative of dPc/dS " << dPcdS
            << " and numeric derivative of dPc/dS." << numerical_dPcdS;

        const double dPcdS1 =
            capillary_pressure_property.template dValue<double>(
                variable_array,
                MaterialPropertyLib::Variable::liquid_saturation, pos, t, dt);

        const double numerical_d2PcdS2 = (dPcdS1 - dPcdS) / (S1 - S);
        ASSERT_LE(std::fabs((d2PcdS2 - numerical_d2PcdS2) / d2PcdS2), 1e-4)
            << "for the analytic second order derivative of d2Pc/dS2 "
            << d2PcdS2 << " and numeric  second order  derivative of d2Pc/dS2."
            << numerical_d2PcdS2;
    }
}

TEST(MaterialPropertyLib, RegularizedCapillaryVanGenuchten)
{
    double const residual_liquid_saturation = 0.1;
    double const residual_gas_saturation = 0.05;
    double const exponent = 0.79;
    double const entry_pressure = 5000;

    MPL::Property const& saturation_property = MPL::SaturationVanGenuchten{
        residual_liquid_saturation, residual_gas_saturation, exponent,
        entry_pressure};

    double const max_capillary_pressure = std::numeric_limits<double>::max();
    MPL::Property const& regularized_capillary_pressure_property =
        MPL::CapillaryPressureRegularizedVanGenuchten{
            residual_liquid_saturation, residual_gas_saturation, exponent,
            entry_pressure, max_capillary_pressure};

    MPL::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    std::random_device rd;
    std::mt19937 mt(rd());
    const double offset = std::sqrt(std::numeric_limits<double>::epsilon());
    std::uniform_real_distribution<double> distributor(
        residual_liquid_saturation + offset,
        1.0 - residual_gas_saturation - offset);

    const int n = 20;
    for (int i = 0; i <= n; ++i)
    {
        double const S = distributor(mt);
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] = S;

        const double pc =
            regularized_capillary_pressure_property.template value<double>(
                variable_array, pos, t, dt);
        variable_array[static_cast<int>(MPL::Variable::capillary_pressure)] =
            pc;

        const double dPcdS =
            regularized_capillary_pressure_property.template dValue<double>(
                variable_array,
                MaterialPropertyLib::Variable::liquid_saturation, pos, t, dt);

        const double d2PcdS2 =
            regularized_capillary_pressure_property.template d2Value<double>(
                variable_array,
                MaterialPropertyLib::Variable::liquid_saturation,
                MaterialPropertyLib::Variable::liquid_saturation, pos, t, dt);

        const double S1 =
            std::clamp(S + offset, residual_liquid_saturation + 2 * offset,
                       1.0 - residual_gas_saturation - 2 * offset);
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] = S1;

        if (std::fabs(S1 - S) < std::numeric_limits<double>::epsilon())
        {
            continue;
        }

        const double pc1 =
            regularized_capillary_pressure_property.template value<double>(
                variable_array, pos, t, dt);
        const double numerical_dPcdS = (pc1 - pc) / (S1 - S);
        ASSERT_LE(std::fabs((dPcdS - numerical_dPcdS) / dPcdS), 1e-4)
            << "for the analytic derivative of dPc/dS " << dPcdS
            << " and numeric derivative of dPc/dS." << numerical_dPcdS;

        const double dPcdS1 =
            regularized_capillary_pressure_property.template dValue<double>(
                variable_array,
                MaterialPropertyLib::Variable::liquid_saturation, pos, t, dt);

        const double numerical_d2PcdS2 = (dPcdS1 - dPcdS) / (S1 - S);
        ASSERT_LE(std::fabs((d2PcdS2 - numerical_d2PcdS2) / d2PcdS2), 1e-4)
            << "for the analytic second order derivative of d2Pc/dS2 "
            << d2PcdS2 << " and numeric  second order  derivative of d2Pc/dS2."
            << numerical_d2PcdS2;
    }
}
