/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  Created on March 27, 2020, 3:01 PM
 *
 */

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <limits>

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/CapillaryPressureRegularizedVanGenuchten.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/CreateCapillaryPressureRegularizedVanGenuchten.h"
#include "TestMPL.h"
#include "Tests/TestTools.h"

namespace MPL = MaterialPropertyLib;

TEST(MaterialPropertyLib, CapillaryPressureRegularizedVanGenuchten)
{
    const char xml_pc_S[] =
        "<property>"
        "   <name>capillary_pressure</name>"
        "   <type>CapillaryPressureRegularizedVanGenuchten</type>"
        "   <residual_liquid_saturation>0.1</residual_liquid_saturation>"
        "   <residual_gas_saturation>0.05</residual_gas_saturation>"
        "   <exponent>0.6</exponent>"
        "   <p_b>1.e+4</p_b>"
        "</property>";
    auto const capillary_pressure_property_ptr = Tests::createTestProperty(
        xml_pc_S, MPL::createCapillaryPressureRegularizedVanGenuchten);

    MPL::Property const& capillary_pressure_property =
        *capillary_pressure_property_ptr;

    double const Sl_r = 0.1;

    MPL::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    std::array<double, 31> S{
        {0,        0.0333333, 0.0666667, 0.1,      0.133333, 0.166667, 0.2,
         0.233333, 0.266667,  0.3,       0.333333, 0.366667, 0.4,      0.433333,
         0.466667, 0.5,       0.533333,  0.566667, 0.6,      0.633333, 0.666667,
         0.7,      0.733333,  0.766667,  0.8,      0.833333, 0.866667, 0.9,
         0.933333, 0.966667,  1}};
    // Data for S in [0, S_max] calculated by the old implementation.
    // When S in [Smax, 1], pc=0.
    std::array<double, 31> pc{{5.3649187738e+11,
                               3.5767283022e+11,
                               1.7885324659e+11,
                               3.4199425947e+07,
                               8.6378746283e+04,
                               5.4166427767e+04,
                               4.1081248712e+04,
                               3.3651450177e+04,
                               2.8734971602e+04,
                               2.5176752591e+04,
                               2.2443387810e+04,
                               2.0251583348e+04,
                               1.8435732700e+04,
                               1.6891590993e+04,
                               1.5549880275e+04,
                               1.4362497435e+04,
                               1.3294509680e+04,
                               1.2319676236e+04,
                               1.1417620350e+04,
                               1.0571723353e+04,
                               9.7678323835e+03,
                               8.9932086994e+03,
                               8.2353166956e+03,
                               7.4806165205e+03,
                               6.7126607840e+03,
                               5.9081641105e+03,
                               5.0279021651e+03,
                               3.9870638390e+03,
                               2.4795040664e+03,
                               0.0,
                               0.0}};

    const int n = 30;
    for (int i = 0; i <= n; ++i)
    {
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S[i];

        const double calculated_pc =
            capillary_pressure_property.template value<double>(variable_array,
                                                               pos, t, dt);

        if (std::fabs(calculated_pc) < std::numeric_limits<double>::epsilon())
        {
            ASSERT_LE(std::fabs(pc[i] - calculated_pc), 1e-16)
                << "for the reference pc " << pc[i] << " and the calculated pc."
                << calculated_pc;
        }
        else
        {
            ASSERT_LE(std::fabs((pc[i] - calculated_pc) / calculated_pc), 1e-10)
                << "for the reference pc " << pc[i] << " and the calculated pc."
                << calculated_pc;
        }

        // The regularized function in not smooth at Sl_r
        if (S[i] == Sl_r)
        {
            continue;
        }

        if (i == 0 || i == n)
        {
            continue;
        }

        const double dPcdS =
            capillary_pressure_property.template dValue<double>(
                variable_array,
                MaterialPropertyLib::Variable::liquid_saturation, pos, t, dt);

        const double offset = (S[i] < Sl_r) ? 1.0e-2 : 1.0e-7;
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S[i] - offset;

        const double pc0 = capillary_pressure_property.template value<double>(
            variable_array, pos, t, dt);
        variable_array[static_cast<int>(MPL::Variable::liquid_saturation)] =
            S[i] + offset;
        const double pc1 = capillary_pressure_property.template value<double>(
            variable_array, pos, t, dt);

        const double numerical_dPcdS = 0.5 * (pc1 - pc0) / offset;

        if (std::fabs(numerical_dPcdS) < std::numeric_limits<double>::epsilon())
        {
            ASSERT_LE(std::fabs((dPcdS - numerical_dPcdS)), 1e-16)
                << "for the analytic derivative of dPc/dS " << dPcdS
                << " and numeric derivative of dPc/dS." << numerical_dPcdS;
        }
        else
        {
            ASSERT_LE(std::fabs((dPcdS - numerical_dPcdS) / dPcdS), 1e-8)
                << "for the analytic derivative of dPc/dS " << dPcdS
                << " and numeric derivative of dPc/dS." << numerical_dPcdS;
        }
    }
}
