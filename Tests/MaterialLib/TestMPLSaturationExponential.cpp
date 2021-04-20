/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include <gtest/gtest.h>
#include <iomanip>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/CapillaryPressureSaturation/SaturationExponential.h"
#include "TestMPL.h"
#include "Tests/TestTools.h"

namespace MPL = MaterialPropertyLib;

TEST(MaterialPropertyLib, SaturationExponential)
{
    double const residual_liquid_saturation = 0.1;
    double const residual_gas_saturation = 0.05;
    double const exponent = 1.5;
    double const p_cap_ref = 30000;

    MPL::Property const& saturation = MPL::SaturationExponential{
        "saturation", residual_liquid_saturation, residual_gas_saturation,
        p_cap_ref, exponent};

    MPL::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    double const p_0 = -p_cap_ref;
    double const p_max = 1.5 * p_cap_ref;
    int const n_steps = 9999;
    for (int i = 0; i <= n_steps; ++i)
    {
        double const p_cap = p_0 + i * (p_max - p_0) / n_steps;
        variable_array[static_cast<int>(MPL::Variable::capillary_pressure)] =
            p_cap;

        double const s_res = residual_liquid_saturation;
        double const s_max = 1. - residual_gas_saturation;
        double const s_eff =
            1. -
            std::pow(std::clamp(p_cap, 0., p_cap_ref) / p_cap_ref, exponent);

        double s_ref = s_eff * (s_max - s_res) + s_res;

        double const S =
            saturation.template value<double>(variable_array, pos, t, dt);
        double const dS = saturation.template dValue<double>(
            variable_array, MPL::Variable::capillary_pressure, pos, t, dt);

        double const eps = 1e-2;
        variable_array[static_cast<int>(MPL::Variable::capillary_pressure)] =
            p_cap - eps;

        double const S_minus =
            saturation.template value<double>(variable_array, pos, t, dt);

        variable_array[static_cast<int>(MPL::Variable::capillary_pressure)] =
            p_cap + eps;

        double const S_plus =
            saturation.template value<double>(variable_array, pos, t, dt);

        double const dS_ref = p_cap > 0 ? (S_plus - S_minus) / 2 / eps : 0.;

        ASSERT_LE(std::abs(s_ref - S), 1e-9)
            << std::setprecision(16) << "with: \ncapillary pressure: " << p_cap
            << "\nsaturation: " << S << "\nsaturation_ref = " << s_ref
            << "\nat point: " << i;
        ASSERT_LE(std::abs(dS - dS_ref), 1e-9)
            << std::setprecision(16) << "with: \np_cap: " << p_cap
            << "\nS: " << S << "\nS_ref = " << s_ref << "\ndS = " << dS
            << "\ndS_ref = " << dS_ref << "\nat point: " << i;
    }
}
