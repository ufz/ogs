/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on January 12, 2021, x:xx AM
 */

#include <gtest/gtest.h>

#include <Eigen/Eigen>
#include <boost/math/constants/constants.hpp>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/CreateGasPressureDependentPermeability.h"
#include "MaterialLib/MPL/Properties/GasPressureDependentPermeability.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "ParameterLib/ConstantParameter.h"
#include "TestMPL.h"
#include "Tests/TestTools.h"

TEST(MaterialPropertyLib, GasPressureDependentPermeability)
{
    ParameterLib::ConstantParameter<double> const k0("k0", 1.e-20);
    double const a1 = 0.125;
    double const a2 = 152.0;
    double const pressure_threshold = 3.2;
    double const min_permeability = 1.e-22;
    double const max_permeability = 1.e-10;

    auto const k_model = MPL::GasPressureDependentPermeability<3>(
        "k_gas", k0, a1, a2, pressure_threshold, min_permeability,
        max_permeability, nullptr);

    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    MPL::VariableArray vars;

    /// For gas pressure smaller than threshold value.
    {
        double const p_gas = 2.5;

        vars[static_cast<int>(MPL::Variable::phase_pressure)] = p_gas;
        auto const k = MPL::formEigenTensor<3>(k_model.value(vars, pos, t, dt));

        double const k_expected = 1.312500000000000000e-20;

        ASSERT_LE(std::fabs(k_expected - k(0, 0)) / k_expected, 1e-10)
            << "for expected permeability with gas pressure below threshold"
            << k_expected
            << " and for computed permeability with gas pressure below "
               "threshold "
            << k(0, 0);
    }
    /// For gas pressure bigger than threshold value.
    {
        double const p_gas = 4.5;

        vars[static_cast<int>(MPL::Variable::phase_pressure)] = p_gas;
        auto const k = MPL::formEigenTensor<3>(k_model.value(vars, pos, t, dt));

        double const k_expected = 1.990000000000000000000e-18;

        ASSERT_LE(std::fabs(k_expected - k(0, 0)) / k_expected, 1e-10)
            << "for expected permeability with gas pressure above threshold "
            << k_expected
            << " and for computed permeability with gas pressure above "
               "threshold "
            << k(0, 0);
    }
}
