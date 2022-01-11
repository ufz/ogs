/**
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
 *
 * Created on November 2, 2016, 3:09 PM
 */

#include <gtest/gtest.h>

#include <vector>

#include "MaterialLib/MPL/Properties/CreateVermaPruessModel.h"
#include "MaterialLib/MPL/Properties/VermaPruessModel.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "ParameterLib/ConstantParameter.h"
#include "TestMPL.h"
#include "Tests/TestTools.h"

TEST(MaterialPropertyLib, PermeabilityVermaPruessModel)
{
    ParameterLib::ConstantParameter<double> const k0("permeability_ini",
                                                     1.5e-14);
    ParameterLib::ConstantParameter<double> const phi0("porosity_ini", 0.46);
    ParameterLib::ConstantParameter<double> const phi_c("critical_porosity",
                                                        0.3956);
    ParameterLib::ConstantParameter<double> const n("exponent", 5);

    auto const k_model = MPL::VermaPruessModel(k0, phi0, phi_c, n);

    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    std::vector<double> const porosity{0.44, 0.38};
    std::vector<double> const k_expected{2.3365591132239605e-15, 0.};

    MPL::VariableArray vars;

    for (std::size_t i = 0; i < porosity.size(); ++i)
    {
        vars[static_cast<int>(MaterialPropertyLib::Variable::porosity)] =
            porosity[i];

        auto const k = MPL::formEigenTensor<1>(k_model.value(vars, pos, t, dt));

        // Compare the calculated value with the expected value.
        if (k_expected[i] == 0.)
        {
            ASSERT_EQ(k_expected[i], k(0, 0));
        }
        else
        {
            ASSERT_LE(std::fabs(k_expected[i] - k(0, 0)) / k_expected[i],
                      1e-10);
        }
    }
}
