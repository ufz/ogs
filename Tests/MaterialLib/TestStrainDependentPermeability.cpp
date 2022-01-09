/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on November 10, 2020, 10:34 AM
 */

#include <gtest/gtest.h>

#include <Eigen/Eigen>
#include <boost/math/constants/constants.hpp>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/CreateStrainDependentPermeability.h"
#include "MaterialLib/MPL/Properties/StrainDependentPermeability.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "ParameterLib/ConstantParameter.h"
#include "TestMPL.h"
#include "Tests/TestTools.h"

TEST(MaterialPropertyLib, StrainDependentPermeability)
{
    ParameterLib::ConstantParameter<double> const k0("k0", 1.e-20);
    double const b1 = 4000.0;
    double const b2 = 300.0;
    double const b3 = 50.0;
    double const min_permeability = 1.e-22;
    double const max_permeability = 1.e-10;

    auto const k_model = MPL::StrainDependentPermeability<3>(
        "k_strain", k0, b1, b2, b3, min_permeability, max_permeability,
        nullptr);

    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    MPL::VariableArray vars;

    /// Extensive volumetric strain, which increases permeability.
    {
        double const e_vol = 1.e-4;
        double const e_vol_pls = 1.e-5;

        vars[static_cast<int>(MPL::Variable::volumetric_strain)] = e_vol;
        vars[static_cast<int>(MPL::Variable::equivalent_plastic_strain)] =
            e_vol_pls;
        auto const k = MPL::formEigenTensor<3>(k_model.value(vars, pos, t, dt));

        double const k_expected = 1.0528627949453586e-20;

        ASSERT_LE(std::fabs(k_expected - k(0, 0)) / k_expected, 1e-10)
            << "for expected increased permeability " << k_expected
            << " and for computed increased permeability " << k(0, 0);
    }
    /// Compact volumetric strain, which decreases permeability.
    {
        double const e_vol = -1.e-4;
        double const e_vol_pls = 1.e-5;

        vars[static_cast<int>(MPL::Variable::volumetric_strain)] = e_vol;
        vars[static_cast<int>(MPL::Variable::equivalent_plastic_strain)] =
            e_vol_pls;
        auto const k = MPL::formEigenTensor<3>(k_model.value(vars, pos, t, dt));

        double const k_expected = 9.7134113133089206e-21;

        ASSERT_LE(std::fabs(k_expected - k(0, 0)) / k_expected, 1e-10)
            << "for expected decreased permeability " << k_expected
            << " and for computed decreased permeability " << k(0, 0);
    }
}
