/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on June 8, 2020, 9:17 AM
 */

#include <gtest/gtest.h>

#include <Eigen/Eigen>
#include <boost/math/constants/constants.hpp>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Properties/CreatePermeabilityMohrCoulombFailureIndexModel.h"
#include "MaterialLib/MPL/Properties/PermeabilityMohrCoulombFailureIndexModel.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"
#include "MathLib/KelvinVector.h"
#include "ParameterLib/ConstantParameter.h"
#include "TestMPL.h"
#include "Tests/TestTools.h"

TEST(MaterialPropertyLib, PermeabilityMohrCoulombFailureIndexModel)
{
    ParameterLib::ConstantParameter<double> const k0("k0", 1.e-20);
    double const kr = 1.0e-19;
    double const b = 3.0;
    double const c = 1.0e+6;
    double const phi = 15.0;
    double const k_max = 1.e-10;
    double const t_sigma_max = c;

    auto const k_model = MPL::PermeabilityMohrCoulombFailureIndexModel<3>(
        "k_f", k0, kr, b, c, phi, k_max, t_sigma_max, nullptr);

    const int symmetric_tensor_size = 6;
    using SymmetricTensor = Eigen::Matrix<double, symmetric_tensor_size, 1>;

    // Under failure, i,e stress beyond the yield.
    SymmetricTensor stress;
    stress[0] = -1.36040e+7;
    stress[1] = -1.78344e+7;
    stress[2] = -1.36627e+7;
    stress[3] = -105408;
    stress[4] = -25377.2;
    stress[5] = -1.39944e+7;

    MPL::VariableArray vars;
    vars[static_cast<int>(MPL::Variable::stress)].emplace<SymmetricTensor>(
        stress);

    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    auto const k = MPL::formEigenTensor<3>(k_model.value(vars, pos, t, dt));

    double const k_expected = 1.1398264890628033e-15;

    ASSERT_LE(std::fabs(k_expected - k(0, 0)) / k_expected, 1e-10)
        << "for expected changed permeability " << k_expected
        << " and for computed changed permeability." << k(0, 0);

    // Stress outside of the yield surface. No change in permeability
    stress[0] = -1.2e+7;
    stress[1] = -1.5e+7;
    stress[2] = -1.2e+7;
    stress[3] = -1e+5;
    stress[4] = -2200.2;
    stress[5] = -8e+5;
    vars[static_cast<int>(MPL::Variable::stress)] = stress;
    auto const k_non_f =
        MPL::formEigenTensor<3>(k_model.value(vars, pos, t, dt));

    auto const k_non_f_expected = k0(t, pos)[0];

    ASSERT_LE(std::fabs(k_non_f_expected - k_non_f(0, 0)) / k_non_f_expected,
              1e-19)
        << "for expected untouched permeability" << k_non_f_expected
        << " and for computed untouched permeability." << k_non_f(0, 0);

    // Stress at the apex. No change in permeability.
    const double val =
        2.0 * c / std::tan(phi * boost::math::constants::degree<double>());
    stress[0] = 0.7 * val;
    stress[1] = 0.3 * val;
    stress[2] = 0.3 * val;
    stress[3] = 0.0;
    stress[4] = 0.0;
    stress[5] = 0.0;
    vars[static_cast<int>(MPL::Variable::stress)] = stress;
    auto const k_apex_f =
        MPL::formEigenTensor<3>(k_model.value(vars, pos, t, dt));

    double const k_apex_expected = 7.2849707418474819e-15;
    ASSERT_LE(std::fabs(k_apex_expected - k_apex_f(0, 0)) / k_apex_expected,
              5e-16)
        << "for expected untouched permeability" << k_non_f_expected
        << " and for computed untouched permeability." << k_apex_f(0, 0);
}
