/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 19, 2021, 4:28 PM
 */

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <limits>

#include "MaterialLib/MPL/Properties/Enthalpy/CreateGeneralWaterVapourLatentHeat.h"
#include "MaterialLib/PhysicalConstant.h"
#include "TestMPL.h"

TEST(MaterialPropertyLib, GeneralWaterVapourLatentHeat)
{
    char const xml[] =
        "<property>"
        "   <name>latent_heat</name>"
        "   <type>GeneralWaterVapourLatentHeat</type>"
        "</property>";

    std::unique_ptr<MaterialPropertyLib::Property> const property_ptr =
        Tests::createTestProperty(
            xml, MaterialPropertyLib::createGeneralWaterVapourLatentHeat);
    MaterialPropertyLib::Property const& property = *property_ptr;

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    double const T_0 = MaterialLib::PhysicalConstant::CelsiusZeroInKelvin;
    double const T_max =
        400.0 + MaterialLib::PhysicalConstant::CelsiusZeroInKelvin;

    int const n = 100;
    double const dT = (T_max - T_0) / static_cast<double>(n);

    std::array const L_wv_expected = {
        2.494034e+06, 2.484795e+06, 2.475606e+06, 2.466455e+06, 2.457332e+06,
        2.448226e+06, 2.439128e+06, 2.430027e+06, 2.420916e+06, 2.411784e+06,
        2.402625e+06, 2.393429e+06, 2.384189e+06, 2.374897e+06, 2.365546e+06,
        2.356129e+06, 2.346640e+06, 2.337072e+06, 2.327418e+06, 2.317673e+06,
        2.307830e+06, 2.297884e+06, 2.287829e+06, 2.277659e+06, 2.267369e+06,
        2.256954e+06, 2.246408e+06, 2.235726e+06, 2.224902e+06, 2.213932e+06,
        2.202809e+06, 2.191530e+06, 2.180087e+06, 2.168476e+06, 2.156691e+06,
        2.144726e+06, 2.132576e+06, 2.120235e+06, 2.107695e+06, 2.094951e+06,
        2.081996e+06, 2.068823e+06, 2.055425e+06, 2.041793e+06, 2.027921e+06,
        2.013800e+06, 1.999421e+06, 1.984775e+06, 1.969852e+06, 1.954643e+06,
        1.939135e+06, 1.923319e+06, 1.907182e+06, 1.890711e+06, 1.873894e+06,
        1.856715e+06, 1.839159e+06, 1.821211e+06, 1.802854e+06, 1.784069e+06,
        1.764837e+06, 1.745138e+06, 1.724948e+06, 1.704246e+06, 1.683004e+06,
        1.661196e+06, 1.638792e+06, 1.615760e+06, 1.592066e+06, 1.567672e+06,
        1.542537e+06, 1.516615e+06, 1.489858e+06, 1.462210e+06, 1.433610e+06,
        1.403989e+06, 1.373270e+06, 1.341365e+06, 1.308172e+06, 1.273575e+06,
        1.237434e+06, 1.199586e+06, 1.159834e+06, 1.117934e+06, 1.073580e+06,
        1.026382e+06, 9.758218e+05, 9.211899e+05, 8.614708e+05, 7.951216e+05,
        7.195987e+05, 6.301645e+05, 5.159389e+05, 3.360458e+05, 0.0,
        0.0,          0.0,          0.0,          0.0,          0.0};

    for (std::size_t i = 0; i < n; ++i)
    {
        double const T = T_0 + i * dT;
        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::temperature)] = T;

        double const L_wv =
            property.template value<double>(variable_array, pos, t, dt);

        double const relative_reference =
            (L_wv_expected[i] > 0) ? L_wv_expected[i] : 1.0;
        ASSERT_LE(std::fabs(L_wv_expected[i] - L_wv) / relative_reference, 1e-6)
            << "for expected water vapour latent heat " << L_wv_expected[i]
            << " and for computed water vapour latent heat " << L_wv;

        double const dT_i = 1.0e-6;

        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::temperature)] = T - dT_i;
        double const L_wv0 =
            property.template value<double>(variable_array, pos, t, dt);

        variable_array[static_cast<int>(
            MaterialPropertyLib::Variable::temperature)] = T + dT_i;
        double const L_wv1 =
            property.template value<double>(variable_array, pos, t, dt);

        double const approximated_dLw_dT = 0.5 * (L_wv1 - L_wv0) / dT_i;

        double const analytic_dLw_dT = property.template dValue<double>(
            variable_array, MaterialPropertyLib::Variable::temperature, pos, t,
            dt);

        double const drelative_reference_dT =
            (std::fabs(analytic_dLw_dT) > 0) ? analytic_dLw_dT : 1.0;
        ASSERT_LE(std::fabs((approximated_dLw_dT - analytic_dLw_dT) /
                            drelative_reference_dT),
                  1e-6)
            << "for the approximated derivative of water vapour latent heat "
               "with "
               "respect to temperature "
            << approximated_dLw_dT
            << " and for the analytical computed derivative of water vapour "
               "latent heat "
               "with respect to temperature."
            << analytic_dLw_dT;
    }
}
