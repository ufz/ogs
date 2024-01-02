/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on Feb 8, 2023, 3:05 PM
 */

#pragma once

#include <array>
#include <cmath>

namespace MaterialPropertyLib::IAPWSIF97Region4
{

/// The saturation-temperature equation function in region 4, from
/// "The International Association for the Properties of Water and Steam"
/// (see http://www.iapws.org/relguide/IF97-Rev.pdf) p.35, section 8.2.
inline double waterSaturationTemperature(double const pressure)
{
    static constexpr std::array n = {0.11670521452767e4,  -0.72421316703206e6,
                                     -0.17073846940092e2, 0.12020824702470e5,
                                     -0.32325550322333e7, 0.14915108613530e2,
                                     -0.48232657361591e4, 0.40511340542057e6,
                                     -0.23855557567849,   0.65017534844798e3};

    static constexpr double p_c = 1e6;

    double const beta2 = std::sqrt(pressure / p_c);
    double const beta = std::sqrt(beta2);

    double const E = beta2 + n[2] * beta + n[5];
    double const F = n[0] * beta2 + n[3] * beta + n[6];
    double const G = n[1] * beta2 + n[4] * beta + n[7];

    double const D = 2 * G / (-F - std::sqrt(F * F - 4 * E * G));

    double const n10pD = n[9] + D;

    return (n10pD - std::sqrt(n10pD * n10pD - 4 * (n[8] + n[9] * D))) / 2;
}
}  // namespace MaterialPropertyLib::IAPWSIF97Region4
