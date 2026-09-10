// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <spdlog/fmt/fmt.h>

#include <array>
#include <cmath>
#include <string_view>

#include "MaterialLib/PhysicalConstant.h"
#include "NumLib/Exceptions.h"

namespace MaterialPropertyLib::IAPWSIF97Region4
{
/// Lower bound of the pressure range of the region 4 saturation line, the
/// saturation pressure at 273.15 K. It is fixed by the formulation itself,
/// see http://www.iapws.org/relguide/IF97-Rev.pdf p.33, section 8.1, and is
/// neither the triple point pressure of ordinary water nor to be updated to
/// it. The upper bound of the range is the critical pressure.
constexpr double minimum_saturation_pressure = 611.213;  ///< Pa

/// Aborts the assembly if \c pressure lies outside the pressure range of the
/// region 4 saturation line. Outside it the correlations below are
/// extrapolated and no longer describe water, and they do not say so
/// themselves: they return plausible looking values well outside the range,
/// 998 kg/m^3 at half the lower bound and 618 kg/m^3 at twice the critical
/// pressure, where no saturation state exists at all. Nothing downstream can
/// catch what this check lets through.
///
/// The pressure is a solution iterate, so leaving the range is a diverging
/// step rather than a broken input, and an AssemblyException lets the time
/// stepper repeat the step with a smaller step size, as for every other
/// inadmissible state the closure reports.
/// \param pressure  the pressure to check, in Pa.
/// \param quantity  the quantity being evaluated, named in the message.
inline void checkPressureInRange(double const pressure,
                                 std::string_view const quantity)
{
    constexpr double maximum_saturation_pressure =
        MaterialLib::PhysicalConstant::CriticalPoint::PressureWater;

    if ((pressure < minimum_saturation_pressure) ||
        (pressure > maximum_saturation_pressure))
    {
        throw NumLib::AssemblyException(fmt::format(
            "Pressure {:g} Pa is out of the range [{:g}, {:g}] Pa for {}.",
            pressure, minimum_saturation_pressure, maximum_saturation_pressure,
            quantity));
    }
}

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
