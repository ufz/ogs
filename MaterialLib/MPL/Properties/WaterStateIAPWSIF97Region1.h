// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <spdlog/fmt/fmt.h>

#include <string_view>

#include "NumLib/Exceptions.h"

namespace MaterialPropertyLib::IAPWSIF97Region1
{
/// Upper bound of the pressure range of region 1, see
/// http://www.iapws.org/relguide/IF97-Rev.pdf p.5, section 4.
constexpr double maximum_pressure = 100e6;  ///< Pa

/// Upper bound of the temperature range of region 1. Above it the formulation
/// describes region 3 rather than the compressed liquid, see
/// http://www.iapws.org/relguide/IF97-Rev.pdf p.5, section 4.
constexpr double maximum_temperature = 623.15;  ///< K

/// Aborts the assembly if the state lies outside the region 1 range, that is
/// outside the compressed liquid the region 1 correlations describe. Like the
/// region 4 correlations, they are extrapolated beyond their range without
/// saying so, and a caller that reaches this check has already established
/// that no saturation state exists at this pressure, so there is nothing left
/// to fall back on.
///
/// The state follows the solution iterate, so an AssemblyException lets the
/// time stepper repeat the step with a smaller step size. A section that is
/// genuinely at such a state is not one this process can describe: it has no
/// region 2 or region 3 properties.
/// \param pressure     the pressure to check, in Pa.
/// \param temperature  the temperature to check, in K.
/// \param quantity     what is being evaluated, named in the message.
inline void checkStateInRange(double const pressure, double const temperature,
                              std::string_view const quantity)
{
    if ((pressure > maximum_pressure) || (temperature > maximum_temperature))
    {
        throw NumLib::AssemblyException(fmt::format(
            "Pressure {:g} Pa and temperature {:g} K are out of the range of "
            "the IAPWS-IF97 region 1 correlations, p <= {:g} Pa and T <= {:g} "
            "K, for {}.",
            pressure, temperature, maximum_pressure, maximum_temperature,
            quantity));
    }
}
}  // namespace MaterialPropertyLib::IAPWSIF97Region1
