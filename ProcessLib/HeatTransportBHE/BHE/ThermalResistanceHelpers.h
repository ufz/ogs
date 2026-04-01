// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cmath>

#include "BaseLib/Error.h"

namespace ProcessLib
{
namespace HeatTransportBHE
{
namespace BHE
{
/// Grout-soil thermal resistance: R_gs = (1 - chi) * R_g.
inline double computeRgs(double const chi, double const R_g)
{
    return (1 - chi) * R_g;
}

/// Inter-grout thermal resistance.
/// Eq. (in Diersch_2011_CG).
inline double computeRgg(double const chi, double const R_gs, double const R_ar,
                         double const R_g)
{
    double const R_gg = 2.0 * R_gs * (R_ar - 2.0 * chi * R_g) /
                        (2.0 * R_gs - R_ar + 2.0 * chi * R_g);
    if (!std::isfinite(R_gg))
    {
        OGS_FATAL(
            "Error!!! Grout Thermal Resistance is an infinite number! The "
            "simulation will be stopped!");
    }
    return R_gg;
}
}  // namespace BHE
}  // namespace HeatTransportBHE
}  // namespace ProcessLib
