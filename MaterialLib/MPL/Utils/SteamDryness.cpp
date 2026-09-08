// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "SteamDryness.h"

#include "BaseLib/Logging.h"

namespace MaterialPropertyLib
{
double steamDryness(double const enthalpy, double const h_sat_liquid,
                    double const h_sat_vapour)
{
    double const dryness =
        (enthalpy - h_sat_liquid) / (h_sat_vapour - h_sat_liquid);

    if (dryness > 1)
    {
        WARN(
            "Specific enthalpy {:g} J/kg exceeds the saturation enthalpy of "
            "the vapour phase, {:g} J/kg. The steam is superheated, which the "
            "saturated mixture cannot represent: the section is evaluated at "
            "the saturation temperature and with the saturation vapour "
            "density.",
            enthalpy, h_sat_vapour);
        return 1.;
    }
    if (dryness < 0)
    {
        return 0.;
    }
    return dryness;
}
}  // namespace MaterialPropertyLib
