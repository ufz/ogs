// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CheckVanGenuchtenExponentRange.h"

#include "BaseLib/Error.h"

namespace MaterialPropertyLib
{
void checkVanGenuchtenExponentRange(const double m)
{
    if (m <= 0 || m >= 1)
    {
        OGS_FATAL(
            "The exponent value m = {:e} of van Genuchten saturation model, is "
            "out of its range of(0, 1) ",
            m);
    }
}
}  // namespace MaterialPropertyLib
