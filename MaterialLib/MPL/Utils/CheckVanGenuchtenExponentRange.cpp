/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on April 20, 2020, 10:47 AM
 */

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
