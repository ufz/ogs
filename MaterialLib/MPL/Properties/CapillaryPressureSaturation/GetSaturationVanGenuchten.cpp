/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on May 5, 2020, 10:11 AM
 */

#include "GetSaturationVanGenuchten.h"

#include <cmath>

namespace MaterialPropertyLib
{
double getSaturationVanGenuchten(double const p_c, double const p_b,
                                 double const S_L_res, double const S_L_max,
                                 double const m)
{
    double const p = p_c / p_b;
    double const n = 1. / (1. - m);
    double const p_to_n = std::pow(p, n);

    double const S_eff = std::pow(p_to_n + 1., -m);
    return S_eff * S_L_max - S_eff * S_L_res + S_L_res;
}

}  // namespace MaterialPropertyLib
