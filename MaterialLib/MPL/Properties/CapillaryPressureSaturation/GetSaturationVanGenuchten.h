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

#pragma once

namespace MaterialPropertyLib
{
/**
 * \fn getSaturationVanGenuchten(double const p_c, double const p_b,
                                 double const S_L_res, double const S_L_max,
                                 double const m)
 * \brief A common function of class SaturationVanGenuchten and class
 * CapillaryPressureVanGenuchten.
 * It is used to compute saturation via capillary pressure.
 * \sa MaterialPropertyLib::SaturationVanGenuchten
 */
double getSaturationVanGenuchten(double const p_c, double const p_b,
                                 double const S_L_res, double const S_L_max,
                                 double const m);

}  // namespace MaterialPropertyLib
