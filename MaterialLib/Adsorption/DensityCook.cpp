/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DensityCook.h"

namespace
{

// NaX_Dean_polyfrac_CC.pickle
// date extracted 2015-06-23 15:38:35 file mtime 2015-06-23 15:19:42
const double c[] = {
    0.3632627555646154,        /* a0 */
    -0.0014090624975800715,    /* a1 */
    -0.0007717609035743321,    /* a2 */
    5.03634836561135e-09,    /* a3 */
    5.478509959282738e-07,    /* a4 */
    6.36458510620815e-10,    /* a5 */
    -1.037977321231462e-10    /* a6 */
};

}

namespace Adsorption
{

double DensityCook::getAdsorbateDensity(const double T_Ads) const
{
    return rhoWaterDean(T_Ads);
}

double DensityCook::getAlphaT(const double T_Ads) const
{
    return alphaTWaterDean(T_Ads);
}

// Characteristic curve. Return W (A)
double DensityCook::characteristicCurve(const double A) const
{
    double W = curvePolyfrac(c, A); //cm^3/g

    if (W < 0.0) {
        W = 0.0; // TODO [CL] debug output
    }

    return W/1.e3; //m^3/kg
}

double DensityCook::dCharacteristicCurve(const double A) const
{
    return dCurvePolyfrac(c, A);
}

}
