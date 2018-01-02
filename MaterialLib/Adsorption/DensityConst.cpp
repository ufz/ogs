/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DensityConst.h"
#include "DensityHauer.h"

namespace
{

// NaX_Constant_polyfrac_CC.pickle
// date extracted 2015-06-23 15:38:35 file mtime 2015-06-23 15:20:05
const double c[] = {
    0.3824098506898007,            /* a0 */
    -0.001316857559708455,        /* a1 */
    -0.0007935756090263691,        /* a2 */
    -1.1600036977157845e-07,    /* a3 */
    5.610354459181838e-07,        /* a4 */
    7.113664938298873e-10,        /* a5 */
    -1.0668790477629686e-10        /* a6 */
};

}

namespace Adsorption
{

double DensityConst::getAdsorbateDensity(const double /*T_Ads*/) const
{
    return rhoWaterHauer(150.0+273.15);
}

double DensityConst::getAlphaT(const double /*T_Ads*/) const
{
    return 0.0;
}

// Characteristic curve. Return W (A)
double DensityConst::characteristicCurve(const double A) const
{
    double W = curvePolyfrac(c, A); //cm^3/g

    if (W < 0.0) {
        W = 0.0; // TODO [CL] debug output
    }

    return W/1.e3; // m^3/kg
}

double DensityConst::dCharacteristicCurve(const double A) const
{
    return dCurvePolyfrac(c, A);
}

}
