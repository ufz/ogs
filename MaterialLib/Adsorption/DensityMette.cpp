/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DensityMette.h"
#include "DensityCook.h"

namespace
{
// NaX_Mette_polyfrac_CC.pickle
// date extracted 2015-06-23 15:38:35 file mtime 2015-06-23 15:19:26
const double c[] = {
    0.36340572890087813,    /* a0 */
    -0.0013449597038375108,    /* a1 */
    -0.0007581210111121073,    /* a2 */
    -7.331279615575401e-08,    /* a3 */
    5.365656973806218e-07,    /* a4 */
    6.854673678427112e-10,    /* a5 */
    -1.0197050219481966e-10    /* a6 */
};
}

namespace Adsorption
{

double DensityMette::getAdsorbateDensity(const double T_Ads) const
{
    const double T0 = 293.15;
    const double rho0 = rhoWaterDean(T0);
    const double alpha20 = alphaTWaterDean(T0);
    return rho0 / (1. + alpha20*(T_Ads-T0));
}


// Thermal expansivity model for water found in the works of Hauer
double DensityMette::getAlphaT(const double T_Ads) const
{
    const double T0 = 293.15;
    const double alpha20 = alphaTWaterDean(T0);
    return alpha20 / (1. + alpha20 * (T_Ads-T0));
}


// Characteristic curve. Return W (A)
double DensityMette::characteristicCurve(const double A) const
{
    double W = curvePolyfrac(c, A); // cm^3/g

    if (W < 0.0) {
        W = 0.0; // TODO [CL] debug output
    }

    return W/1.e3; // m^3/kg
}

double DensityMette::dCharacteristicCurve(const double A) const
{
    return dCurvePolyfrac(c, A);
}

}
