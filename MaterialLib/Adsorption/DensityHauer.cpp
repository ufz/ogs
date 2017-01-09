/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DensityHauer.h"

namespace
{

// NaX_Hauer_polyfrac_CC.pickle
// date extracted 2015-06-23 15:38:35 file mtime 2015-06-23 15:19:19
const double c[] = {
    0.36490158988356747,    /* a0 */
    -0.0013723270478333963,    /* a1 */
    -0.0007655780628099964,    /* a2 */
    -3.353324854315774e-08,    /* a3 */
    5.424357157710913e-07,    /* a4 */
    6.613430586648678e-10,    /* a5 */
    -1.0300151379421499e-10    /* a6 */
};

}

namespace Adsorption
{

double DensityHauer::getAdsorbateDensity(const double T_Ads) const
{
    return rhoWaterHauer(T_Ads);
}

// Thermal expansivity model for water found in the works of Hauer
double DensityHauer::getAlphaT(const double T_Ads) const
{
    // data like in python script
    const double T0 = 283.15, alpha0 = 3.781e-4; //K; 1/K

    return alpha0/(1. - alpha0 * (T_Ads-T0)); //in 1/K
}

// Characteristic curve. Return W (A)
double DensityHauer::characteristicCurve(const double A) const
{
    double W = curvePolyfrac(c, A); // cm^3/g

    if (W < 0.0) {
        W = 0.0; // TODO [CL] debug output
    }

    return W/1.e3; // m^3/kg
}

double DensityHauer::dCharacteristicCurve(const double A) const
{
    return dCurvePolyfrac(c, A);
}

}
