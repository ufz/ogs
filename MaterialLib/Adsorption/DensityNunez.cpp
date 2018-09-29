/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DensityNunez.h"

namespace
{

// NaX_Nunez_polyfrac_CC.pickle
// date extracted 2015-06-23 15:38:35 file mtime 2015-06-23 15:19:34
const double c[] = {
    0.3631900485031771,        /* a0 */
    -0.0014242280940080726,    /* a1 */
    -0.0007751726942386291,    /* a2 */
    2.1775655036811842e-08,    /* a3 */
    5.488166913667265e-07,    /* a4 */
    6.204064716725214e-10,    /* a5 */
    -1.0345385018952998e-10    /* a6 */
};

}

namespace Adsorption
{

double DensityNunez::getAdsorbateDensity(const double T_Ads) const
{
    // TODO admissable T range: 273.16 K <= T_Ads <= 633.15 K
    const double a[] = { 1.0644e3,-8.01905,1.445348e-2,-4.19589e-6,-4.5294e-9 };
    const double b[] = { -8.039e-3,1.8698e-5,-2.3015e-8,2.3809e-11,-1.388e-14 };
    const double u = a[0] + T_Ads * (a[1] + T_Ads * (a[2] + T_Ads * (a[3] + T_Ads * a[4]) ) );
    const double v = 1.0 + T_Ads * (b[0] + T_Ads * (b[1] + T_Ads * (b[2] + T_Ads * (b[3] + T_Ads * b[4]) ) ) );
    return u/v;
}


// Thermal expansivity model for water found in the works of Hauer
double DensityNunez::getAlphaT(const double T_Ads) const
{
    // TODO admissable T range: 273.16 K <= T_Ads <= 633.15 K
    const double a[] = { 1.0644e3,-8.01905,1.445348e-2,-4.19589e-6,-4.5294e-9 };
    const double b[] = { -8.039e-3,1.8698e-5,-2.3015e-8,2.3809e-11,-1.388e-14 };
    const double u = a[0] + T_Ads * (a[1] + T_Ads * (a[2] + T_Ads * (a[3] + T_Ads * a[4]) ) );
    const double v = 1.0 + T_Ads * (b[0] + T_Ads * (b[1] + T_Ads * (b[2] + T_Ads * (b[3] + T_Ads * b[4]) ) ) );
    const double du = a[1] + T_Ads * (2.0*a[2] + T_Ads * (3.0*a[3] + T_Ads * 4.0*a[4]) );
    const double dv = b[0] + T_Ads * (2.0*b[1] + T_Ads * (3.0*b[2] + T_Ads * (4.0*b[3] + T_Ads * 5.0*b[4]) ) );
    return dv/v - du/u;
}


// Characteristic curve. Return W (A)
double DensityNunez::characteristicCurve(const double A) const
{
    double W = curvePolyfrac(c, A); // cm^3/g

    if (W < 0.0) {
        W = 0.0; // TODO [CL] debug output
    }

    return W/1.e3; // m^3/kg
}

double DensityNunez::dCharacteristicCurve(const double A) const
{
    return dCurvePolyfrac(c, A);
}

}
