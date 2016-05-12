/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
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

double DensityConst::get_adsorbate_density(const double /*T_Ads*/) const
{
    return rho_water_Hauer(150.0+273.15);
}


//Thermal expansivity model for water found in the works of Hauer
double DensityConst::get_alphaT(const double /*T_Ads*/) const
{
    return 0.0;
}


//Characteristic curve. Return W (A)
double DensityConst::characteristic_curve(const double A) const
{
    double W = curve_polyfrac(c, A); //cm^3/g

    if (W < 0.0) {
        W = 0.0; // TODO [CL] debug output
    }

    return W/1.e3; //m^3/kg
}

double DensityConst::d_characteristic_curve(const double A) const
{
    return d_curve_polyfrac(c, A);
}

}
