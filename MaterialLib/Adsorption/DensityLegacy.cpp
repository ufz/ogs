/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DensityLegacy.h"

namespace
{

//parameters from least squares fit (experimental data)
const double c[] = { 0.34102920966608297,
                     -0.0013106032830951296,
                     -0.00060754147575378876,
                     3.7843404172683339e-07,
                     4.0107503869519016e-07,
                     3.1274595098338057e-10,
                     -7.610441241719489e-11
                   };

}

namespace Adsorption
{

double DensityLegacy::getAdsorbateDensity(const double T_Ads) const
{
    //set reference state for adsorbate EOS in Hauer
    const double T0 = 293.15, rho0 = 998.084, alpha0 = 2.06508e-4; // K; kg/m^3; 1/K

    return (rho0 * (1. - alpha0 * (T_Ads-T0))); // in kg/m^3
}

// Thermal expansivity model for water found in the works of Hauer
double DensityLegacy::getAlphaT(const double T_Ads) const
{
    //set reference state for adsorbate EOS in Hauer
    const double T0 = 293.15, alpha0 = 2.06508e-4; // K; 1/K

    return (alpha0/(1. - alpha0 * (T_Ads-T0))); // in 1/K
}

// Characteristic curve. Return W (A)
double DensityLegacy::characteristicCurve(const double A) const
{
    double W = curvePolyfrac(c, A); // cm^3/g

    if (W < 0.0) {
        W = 0.0; // TODO [CL] debug output
    }

    return W/1.e3; // m^3/kg
}

double DensityLegacy::dCharacteristicCurve(const double A) const
{
    return dCurvePolyfrac(c, A);
}

}
