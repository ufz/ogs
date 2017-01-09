/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DensityDubinin.h"
#include "MaterialLib/PhysicalConstant.h"
#include "DensityCook.h"
#include "Adsorption.h"

namespace
{

// NaX_Dubinin_polyfrac_CC.pickle
// date extracted 2015-06-23 16:47:50 file mtime 2015-06-23 16:47:23
const double c[] = {
    0.3635538371322433,        /* a0 */
    -0.0014521033261199435,    /* a1 */
    -0.0007855160157616825,    /* a2 */
    4.385666000850253e-08,    /* a3 */
    5.567776459188524e-07,    /* a4 */
    6.026002134230559e-10,    /* a5 */
    -1.0477401124006098e-10    /* a6 */
};

}

namespace Adsorption
{

double DensityDubinin::getAdsorbateDensity(const double T_Ads) const
{
    const double Tb = 373.1;

    if (T_Ads < Tb) {
        return rhoWaterDean(T_Ads);
    } else {
        const double Tc = 647.3; // K
        const double pc = 221.2e5; // Pa
        // boiling point density
        const double rhob = rhoWaterDean(Tb);
        // state values
        const double R = MaterialLib::PhysicalConstant::IdealGasConstant;
        const double M = MaterialLib::PhysicalConstant::MolarMass::Water;
        const double b = R * Tc/(8. * pc); // m^3/mol
        const double rhom = M/b; // kg/m^3
        const double rho = rhob - (rhob-rhom)/(Tc-Tb)*(T_Ads-Tb);
        return rho;
    }
}

//Thermal expansivity model for water found in the works of Hauer
double DensityDubinin::getAlphaT(const double T_Ads) const
{
    const double Tb = 373.1;
    if (T_Ads <= Tb) {
        return alphaTWaterDean(T_Ads);
    } else {
        // critical T and p
        const double Tc = 647.3; // K
        // const double rhoc = 322.; // kg/m^3
        const double pc = 221.2e5; // Pa
        // boiling point density
        const double rhob = rhoWaterDean(Tb);
        // state values
        const double R = MaterialLib::PhysicalConstant::IdealGasConstant;
        const double M = MaterialLib::PhysicalConstant::MolarMass::Water;
        const double b = R * Tc/(8. * pc); // m^3/mol
        const double rhom = M/(b); // kg/m^3
        const double rho = rhob - (rhob-rhom)/(Tc-Tb)*(T_Ads-Tb);
        return ((rhob-rhom)/(Tc-Tb)*1./rho);
    }
}

// Characteristic curve. Return W (A)
double DensityDubinin::characteristicCurve(const double A) const
{
    double W = curvePolyfrac(c, A); // cm^3/g

    if (W < 0.0) {
        W = 0.0; // TODO [CL] debug output
    }

    return W/1.e3; // m^3/kg
}

double DensityDubinin::dCharacteristicCurve(const double A) const
{
    return dCurvePolyfrac(c, A);
}

}
