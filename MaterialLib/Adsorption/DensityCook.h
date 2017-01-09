/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATERIALSLIB_ADSORPTION_DENSITYCOOK_H
#define MATERIALSLIB_ADSORPTION_DENSITYCOOK_H

#include "Adsorption.h"

namespace Adsorption
{

class DensityCook : public AdsorptionReaction
{
public:
    double getAdsorbateDensity(const double T_Ads) const;
    double getAlphaT(const double T_Ads) const;
    double characteristicCurve(const double A) const;
    double dCharacteristicCurve(const double A) const;
};

inline double rhoWaterDean(const double T_Ads)
{
    const double Tcel = T_Ads - 273.15;
    const double b[] = { 999.9,2.03E-02,-6.16E-03,2.26E-05,-4.68E-08 };
    if (Tcel <= 100.) {
        return b[0] + Tcel * (b[1] + Tcel * (b[2] + Tcel * (b[3] + Tcel * b[4]) ) );
    }
    else {
        const double rho_100 = b[0] + b[1]*1.e2 + b[2]*1.e4 + b[3]*1.e6 + b[4]*1.e8;
        const double aT_100  = -1./rho_100 * (b[1] + 2.*b[2]*1.e2 + 3.*b[3]*1.e4 + 4.*b[4]*1.e6);
        return rho_100 * (1. - aT_100*(Tcel-100.));
    }
}

inline double alphaTWaterDean(const double T_Ads)
{
    const double Tcel = T_Ads - 273.15;
    const double b[] = { 999.9,2.03E-02,-6.16E-03,2.26E-05,-4.68E-08 };
    if (Tcel <= 100.) {
        const double r = b[0] + Tcel * (b[1] + Tcel * (b[2] + Tcel * (b[3] + Tcel * b[4]) ) );
        return -1.0/r * ( b[1] + Tcel * (2.0*b[2] + Tcel * (3.0*b[3] + Tcel * 4.0*b[4]) ) );
    }
    else {
        const double rho_100 = b[0] + b[1]*1.e2 + b[2]*1.e4 + b[3]*1.e6 + b[4]*1.e8;
        const double aT_100  = -1./rho_100 * (b[1] + 2.*b[2]*1.e2 + 3.*b[3]*1.e4 + 4.*b[4]*1.e6);
        return aT_100 / (1. - aT_100*(Tcel-100.));
    }
}

}
#endif // MATERIALSLIB_ADSORPTION_DENSITYCOOK_H
