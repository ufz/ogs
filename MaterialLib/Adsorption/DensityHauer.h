/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Adsorption.h"
#include "DensityCook.h"

namespace Adsorption
{

class DensityHauer : public AdsorptionReaction
{
public:
    double getAdsorbateDensity(const double T_Ads) const override;
    double getAlphaT(const double T_Ads) const override;
    double characteristicCurve(const double A) const override;
    double dCharacteristicCurve(const double A) const override;
};

inline double rhoWaterHauer(const double T_Ads)
{
    // data like in python script
    const double T0 = 283.15, rho0 = rhoWaterDean(T0), alpha0 = 3.781e-4; // K; kg/m^3; 1/K

    return rho0 * (1. - alpha0 * (T_Ads-T0)); // in kg/m^3
}

}  // namespace Adsorption
